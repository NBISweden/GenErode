##########################################################################
##################### WORKFLOW VARIABLES AND CODE ########################
##########################################################################

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate
import os
import pandas as pd

min_version("5.19.0")

configfile: "config/config.yaml"

report: "../report/workflow.rst"

# reference assembly variables
REF_DIR = os.path.dirname(config["ref_path"])
REF_FASTA = os.path.basename(config["ref_path"])
REF_NAME, REF_EXT = os.path.splitext(REF_FASTA)


### Global wildcard contraints that apply to all rules
wildcard_constraints:
    sample=r"[A-Za-z0-9]+",
    DP=r"[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?",  # avoid extension with "rm.autos" when running mlRho
    fmiss=r"[0-1].[0-9]+", # avoid extension with ".autos" when filtering the merged BCF file
    CpG_method=r"CpG_[vcfre]{3,6}",
    chunk=r"chunk[0-9]+",


### Code to generate lists and dictionaries from the metadata tables
def check_metadata_file(dataframe):
    if "samplename" in dataframe.columns and "library_id" in dataframe.columns and "lane" in dataframe.columns:
        # concatenate the sample name, library id and lane number with "_" to create a unique identifier for each fastq file
        dataframe["samplename_index_lane"] = dataframe["samplename"] + "_" + dataframe["library_id"] + "_" + dataframe["lane"]
    # raise an error if a sample exists with both fastq and bam files
    if "path_to_processed_bam_file" in dataframe.columns and "path_to_R1_fastq_file" in dataframe.columns and "path_to_R2_fastq_file" in dataframe.columns:
        # check if any sample has both fastq and bam files
        fastq_samples = dataframe["samplename"][dataframe["path_to_R1_fastq_file"].notnull()].unique()
        bam_samples = dataframe["samplename"][dataframe["path_to_processed_bam_file"].notnull()].unique()
        # check if any sample has both fastq and bam files
        if set(fastq_samples) & set(bam_samples):
            # raise an error if any sample has both fastq and bam files
            raise WorkflowError("Both fastq files and bam files are provided for some samples. Please check your metadata file.")
    return dataframe

# Create sample lists
# sample list for each fastq-file (["samplename_index_lane"])
def samplename_index_lane_func(dataframe):
    if "samplename" in dataframe.columns and "library_id" in dataframe.columns and "lane" in dataframe.columns:
        # concatenate the sample name, library id and lane number with "_" to create a unique identifier for each fastq file
        dataframe["samplename_index_lane"] = dataframe["samplename"] + "_" + dataframe["library_id"] + "_" + dataframe["lane"]
        if dataframe["samplename_index_lane"].duplicated().any():
            raise WorkflowError("Samples found with duplicate library ID and lane number. Please check your metadata file.")
        # convert the first column into a list
        return list(dataframe["samplename_index_lane"].drop_duplicates())
    else:
        return []


# sample list for merging of files per lane (["samplename_index"])
def samplename_index_func(dataframe):
    if "samplename" in dataframe.columns and "library_id" in dataframe.columns and "lane" in dataframe.columns:
        dataframe["samplename_index"] = dataframe["samplename"] + "_" + dataframe["library_id"]
        return list(dataframe["samplename_index"].drop_duplicates())
    else:
        return []


# sample list (["samplename"])
def samplename_func(dataframe):
    return list(dataframe["samplename"].drop_duplicates())


# Functions to create symbolic links to fastq files, to look up read group information and for merging of bam files
# symbolic links dictionary
def fastq_symlinks_dict_func(dataframe):
    fastq_symlinks_dict = {}
    if "path_to_R1_fastq_file" in dataframe.columns and "path_to_R2_fastq_file" in dataframe.columns:
        for index, row in dataframe.iterrows():
            fastq_symlinks_dict[row["samplename_index_lane"]] = {"R1": os.path.abspath(row["path_to_R1_fastq_file"]), "R2": os.path.abspath(row["path_to_R2_fastq_file"])}
    return fastq_symlinks_dict


# read group dictionary
def rg_dict_func(dataframe):
    rg_dict = {}
    if "readgroup_id" in dataframe.columns and "readgroup_platform" in dataframe.columns and "library_id" in dataframe.columns:
        for index, row in dataframe.iterrows():
            rg_dict[row["samplename_index_lane"]] = {"ID": row["readgroup_id"], "SM": row["samplename"], "PL": row["readgroup_platform"], "LB": row["library_id"]}
    return rg_dict


# dictionary for bam file merging per lane
def sampleidxln_dict_func(dataframe):
    sampleidxln_dict = {}
    if "samplename" in dataframe.columns and "library_id" in dataframe.columns and "lane" in dataframe.columns:
        for index, row in dataframe.iterrows():
            smid = row["samplename"] + "_" + row["library_id"]  # take the sample name plus index
            smidln = row["samplename"] + "_" + row["library_id"] + "_" + row["lane"]
            if (smid in sampleidxln_dict):  # if "sample_index" is already in the dictionary
                if (smidln not in sampleidxln_dict[smid]):  # if "sample_index_lane" is not in the list for "sample_index"
                    sampleidxln_dict[smid].append(smidln)  # add "sample_index_lane" for "sample_index"
            else:  # if "sample_index" is not yet in the dictionary
                sampleidxln_dict[smid] = [smidln]  # add "sample_index_lane" for "sample_index"
    return sampleidxln_dict


# dictionary for bam file merging per library ID
def sampleidx_dict_func(dataframe):
    sampleidx_dict = {}
    if "samplename" in dataframe.columns and "library_id" in dataframe.columns and "lane" in dataframe.columns:
        for index, row in dataframe.iterrows():
            sm = row["samplename"]  # take the sample name from the first column of each line
            smid = row["samplename"] + "_" + row["library_id"]  # take the sample name plus index
            if sm in sampleidx_dict:  # if "sample" is already in the dictionary
                if (smid not in sampleidx_dict[sm]):  # if "sample_index" is not in the list for "sample"
                    sampleidx_dict[sm].append(smid)  # add "sample_index" for "sample"
            else:  # if "sample" is not yet in the dictionary
                sampleidx_dict[sm] = [smid]  # add "sample_index" for "sample"
    return sampleidx_dict


# dictionary for bam file merging per sample (mitochondrial genomes)
def sample_dict_func(dataframe):
    sample_dict = {}
    if "samplename" in dataframe.columns and "library_id" in dataframe.columns and "lane" in dataframe.columns:
        for index, row in dataframe.iterrows():
            sm = row["samplename"]  # take the sample name from the first column of each line
            smidln = row["samplename"] + "_" + row["library_id"] + "_" + row["lane"]
            if sm in sample_dict:  # if "sample" is already in the dictionary
                if (smidln not in sample_dict[sm]):  # if "sample_index_lane" is not in the list for "sample"
                    sample_dict[sm].append(smidln)  # add "sample_index_lane" for "sample"
            else:  # if "sample" is not yet in the dictionary
                sample_dict[sm] = [smidln]  # add "sample_index_lane" for "sample"
    return sample_dict

# Functions to collect user-provided bam files
# symbolic links dictionaries
def user_bam_symlinks_dict_func(dataframe):
    if "path_to_processed_bam_file" in dataframe.columns:
        user_bam_symlinks_dict = {}
        for index, row in dataframe.iterrows():
            user_bam_symlinks_dict[row["samplename"]] = {"bam": os.path.abspath(row["path_to_processed_bam_file"])}
        return user_bam_symlinks_dict

# lists of samples with user-provided bams and pipeline-generated bams
def user_bam_samples_func(dataframe):
    if "path_to_processed_bam_file" in dataframe.columns:
        user_bam_samples = list(dataframe["samplename"][dataframe["path_to_processed_bam_file"].notnull()])
        if len(user_bam_samples) != len(set(user_bam_samples)):
            raise WorkflowError("Samples found with duplicate user-provided bam files. Please check your metadata file.")
    else:
        user_bam_samples = []
    return user_bam_samples

def pipeline_bam_samples_func(dataframe):
    if "path_to_processed_bam_file" in dataframe.columns:
        pipeline_bam_samples = list(dataframe["samplename"][dataframe["path_to_processed_bam_file"].isnull()].unique())
    else:
        pipeline_bam_samples = list(dataframe["samplename"].drop_duplicates())
    return pipeline_bam_samples

# Apply the functions to metadata tables for historical and modern samples
if os.path.exists(config["historical_samples"]):
    historical_df = pd.read_csv(config["historical_samples"], sep=";|,| |\t", engine='python', dtype=str)  # read in the metadata as dataframe
    validate(historical_df, schema="../schemas/metadata.schema.yaml")  # validate metadata file format with JSON schema
    historical_df = check_metadata_file(historical_df)  # check metadata file format
    hist_sm = samplename_func(historical_df)
    hist_sm_idx = samplename_index_func(historical_df)
    hist_sm_idx_ln = samplename_index_lane_func(historical_df)
    hist_fastq_symlinks_dict = fastq_symlinks_dict_func(historical_df)
    hist_sample_dict = sample_dict_func(historical_df)
    hist_rg_dict = rg_dict_func(historical_df)
    hist_sampleidxln_dict = sampleidxln_dict_func(historical_df)
    hist_sampleidx_dict = sampleidx_dict_func(historical_df)
    hist_user_bam_symlinks_dict = user_bam_symlinks_dict_func(historical_df)
    hist_user_bam_sm = user_bam_samples_func(historical_df)
    hist_pipeline_bam_sm = pipeline_bam_samples_func(historical_df)
else:
    hist_sm = []
    hist_sm_idx = []
    hist_sm_idx_ln = []
    hist_fastq_symlinks_dict = {}
    hist_sample_dict = {}
    hist_rg_dict = {}
    hist_sampleidxln_dict = {}
    hist_sampleidx_dict = {}
    hist_user_bam_symlinks_dict = {}
    hist_user_bam_sm = []
    hist_pipeline_bam_sm = []


if os.path.exists(config["modern_samples"]):
    modern_df = pd.read_csv(config["modern_samples"], sep=";|,| |\t", engine='python', dtype=str)  # read in the metadata as dataframe
    validate(modern_df, schema="../schemas/metadata.schema.yaml") # validate metadata file format with JSON schema
    modern_df = check_metadata_file(modern_df)  # check metadata file format
    mod_sm = samplename_func(modern_df)
    mod_sm_idx = samplename_index_func(modern_df)
    mod_sm_idx_ln = samplename_index_lane_func(modern_df)
    mod_fastq_symlinks_dict = fastq_symlinks_dict_func(modern_df)
    mod_rg_dict = rg_dict_func(modern_df)
    mod_sampleidxln_dict = sampleidxln_dict_func(modern_df)
    mod_sampleidx_dict = sampleidx_dict_func(modern_df)
    mod_user_bam_symlinks_dict = user_bam_symlinks_dict_func(modern_df)
    mod_user_bam_sm = user_bam_samples_func(modern_df)
    mod_pipeline_bam_sm = pipeline_bam_samples_func(modern_df)
else:
    mod_sm = []
    mod_sm_idx = []
    mod_sm_idx_ln = []
    mod_fastq_symlinks_dict = {}
    mod_rg_dict = {}
    mod_sampleidxln_dict = {}
    mod_sampleidx_dict = {}
    mod_user_bam_symlinks_dict = {}
    mod_user_bam_sm = []
    mod_pipeline_bam_sm = []


### Parameters regarding optional steps

# mitochondrial genome fasta variables for mapping analysis
if config["map_historical_to_mitogenomes"]:
    # mitogenome fasta file variables for target species
    MITO_IN_DIR = os.path.dirname(config["species_mt_path"])
    MITO_IN_FASTA = os.path.basename(config["species_mt_path"])
    MITO_NAME, MITO_EXT = os.path.splitext(MITO_IN_FASTA)
    # 1.2 additional mitochondrial genomes to check how many reads map to them
    HUMAN_MT = "data/mitogenomes/human_NC_012920.fasta"
    CHICK_MT = "data/mitogenomes/chicken_NC_001323.fasta"
    COW_MT = "data/mitogenomes/cow_NC_006853.fasta"
    PIG_MT = "data/mitogenomes/pig_NC_000845.fasta"
    MOUSE_MT = "data/mitogenomes/mouse_NC_005089.fasta"
    MITO_DIR = "data/mitogenomes/"
    HUMAN_FASTA = os.path.basename(HUMAN_MT)
    HUMAN_NAME = os.path.splitext(HUMAN_FASTA)[0]
    CHICK_FASTA = os.path.basename(CHICK_MT)
    CHICK_NAME = os.path.splitext(CHICK_FASTA)[0]
    COW_FASTA = os.path.basename(COW_MT)
    COW_NAME = os.path.splitext(COW_FASTA)[0]
    PIG_FASTA = os.path.basename(PIG_MT)
    PIG_NAME = os.path.splitext(PIG_FASTA)[0]
    MOUSE_FASTA = os.path.basename(MOUSE_MT)
    MOUSE_NAME = os.path.splitext(MOUSE_FASTA)[0]


# Base quality rescaling, subsampling to common depth and CpG site removal
# Lists of historical samples
###
# rescaled (prior to subsampling and CpG filtering)
HIST_RESCALED_SAMPLES = list(
    set(hist_sm) & 
    set(config["historical_rescaled_samplenames"]))

# not rescaled (prior to subsampling and CpG filtering)
HIST_NOT_RESCALED_SAMPLES = list(
    set(hist_sm) - 
    set(HIST_RESCALED_SAMPLES))


###
# subsampled (prior to CpG filtering)
HIST_SUBSAMPLED_SAMPLES = list(
    set(hist_sm) & 
    set(config["subsampling_samplenames"]))

# not subsampled (prior to CpG filtering)
HIST_NOT_SUBSAMPLED_SAMPLES = list(
    set(hist_sm) - 
    set(config["subsampling_samplenames"]))


###
# rescaled and subsampled (prior to CpG filtering)
HIST_RESCALED_SUBSAMPLED_SAMPLES = list(
    set(HIST_RESCALED_SAMPLES) & 
    set(HIST_SUBSAMPLED_SAMPLES))

# rescaled, but not subsampled (prior to CpG filtering)
HIST_RESCALED_NOT_SUBSAMPLED_SAMPLES = list(
    set(HIST_RESCALED_SAMPLES) & 
    set(HIST_NOT_SUBSAMPLED_SAMPLES))

# not rescaled, but subsampled (prior to CpG filtering)
HIST_NOT_RESCALED_SUBSAMPLED_SAMPLES = list(
    set(HIST_NOT_RESCALED_SAMPLES) & 
    set(HIST_SUBSAMPLED_SAMPLES))

# neither rescaled nor subsampled (prior to CpG filtering)
HIST_NOT_RESCALED_NOT_SUBSAMPLED_SAMPLES = list(
    set(HIST_NOT_RESCALED_SAMPLES) & 
    set(HIST_NOT_SUBSAMPLED_SAMPLES))


###
# CpG filtered
HIST_CpG_SAMPLES = list(
    set(hist_sm) & 
    set(config["CpG_samplenames"]))

# not CpG filtered
HIST_NOT_CpG_SAMPLES = list(
    set(hist_sm) - 
    set(config["CpG_samplenames"]))


###
# rescaled, subsampled, CpG filtered
HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES = list(
    set(HIST_RESCALED_SAMPLES)
    & set(HIST_SUBSAMPLED_SAMPLES)
    & set(HIST_CpG_SAMPLES)
)

# rescaled, not subsampled, CpG filtered
HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES = list(
    set(HIST_RESCALED_SAMPLES)
    & set(HIST_NOT_SUBSAMPLED_SAMPLES)
    & set(HIST_CpG_SAMPLES)
)

# not rescaled, subsampled, CpG filtered
HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES = list(
    set(HIST_NOT_RESCALED_SAMPLES)
    & set(HIST_SUBSAMPLED_SAMPLES)
    & set(HIST_CpG_SAMPLES)
)

# not rescaled, not subsampled, CpG filtered
HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES = list(
    set(HIST_NOT_RESCALED_SAMPLES)
    & set(HIST_NOT_SUBSAMPLED_SAMPLES)
    & set(HIST_CpG_SAMPLES)
)


###
# rescaled, subsampled, not CpG filtered
HIST_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES = list(
    set(HIST_RESCALED_SAMPLES)
    & set(HIST_SUBSAMPLED_SAMPLES)
    & set(HIST_NOT_CpG_SAMPLES)
)

# rescaled, not subsampled, not CpG filtered
HIST_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES = list(
    set(HIST_RESCALED_SAMPLES)
    & set(HIST_NOT_SUBSAMPLED_SAMPLES)
    & set(HIST_NOT_CpG_SAMPLES)
)

# not rescaled, subsampled, not CpG filtered
HIST_NOT_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES = list(
    set(HIST_NOT_RESCALED_SAMPLES)
    & set(HIST_SUBSAMPLED_SAMPLES)
    & set(HIST_NOT_CpG_SAMPLES)
)

# not rescaled, not subsampled, not CpG filtered
HIST_NOT_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES = list(
    set(HIST_NOT_RESCALED_SAMPLES)
    & set(HIST_NOT_SUBSAMPLED_SAMPLES)
    & set(HIST_NOT_CpG_SAMPLES)
)


# Lists of modern samples
###
# subsampled
MODERN_SUBSAMPLED_SAMPLES = list(
    set(mod_sm) & 
    set(config["subsampling_samplenames"]))

# not subsampled
MODERN_NOT_SUBSAMPLED_SAMPLES = list(
    set(mod_sm) - 
    set(config["subsampling_samplenames"]))

###
# CpG filtered
MODERN_CpG_SAMPLES = list(
    set(mod_sm) & 
    set(config["CpG_samplenames"]))

# not CpG filtered
MODERN_NOT_CpG_SAMPLES = list(
    set(mod_sm) - 
    set(config["CpG_samplenames"]))

###
# subsampled and CpG filtered
MODERN_SUBSAMPLED_CpG_SAMPLES = list(
    set(MODERN_SUBSAMPLED_SAMPLES) & 
    set(MODERN_CpG_SAMPLES))

# subsampled but not CpG filtered
MODERN_SUBSAMPLED_NOT_CpG_SAMPLES = list(
    set(MODERN_SUBSAMPLED_SAMPLES) & 
    set(MODERN_NOT_CpG_SAMPLES))

# not subsampled but CpG filtered
MODERN_NOT_SUBSAMPLED_CpG_SAMPLES = list(
    set(MODERN_NOT_SUBSAMPLED_SAMPLES) & 
    set(MODERN_CpG_SAMPLES))

# not subsampled not CpG filtered
MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES = list(
    set(MODERN_NOT_SUBSAMPLED_SAMPLES) & 
    set(MODERN_NOT_CpG_SAMPLES))

# VCF file merging
ALL_SAMPLES = list(hist_sm + mod_sm)


###
# mlRho, merge VCFs: sex chromosomal scaffolds
sexchromosomeList = []  # fill the list with scaffold/contig names from the list of sex chromosome-linked scaffolds/contigs, if available
if os.path.exists(config["sexchromosomes"]):
    with open(config["sexchromosomes"], "r") as file:
        for line in file:
            sexchromosomeList.append(line.strip())

if len(sexchromosomeList) > 0:
    CHR = "autos"
elif len(sexchromosomeList) == 0:
    CHR = "genome"

###
# snpEff
if config["snpEff"]:
    if os.path.exists(config["gtf_path"]):
        GTF_DIR = os.path.dirname(config["gtf_path"])
        GTF_FILE = os.path.basename(config["gtf_path"])

###
# GERP
if config["gerp"]:
    # GERP input fasta path
    if config["gerp_ref_path"].endswith("/"):
        GERP_REF_PATH = config["gerp_ref_path"].rstrip("/")
    else:
        GERP_REF_PATH = config["gerp_ref_path"]

    # GERP input fasta file lists
    GERP_REF_FASTA_GZIP = [file for file in os.listdir(GERP_REF_PATH) if file.endswith(".gz")]
    GERP_REF_FASTA = [fasta.replace(".gz", "") for fasta in GERP_REF_FASTA_GZIP]
    GERP_REF_NAMES = [os.path.splitext(name)[0] for name in GERP_REF_FASTA]  # names of reference genomes of outgroup species
    ALL_GERP_REF_NAMES = GERP_REF_NAMES[:]
    ALL_GERP_REF_NAMES.append(REF_NAME)  # names of all genomes in the analysis, incl. the target species

    # Create list of chunk names for parallelization of GERP step
    # Adjust zero padding depending on the number of chunks set in the config file
    if config["gerp_chunks"] < 100:
        CHUNKS = ["chunk" + str(i).zfill(2) for i in range(1,config["gerp_chunks"]+1)]  # list of chunk names
    elif config["gerp_chunks"] >= 100:
        CHUNKS = ["chunk" + str(i).zfill(3) for i in range(1,config["gerp_chunks"]+1)]  # list of chunk names
    elif config["gerp_chunks"] >= 1000:
        CHUNKS = ["chunk" + str(i).zfill(4) for i in range(1,config["gerp_chunks"]+1)]  # list of chunk names
