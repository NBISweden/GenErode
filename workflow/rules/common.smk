##########################################################################
##################### WORKFLOW VARIABLES AND CODE ########################
##########################################################################

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate
import os
import pandas as pd

min_version("7.0.0")

configfile: "config/config.yaml"

report: "../report/workflow.rst"

## Variables for software containers for easier version updating
bwa_container = "docker://biocontainers/bwa:v0.7.17-3-deb_cv1"
bwa_samtools_container = "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
picard_container = "docker://quay.io/biocontainers/picard:2.26.6--hdfd78af_0"
repeatmodeler_container = "https://depot.galaxyproject.org/singularity/repeatmodeler:2.0.5--pl5321hdfd78af_0"
bedtools_htslib_container = "oras://community.wave.seqera.io/library/bedtools_htslib:06ed4722f423d939"
fastqc_container = "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
multiqc_container = "oras://community.wave.seqera.io/library/multiqc:1.28--d466e41d58d6d704"
fastp_container = "docker://quay.io/biocontainers/fastp:0.24.0--h125f33a_0"
qualimap_container = "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
samtools_python_container = "oras://community.wave.seqera.io/library/samtools_python:2e56d0f345426c81"
gatk3_container = "docker://broadinstitute/gatk3:3.7-0"
mapdamage_container = "docker://biocontainers/mapdamage:v2.0.9dfsg-1-deb_cv1"
bcftools_container = "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0"
mlrho_container = "docker://nbisweden/generode-mlrho"
plink_container = "https://depot.galaxyproject.org/singularity/plink:1.90b6.12--heea4ae3_0"
vcftools_container = "docker://biocontainers/vcftools:v0.1.16-1-deb_cv1"
snpeff_container = "docker://quay.io/biocontainers/snpeff:4.3.1t--3"
seqtk_container = "oras://community.wave.seqera.io/library/seqtk:1.4--e75a8dec899d1be8"
gerp_container = "https://depot.galaxyproject.org/singularity/gerp:2.1--h1b792b2_2"
# shell_container = "oras://community.wave.seqera.io/library/biopython_matplotlib_numpy_pandas_python:3da9b5da9b1e30c6"


## Reference assembly variables
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
# Create sample lists
# sample list for each fastq-file (["samplename_index_lane"])
def samplename_index_lane_list_func(dataframe):
    # convert the first column into a list
    samplename_index_lane = list(dataframe["samplename_index_lane"])
    sm_idx_ln = list(sorted(set(samplename_index_lane)))  # keep only unique entries and sort the lists
    return sm_idx_ln


# sample list for merging of files per lane (["samplename_index"])
def samplename_index_list_func(dataframe):
    char = "_"
    samplename_index = []
    for i in list(dataframe["samplename_index_lane"]): # loop through the first column of the dataframe (converted to a list)
        smid = char.join(i.split(char)[:2])  # split the string by "_" and join the elements up to the second "_" by "_" to get sample name plus index
        samplename_index.append(smid)  # add samplename_index to the list
    sm_idx = list(sorted(set(samplename_index)))  # keep only unique entries and sort the lists
    return sm_idx


# sample list (["samplename"])
def samplename_list_func(dataframe):
    char = "_"
    samplename = []
    for i in list(dataframe["samplename_index_lane"]):
        sm = i.split(char)[0]  # split the string by "_" and take the first element to get the sample name
        samplename.append(sm)  # add the sample name to the list of sample names
    sm = list(sorted(set(samplename)))  # keep only unique entries and sort the lists
    return sm


# Functions to create symbolic links to fastq files, to look up read group information and for merging of bam files
# symbolic links dictionary
def symlinks_dict_func(dataframe):
    symlinks_dict = {}
    for index, row in dataframe.iterrows():
        symlinks_dict[row["samplename_index_lane"]] = {"R1": os.path.abspath(row["path_to_R1_fastq_file"]), "R2": os.path.abspath(row["path_to_R2_fastq_file"])}
    return symlinks_dict


# read group dictionary
def rg_dict_func(dataframe):
    char = "_"
    rg_dict = {}
    for index, row in dataframe.iterrows():
        sm = row["samplename_index_lane"].split(char)[0] # take the sample name from the first column of each line
        lb = row["samplename_index_lane"].split(char)[1]  # take the library id from the first column of each line
        rg_dict[row["samplename_index_lane"]] = {"ID": row["readgroup_id"], "SM": sm, "PL": row["readgroup_platform"], "LB": lb}
    return rg_dict


# dictionary for bam file merging per lane
def sampleidxln_dict_func(dataframe):
    char = "_"
    sampleidxln_dict = {}
    for index, row in dataframe.iterrows():
        smid = char.join(row["samplename_index_lane"].split(char)[:2])  # take the sample name plus index from the first column of each line
        if (smid in sampleidxln_dict):  # if "sample_index" is already in the dictionary
            if (row["samplename_index_lane"] not in sampleidxln_dict[smid]):  # if "sample_index_lane" is not in the list for "sample_index"
                sampleidxln_dict[smid].append(row["samplename_index_lane"])  # add "sample_index_lane" for "sample_index"
        else:  # if "sample_index" is not yet in the dictionary
            sampleidxln_dict[smid] = [row["samplename_index_lane"]]  # add "sample_index_lane" for "sample_index"
    return sampleidxln_dict


# dictionary for bam file merging per PCR
def sampleidx_dict_func(dataframe):
    char = "_"
    sampleidx_dict = {}
    for index, row in dataframe.iterrows():
        sm = row["samplename_index_lane"].split(char)[0]  # take the sample name from the first column of each line
        smid = char.join(row["samplename_index_lane"].split(char)[:2])  # take the sample name plus index from the first column of each line
        if sm in sampleidx_dict:  # if "sample" is already in the dictionary
            if (smid not in sampleidx_dict[sm]):  # if "sample_index" is not in the list for "sample"
                sampleidx_dict[sm].append(smid)  # add "sample_index" for "sample"
        else:  # if "sample" is not yet in the dictionary
            sampleidx_dict[sm] = [smid]  # add "sample_index" for "sample"
    return sampleidx_dict


# dictionary for bam file merging per sample (mitochondrial genomes)
def sample_dict_func(dataframe):
    char = "_"
    sample_dict = {}
    for index, row in dataframe.iterrows():
        sm = row["samplename_index_lane"].split(char)[0]  # take the sample name from the first column of each line
        if sm in sample_dict:  # if "sample" is already in the dictionary
            if (row["samplename_index_lane"] not in sample_dict[sm]):  # if "sample_index_lane" is not in the list for "sample"
                sample_dict[sm].append(row["samplename_index_lane"])  # add "sample_index_lane" for "sample"
        else:  # if "sample" is not yet in the dictionary
            sample_dict[sm] = [row["samplename_index_lane"]]  # add "sample_index_lane" for "sample"
    return sample_dict


# Apply the functions to metadata tables for historical and modern samples
if os.path.exists(config["historical_samples"]):
    historical_df = pd.read_csv(config["historical_samples"], sep=";|,| |\t", engine='python')  # read in the metadata as dataframe
    validate(historical_df, schema="../schemas/metadata.schema.yaml")  # validate metadata file format with JSON schema
    hist_sm = samplename_list_func(historical_df)
    hist_sm_idx = samplename_index_list_func(historical_df)
    hist_sm_idx_ln = samplename_index_lane_list_func(historical_df)
    hist_symlinks_dict = symlinks_dict_func(historical_df)
    hist_sample_dict = sample_dict_func(historical_df)
    hist_rg_dict = rg_dict_func(historical_df)
    hist_sampleidxln_dict = sampleidxln_dict_func(historical_df)
    hist_sampleidx_dict = sampleidx_dict_func(historical_df)
else:
    hist_sm = []
    hist_sm_idx = []
    hist_sm_idx_ln = []
    hist_symlinks_dict = {}
    hist_sample_dict = {}
    hist_rg_dict = {}
    hist_sampleidxln_dict = {}
    hist_sampleidx_dict = {}


if os.path.exists(config["modern_samples"]):
    modern_df = pd.read_csv(config["modern_samples"], sep=";|,| |\t", engine='python')  # read in the metadata as dataframe
    validate(modern_df, schema="../schemas/metadata.schema.yaml") # validate metadata file format with JSON schema
    mod_sm = samplename_list_func(modern_df)
    mod_sm_idx = samplename_index_list_func(modern_df)
    mod_sm_idx_ln = samplename_index_lane_list_func(modern_df)
    mod_symlinks_dict = symlinks_dict_func(modern_df)
    mod_rg_dict = rg_dict_func(modern_df)
    mod_sampleidxln_dict = sampleidxln_dict_func(modern_df)
    mod_sampleidx_dict = sampleidx_dict_func(modern_df)
else:
    mod_sm = []
    mod_sm_idx = []
    mod_sm_idx_ln = []
    mod_symlinks_dict = {}
    mod_rg_dict = {}
    mod_sampleidxln_dict = {}
    mod_sampleidx_dict = {}


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
