##########################################################################
##################### WORKFLOW VARIABLES AND CODE ########################
##########################################################################

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate
import os
import pandas as pd

min_version("9.0.0")

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

## Fix path to scratch directory
if os.path.exists(config["scratch_dir"]):
    if config["scratch_dir"].endswith("/"):
        scratch_dir = config["scratch_dir"]
    else:
        scratch_dir = config["scratch_dir"] + "/"
elif config["scratch_dir"] == ".":
    scratch_dir = config["scratch_dir"]
else:
    raise WorkflowError('No scratch directory found. Please add path to scratch directory in config.yaml file or set to "." for current directory.')


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
def check_metadata_file(dataframe):
    # raise an error if a sample exists with both fastq and bam files
    if "path_to_processed_bam_file" in dataframe.columns and "path_to_R1_fastq_file" in dataframe.columns:
        # check if any sample has both fastq and bam files
        fastq_samples = dataframe["samplename"][dataframe["path_to_R1_fastq_file"].notnull()].unique()
        bam_samples = dataframe["samplename"][dataframe["path_to_processed_bam_file"].notnull()].unique()
        # check if any sample has both fastq and bam files
        if set(fastq_samples) & set(bam_samples):
            # raise an error if any sample has both fastq and bam files
            raise WorkflowError("Both fastq files and bam files are provided for some samples. Please check your metadata file.")
    # add columns for merging of bam files
    if "samplename" in dataframe.columns and "library_id" in dataframe.columns and "lane" in dataframe.columns:
        # concatenate the sample name, library id and lane number with "_" to create a unique identifier for each fastq file
        dataframe["samplename_index_lane"] = dataframe["samplename"] + "_" + dataframe["library_id"] + "_" + dataframe["lane"]
        # set entries to NA if "samplename_index_lane" ends with "_" (schema validation checks for "_" in original columns)
        dataframe.loc[dataframe["samplename_index_lane"].fillna("").str.endswith("_"), "samplename_index_lane"] = pd.NA
        # concatenate the sample name and library id with "_" to create a unique identifier for merging of bam files
        dataframe["samplename_index"] = dataframe["samplename"] + "_" + dataframe["library_id"]
        # set entries to NA if "samplename_index" ends with "_" (schema validation checks for "_" in original columns)
        dataframe.loc[dataframe["samplename_index"].fillna("").str.endswith("_"), "samplename_index"] = pd.NA
    return dataframe

# Create sample lists
# sample list for each fastq-file (["samplename_index_lane"])
def samplename_index_lane_func(dataframe):
    if "samplename" in dataframe.columns and "library_id" in dataframe.columns and "lane" in dataframe.columns:
        if dataframe["samplename_index_lane"].duplicated().any():
            raise WorkflowError("Samples found with duplicate library ID and lane number. Please check your metadata file.")
        else:
            # convert the column into a list
            return list(dataframe["samplename_index_lane"].dropna())
    else:
        return []


# sample list for merging of files per lane (["samplename_index"])
def samplename_index_func(dataframe):
    if "samplename" in dataframe.columns and "library_id" in dataframe.columns and "lane" in dataframe.columns:
        return list(dataframe["samplename_index"].drop_duplicates().dropna())
    else:
        return []


# sample list (["samplename"])
def samplename_func(dataframe):
    return list(dataframe["samplename"].drop_duplicates().dropna())


# Functions to create symbolic links to fastq files, to look up read group information and for merging of bam files
# symbolic links dictionary
def fastq_symlinks_dict_func(dataframe):
    fastq_symlinks_dict = {}
    if "path_to_R1_fastq_file" in dataframe.columns and "path_to_R2_fastq_file" in dataframe.columns:
        for index, row in dataframe.iterrows():
            if pd.notnull(row["samplename_index_lane"]) and pd.notnull(row["path_to_R1_fastq_file"]) and pd.notnull(row["path_to_R2_fastq_file"]):
                # create a dictionary with the sample name as key and the path to the fastq files as value
                fastq_symlinks_dict[row["samplename_index_lane"]] = {"R1": os.path.abspath(row["path_to_R1_fastq_file"]), "R2": os.path.abspath(row["path_to_R2_fastq_file"])}
    return fastq_symlinks_dict


# read group dictionary
def rg_dict_func(dataframe):
    rg_dict = {}
    if "samplename_index_lane" in dataframe.columns and "readgroup_id" in dataframe.columns and "readgroup_platform" in dataframe.columns and "library_id" in dataframe.columns:
        for index, row in dataframe.iterrows():
            if pd.notnull(row["samplename_index_lane"]) and pd.notnull(row["readgroup_id"]) and pd.notnull(row["readgroup_platform"]) and pd.notnull(row["library_id"]):
                # create a dictionary with the sample name as key and the read group information as value
                rg_dict[row["samplename_index_lane"]] = {"ID": row["readgroup_id"], "SM": row["samplename"], "PL": row["readgroup_platform"], "LB": row["library_id"]}
    return rg_dict


# dictionary for bam file merging per lane
def sampleidxln_dict_func(dataframe):
    sampleidxln_dict = {}
    if "samplename_index_lane" in dataframe.columns and "samplename_index" in dataframe.columns:
        for index, row in dataframe.iterrows():
            if pd.notnull(row["samplename_index"]) and pd.notnull(row["samplename_index_lane"]):
                if (row["samplename_index"] in sampleidxln_dict):  # if "sample_index" is already in the dictionary
                    if (row["samplename_index_lane"] not in sampleidxln_dict[row["samplename_index"]]):  # if "sample_index_lane" is not in the list for "sample_index"
                        sampleidxln_dict[row["samplename_index"]].append(row["samplename_index_lane"])  # add "sample_index_lane" for "sample_index"
                else:  # if "sample_index" is not yet in the dictionary
                    sampleidxln_dict[row["samplename_index"]] = [row["samplename_index_lane"]]  # add "sample_index_lane" for "sample_index"
    return sampleidxln_dict


# dictionary for bam file merging per library ID
def sampleidx_dict_func(dataframe):
    sampleidx_dict = {}
    if "samplename_index" in dataframe.columns:
        for index, row in dataframe.iterrows():
            if pd.notnull(row["samplename_index_lane"]):
                if row["samplename"] in sampleidx_dict:  # if "sample" is already in the dictionary
                    if (row["samplename_index"] not in sampleidx_dict[row["samplename"]]):  # if "sample_index" is not in the list for "sample"
                        sampleidx_dict[row["samplename"]].append(row["samplename_index"])  # add "sample_index" for "sample"
                else:  # if "sample" is not yet in the dictionary
                    sampleidx_dict[row["samplename"]] = [row["samplename_index"]]  # add "sample_index" for "sample"
    return sampleidx_dict


# dictionary for bam file merging per sample (mitochondrial genomes)
def mito_sample_dict_func(dataframe):
    mito_sample_dict = {}
    if "samplename_index_lane" in dataframe.columns:
        fastq_samples = dataframe[dataframe["path_to_R1_fastq_file"].notnull()]
        for index, row in fastq_samples.iterrows():
            if pd.notnull(row["samplename_index_lane"]):
                if row["samplename"] in mito_sample_dict:  # if "sample" is already in the dictionary
                    if (row["samplename_index_lane"] not in mito_sample_dict[row["samplename"]]):  # if "sample_index_lane" is not in the list for "sample"
                        mito_sample_dict[row["samplename"]].append(row["samplename_index_lane"])  # add "sample_index_lane" for "sample"
                else:  # if "sample" is not yet in the dictionary
                    mito_sample_dict[row["samplename"]] = [row["samplename_index_lane"]]  # add "sample_index_lane" for "sample"
    return mito_sample_dict


# Functions to collect user-provided bam files
# symbolic links dictionaries
def user_bam_symlinks_dict_func(dataframe):
    if "path_to_processed_bam_file" in dataframe.columns:
        user_bam_symlinks_dict = {}
        for index, row in dataframe.iterrows():
            if pd.notnull(row["path_to_processed_bam_file"]):
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
    historical_df_checked = check_metadata_file(historical_df)  # check metadata file format
    hist_sm = samplename_func(historical_df_checked) # "samplename" for all samples
    hist_sm_idx = samplename_index_func(historical_df_checked) # "samplename_index" for all samples
    hist_sm_idx_ln = samplename_index_lane_func(historical_df_checked) # "samplename_index_lane" for all samples
    # for pipeline-processed fastq files
    hist_pipeline_bam_sm = pipeline_bam_samples_func(historical_df_checked) # "samplename" for pipeline-processed samples
    hist_pipeline_bam_sm_idx = [smid for smid in hist_sm_idx for sm in hist_pipeline_bam_sm if sm in smid] # "samplename_index" for pipeline-processed samples
    hist_pipeline_bam_sm_idx_ln = [smidln for smidln in hist_sm_idx_ln for sm in hist_pipeline_bam_sm if sm in smidln] # "samplename_index_lane" for pipeline-processed samples
    hist_fastq_symlinks_dict = fastq_symlinks_dict_func(historical_df_checked)
    hist_mito_sample_dict = mito_sample_dict_func(historical_df_checked) # mitogenome mapping
    hist_rg_dict = rg_dict_func(historical_df_checked)
    hist_sampleidxln_dict = sampleidxln_dict_func(historical_df_checked)
    hist_sampleidx_dict = sampleidx_dict_func(historical_df_checked)
    # for user-provided bam files
    hist_user_bam_symlinks_dict = user_bam_symlinks_dict_func(historical_df_checked)
    hist_user_bam_sm = user_bam_samples_func(historical_df_checked)

else:
    hist_sm = []
    hist_sm_idx = []
    hist_sm_idx_ln = []
    hist_pipeline_bam_sm = []
    hist_pipeline_bam_sm_idx = []
    hist_pipeline_bam_sm_idx_ln = []
    hist_fastq_symlinks_dict = {}
    hist_mito_sample_dict = {}
    hist_rg_dict = {}
    hist_sampleidxln_dict = {}
    hist_sampleidx_dict = {}
    hist_user_bam_symlinks_dict = {}
    hist_user_bam_sm = []


if os.path.exists(config["modern_samples"]):
    modern_df = pd.read_csv(config["modern_samples"], sep=";|,| |\t", engine='python', dtype=str)  # read in the metadata as dataframe
    validate(modern_df, schema="../schemas/metadata.schema.yaml") # validate metadata file format with JSON schema
    modern_df_checked = check_metadata_file(modern_df)  # check metadata file format
    mod_sm = samplename_func(modern_df_checked) # "samplename" for all samples
    mod_sm_idx = samplename_index_func(modern_df_checked) # "samplename_index" for all samples
    mod_sm_idx_ln = samplename_index_lane_func(modern_df_checked) # "samplename_index_lane" for all samples
    # for pipeline-processed fastq files
    mod_pipeline_bam_sm = pipeline_bam_samples_func(modern_df_checked) # "samplename" for pipeline-processed samples
    mod_pipeline_bam_sm_idx = [smid for smid in mod_sm_idx for sm in mod_pipeline_bam_sm if sm in smid] # "samplename_index" for pipeline-processed samples
    mod_pipeline_bam_sm_idx_ln = [smidln for smidln in mod_sm_idx_ln for sm in mod_pipeline_bam_sm if sm in smidln] # "samplename_index_lane" for pipeline-processed samples
    mod_fastq_symlinks_dict = fastq_symlinks_dict_func(modern_df_checked)
    mod_rg_dict = rg_dict_func(modern_df_checked)
    mod_sampleidxln_dict = sampleidxln_dict_func(modern_df_checked)
    mod_sampleidx_dict = sampleidx_dict_func(modern_df_checked)

    # for user-provided bam files
    mod_user_bam_symlinks_dict = user_bam_symlinks_dict_func(modern_df_checked)
    mod_user_bam_sm = user_bam_samples_func(modern_df_checked)
else:
    mod_sm = []
    mod_sm_idx = []
    mod_sm_idx_ln = []
    mod_pipeline_bam_sm = []
    mod_pipeline_bam_sm_idx = []
    mod_pipeline_bam_sm_idx_ln = []
    mod_fastq_symlinks_dict = {}
    mod_rg_dict = {}
    mod_sampleidxln_dict = {}
    mod_sampleidx_dict = {}
    mod_user_bam_symlinks_dict = {}
    mod_user_bam_sm = []


### Parameters regarding optional steps

# mitochondrial genome fasta variables for mapping analysis
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

# pipeline-processed BAM files
HIST_PIPELINE_RESCALED_SAMPLES = list(
    set(hist_pipeline_bam_sm) & 
    set(config["historical_rescaled_samplenames"]))
# user-provided BAM files
HIST_USER_RESCALED_SAMPLES = list(
    set(hist_user_bam_sm) & 
    set(config["historical_rescaled_samplenames"]))

# not rescaled (prior to subsampling and CpG filtering)
# pipeline-processed BAM files
HIST_PIPELINE_NOT_RESCALED_SAMPLES = list(
    set(hist_pipeline_bam_sm) - 
    set(HIST_PIPELINE_RESCALED_SAMPLES))
# user-provided BAM files
HIST_USER_NOT_RESCALED_SAMPLES = list(
    set(hist_user_bam_sm) - 
    set(HIST_USER_RESCALED_SAMPLES))

###
# subsampled (prior to CpG filtering)
HIST_SUBSAMPLED_SAMPLES = list(
    set(hist_sm) & 
    set(config["subsampling_samplenames"]))

# pipeline-processed BAM files
HIST_PIPELINE_SUBSAMPLED_SAMPLES = list(
    set(hist_pipeline_bam_sm) & 
    set(config["subsampling_samplenames"]))
# user-provided BAM files
HIST_USER_SUBSAMPLED_SAMPLES = list(
    set(hist_user_bam_sm) & 
    set(config["subsampling_samplenames"]))

# not subsampled (prior to CpG filtering)
# pipeline-processed BAM files
HIST_PIPELINE_NOT_SUBSAMPLED_SAMPLES = list(
    set(hist_pipeline_bam_sm) - 
    set(HIST_PIPELINE_SUBSAMPLED_SAMPLES))
# user-provided BAM files
HIST_USER_NOT_SUBSAMPLED_SAMPLES = list(
    set(hist_user_bam_sm) - 
    set(HIST_USER_SUBSAMPLED_SAMPLES))

###
# rescaled and subsampled (prior to CpG filtering)
# pipeline-processed BAM files
HIST_PIPELINE_RESCALED_SUBSAMPLED_SAMPLES = list(
    set(HIST_PIPELINE_RESCALED_SAMPLES) & 
    set(HIST_PIPELINE_SUBSAMPLED_SAMPLES))
# user-provided BAM files
HIST_USER_RESCALED_SUBSAMPLED_SAMPLES = list(
    set(HIST_USER_RESCALED_SAMPLES) & 
    set(HIST_USER_SUBSAMPLED_SAMPLES))

# rescaled, but not subsampled (prior to CpG filtering)
# pipeline-processed BAM files
HIST_PIPELINE_RESCALED_NOT_SUBSAMPLED_SAMPLES = list(
    set(HIST_PIPELINE_RESCALED_SAMPLES) & 
    set(HIST_PIPELINE_NOT_SUBSAMPLED_SAMPLES))
# user-provided BAM files
HIST_USER_RESCALED_NOT_SUBSAMPLED_SAMPLES = list(
    set(HIST_USER_RESCALED_SAMPLES) & 
    set(HIST_USER_NOT_SUBSAMPLED_SAMPLES))

# not rescaled, but subsampled (prior to CpG filtering)
# pipeline-processed BAM files
HIST_PIPELINE_NOT_RESCALED_SUBSAMPLED_SAMPLES = list(
    set(HIST_PIPELINE_NOT_RESCALED_SAMPLES) & 
    set(HIST_PIPELINE_SUBSAMPLED_SAMPLES))
# user-provided BAM files
HIST_USER_NOT_RESCALED_SUBSAMPLED_SAMPLES = list(
    set(HIST_USER_NOT_RESCALED_SAMPLES) & 
    set(HIST_USER_SUBSAMPLED_SAMPLES))

# neither rescaled nor subsampled (prior to CpG filtering)
# pipeline-processed BAM files
HIST_PIPELINE_NOT_RESCALED_NOT_SUBSAMPLED_SAMPLES = list(
    set(HIST_PIPELINE_NOT_RESCALED_SAMPLES) & 
    set(HIST_PIPELINE_NOT_SUBSAMPLED_SAMPLES))
# user-provided BAM files
HIST_USER_NOT_RESCALED_NOT_SUBSAMPLED_SAMPLES = list(
    set(HIST_USER_NOT_RESCALED_SAMPLES) & 
    set(HIST_USER_NOT_SUBSAMPLED_SAMPLES))

###
# CpG filtered
HIST_CpG_SAMPLES = list(
    set(hist_sm) & 
    set(config["CpG_samplenames"]))

# not CpG filtered
HIST_NOT_CpG_SAMPLES = list(
    set(hist_sm) - 
    set(config["CpG_samplenames"]))


# Lists of modern samples
###
# subsampled
MOD_SUBSAMPLED_SAMPLES = list(
    set(mod_sm) & 
    set(config["subsampling_samplenames"]))

# pipeline-processed BAM files
MOD_PIPELINE_SUBSAMPLED_SAMPLES = list(
    set(mod_pipeline_bam_sm) & 
    set(config["subsampling_samplenames"]))
# user-provided BAM files
MOD_USER_SUBSAMPLED_SAMPLES = list(
    set(mod_user_bam_sm) & 
    set(config["subsampling_samplenames"]))

# pipeline-processed BAM files
MOD_PIPELINE_NOT_SUBSAMPLED_SAMPLES = list(
    set(mod_pipeline_bam_sm) - 
    set(MOD_PIPELINE_SUBSAMPLED_SAMPLES))
# user-provided BAM files
MOD_USER_NOT_SUBSAMPLED_SAMPLES = list(
    set(mod_user_bam_sm) - 
    set(MOD_USER_SUBSAMPLED_SAMPLES))

###
# CpG filtered
MOD_CpG_SAMPLES = list(
    set(mod_sm) & 
    set(config["CpG_samplenames"]))

# not CpG filtered
MOD_NOT_CpG_SAMPLES = list(
    set(mod_sm) - 
    set(config["CpG_samplenames"]))

###
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

### mlRho, genotyping, and VCF quality and repeat filtering input files
def processed_bam_file(wildcards):
    """Select correct bam file for each sample"""
    # pipeline-processed historical samples
    if wildcards.sample in HIST_PIPELINE_NOT_RESCALED_NOT_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam".format(
            sample=wildcards.sample,)
    elif wildcards.sample in HIST_PIPELINE_RESCALED_NOT_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.bam".format(
            sample=wildcards.sample,)
    elif wildcards.sample in HIST_PIPELINE_NOT_RESCALED_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.bam".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])
    elif wildcards.sample in HIST_PIPELINE_RESCALED_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.bam".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])
    # pipeline-processed modern samples
    elif wildcards.sample in MOD_PIPELINE_NOT_SUBSAMPLED_SAMPLES:
        return "results/modern/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam".format(
            sample=wildcards.sample,)
    elif wildcards.sample in MOD_PIPELINE_SUBSAMPLED_SAMPLES:
        return "results/modern/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.bam".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])
    # user-provided historical samples
    elif wildcards.sample in HIST_USER_NOT_RESCALED_NOT_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.userprovided.bam".format(
            sample=wildcards.sample,)
    elif wildcards.sample in HIST_USER_RESCALED_NOT_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.userprovided.rescaled.bam".format(
            sample=wildcards.sample,)
    elif wildcards.sample in HIST_USER_NOT_RESCALED_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.userprovided.mapped_q30.subs_dp{DP}.bam".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])
    elif wildcards.sample in HIST_USER_RESCALED_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.userprovided.rescaled.mapped_q30.subs_dp{DP}.bam".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])
    # user-provided modern samples
    elif wildcards.sample in MOD_USER_NOT_SUBSAMPLED_SAMPLES:
        return "results/modern/mapping/" + REF_NAME + "/{sample}.userprovided.bam".format(
            sample=wildcards.sample,)
    elif wildcards.sample in MOD_USER_SUBSAMPLED_SAMPLES:
        return "results/modern/mapping/" + REF_NAME + "/{sample}.userprovided.mapped_q30.subs_dp{DP}.bam".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])


def depth_file(wildcards):
    """Select correct depth stats file for each sample"""
    # pipeline-processed historical samples
    if wildcards.sample in HIST_PIPELINE_NOT_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dpstats.txt".format(
            sample=wildcards.sample,)
    elif wildcards.sample in HIST_PIPELINE_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])
    # pipeline-processed modern samples
    elif wildcards.sample in MOD_PIPELINE_NOT_SUBSAMPLED_SAMPLES:
        return "results/modern/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dpstats.txt".format(
            sample=wildcards.sample,)
    elif wildcards.sample in MOD_PIPELINE_SUBSAMPLED_SAMPLES:
        return "results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])
    # user-provided historical samples
    elif wildcards.sample in HIST_USER_NOT_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.repma.Q30.bam.dpstats.txt".format(
            sample=wildcards.sample,)
    elif wildcards.sample in HIST_USER_SUBSAMPLED_SAMPLES:
        return "results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.userprovided.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])
    # user-provided modern samples
    elif wildcards.sample in MOD_USER_NOT_SUBSAMPLED_SAMPLES:
        return "results/modern/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.repma.Q30.bam.dpstats.txt".format(
            sample=wildcards.sample,)
    elif wildcards.sample in MOD_USER_SUBSAMPLED_SAMPLES:
        return "results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.userprovided.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(
            sample=wildcards.sample,
            DP=config["subsampling_depth"])

###
# snpEff
if os.path.exists(config["gtf_path"]):
    GTF_DIR = os.path.dirname(config["gtf_path"])
    GTF_FILE = os.path.basename(config["gtf_path"])
else:
    GTF_DIR = config["gtf_path"]
    GTF_FILE = config["gtf_path"]

###
# GERP
# GERP input fasta path
if config["gerp_ref_path"].endswith("/"):
    GERP_REF_PATH = config["gerp_ref_path"].rstrip("/")
else:
    GERP_REF_PATH = config["gerp_ref_path"]

# GERP input fasta file lists
if os.path.exists(GERP_REF_PATH):
    GERP_REF_FASTA_GZIP = [file for file in os.listdir(GERP_REF_PATH) if file.endswith(".gz")]
    GERP_REF_FASTA = [fasta.replace(".gz", "") for fasta in GERP_REF_FASTA_GZIP]
    GERP_REF_NAMES = [os.path.splitext(name)[0] for name in GERP_REF_FASTA]  # names of reference genomes of outgroup species
    ALL_GERP_REF_NAMES = GERP_REF_NAMES[:]
    ALL_GERP_REF_NAMES.append(REF_NAME)  # names of all genomes in the analysis, incl. the target species
else:
    GERP_REF_FASTA_GZIP = []
    GERP_REF_FASTA = []
    GERP_REF_NAMES = []  # names of reference genomes of outgroup species
    ALL_GERP_REF_NAMES = GERP_REF_NAMES

# Create list of chunk names for parallelization of GERP step
# Adjust zero padding depending on the number of chunks set in the config file
if config["gerp_chunks"] < 100:
    CHUNKS = ["chunk" + str(i).zfill(2) for i in range(1,config["gerp_chunks"]+1)]  # list of chunk names
elif config["gerp_chunks"] >= 100:
    CHUNKS = ["chunk" + str(i).zfill(3) for i in range(1,config["gerp_chunks"]+1)]  # list of chunk names
elif config["gerp_chunks"] >= 1000:
    CHUNKS = ["chunk" + str(i).zfill(4) for i in range(1,config["gerp_chunks"]+1)]  # list of chunk names
else:
    CHUNKS = []