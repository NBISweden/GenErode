##########################################################################
# This is the Snakefile of the GenErode pipeline for historical or       #
# ancient and modern samples to study patterns of genome erosion         #
#                                                                        #
# Pipeline version 0.5.1                                                 #
#                                                                        #
# Written by Verena Kutschera, Marcin Kierczak and Tom van der Valk      #
# Email: generode@nbis.se                                                #
##########################################################################

##########################################################################
######################### HOW TO RUN THE WORKFLOW ########################
##########################################################################

###
# The full pipeline documentation can be found on the GitHub wiki page
# (https://github.com/NBISweden/GenErode/wiki).
###


##########################################################################
########################### STEPS TO INCLUDE #############################
##########################################################################

include: "workflow/rules/common.smk"


all_outputs = []

if config["reference_repeat_identification"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"


if config["fastq_processing"]:
    include: "workflow/rules/1.1_fastq_processing.smk"


if config["map_historical_to_mitogenomes"]:
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/1.2_map_to_mitogenomes.smk"


if config["mapping"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"


if config["bam_rmdup_realign_indels"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"


if config["historical_bam_mapDamage"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"


if config["bam_subsampling"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"


if config["genotyping"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"
    include: "workflow/rules/4_genotyping.smk"


if config["CpG_identification"]:
    if config["CpG_from_vcf"] == True:
        include: "workflow/rules/0.1_reference_genome_preps.smk"
        include: "workflow/rules/0.2_repeat_identification.smk"
        include: "workflow/rules/1.1_fastq_processing.smk"
        include: "workflow/rules/2_mapping.smk"
        include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
        include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
        include: "workflow/rules/3.3_bam_subsampling.smk"
        include: "workflow/rules/4_genotyping.smk"
        include: "workflow/rules/5_CpG_identification.smk"


    elif config["CpG_from_reference"] == True:
        include: "workflow/rules/0.1_reference_genome_preps.smk"
        include: "workflow/rules/0.2_repeat_identification.smk"
        include: "workflow/rules/1.1_fastq_processing.smk"
        include: "workflow/rules/2_mapping.smk"
        include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
        include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
        include: "workflow/rules/3.3_bam_subsampling.smk"
        include: "workflow/rules/5_CpG_identification.smk"


    elif config["CpG_from_vcf_and_reference"] == True:
        include: "workflow/rules/0.1_reference_genome_preps.smk"
        include: "workflow/rules/0.2_repeat_identification.smk"
        include: "workflow/rules/1.1_fastq_processing.smk"
        include: "workflow/rules/2_mapping.smk"
        include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
        include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
        include: "workflow/rules/3.3_bam_subsampling.smk"
        include: "workflow/rules/4_genotyping.smk"
        include: "workflow/rules/5_CpG_identification.smk"


###
if config["autosome_sexchromosome_bed_files"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"
    include: "workflow/rules/6_autosome_sexchromosome_bed_files.smk"


if config["mlRho"]:
    if len(config["CpG_samplenames"]) > 0:  # to avoid genotyping if not necessary
        if config["CpG_from_vcf"] == True:
            include: "workflow/rules/0.1_reference_genome_preps.smk"
            include: "workflow/rules/0.2_repeat_identification.smk"
            include: "workflow/rules/1.1_fastq_processing.smk"
            include: "workflow/rules/2_mapping.smk"
            include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
            include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
            include: "workflow/rules/3.3_bam_subsampling.smk"
            include: "workflow/rules/4_genotyping.smk"
            include: "workflow/rules/5_CpG_identification.smk"
            include: "workflow/rules/6_autosome_sexchromosome_bed_files.smk"
            include: "workflow/rules/7_mlRho.smk"

        elif config["CpG_from_reference"] == True:
            include: "workflow/rules/0.1_reference_genome_preps.smk"
            include: "workflow/rules/0.2_repeat_identification.smk"
            include: "workflow/rules/1.1_fastq_processing.smk"
            include: "workflow/rules/2_mapping.smk"
            include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
            include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
            include: "workflow/rules/3.3_bam_subsampling.smk"
            include: "workflow/rules/5_CpG_identification.smk"
            include: "workflow/rules/6_autosome_sexchromosome_bed_files.smk"
            include: "workflow/rules/7_mlRho.smk"

        elif config["CpG_from_vcf_and_reference"] == True:
            include: "workflow/rules/0.1_reference_genome_preps.smk"
            include: "workflow/rules/0.2_repeat_identification.smk"
            include: "workflow/rules/1.1_fastq_processing.smk"
            include: "workflow/rules/2_mapping.smk"
            include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
            include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
            include: "workflow/rules/3.3_bam_subsampling.smk"
            include: "workflow/rules/4_genotyping.smk"
            include: "workflow/rules/5_CpG_identification.smk"
            include: "workflow/rules/6_autosome_sexchromosome_bed_files.smk"
            include: "workflow/rules/7_mlRho.smk"

    elif len(config["CpG_samplenames"]) == 0:
        include: "workflow/rules/0.1_reference_genome_preps.smk"
        include: "workflow/rules/0.2_repeat_identification.smk"
        include: "workflow/rules/1.1_fastq_processing.smk"
        include: "workflow/rules/2_mapping.smk"
        include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
        include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
        include: "workflow/rules/3.3_bam_subsampling.smk"
        include: "workflow/rules/6_autosome_sexchromosome_bed_files.smk"
        include: "workflow/rules/7_mlRho.smk"


###
if config["vcf_CpG_filtering"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"
    include: "workflow/rules/4_genotyping.smk"
    include: "workflow/rules/5_CpG_identification.smk"
    include: "workflow/rules/8.1_vcf_CpG_filtering.smk"


if config["vcf_qual_repeat_filtering"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"
    include: "workflow/rules/4_genotyping.smk"
    include: "workflow/rules/5_CpG_identification.smk"
    include: "workflow/rules/8.1_vcf_CpG_filtering.smk"
    include: "workflow/rules/8.2_vcf_qual_repeat_filtering.smk"


if config["merge_vcfs_per_dataset"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"
    include: "workflow/rules/4_genotyping.smk"
    include: "workflow/rules/5_CpG_identification.smk"
    include: "workflow/rules/8.1_vcf_CpG_filtering.smk"
    include: "workflow/rules/8.2_vcf_qual_repeat_filtering.smk"
    include: "workflow/rules/9_merge_vcfs.smk"


if config["pca"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"
    include: "workflow/rules/4_genotyping.smk"
    include: "workflow/rules/5_CpG_identification.smk"
    include: "workflow/rules/8.1_vcf_CpG_filtering.smk"
    include: "workflow/rules/8.2_vcf_qual_repeat_filtering.smk"
    include: "workflow/rules/9_merge_vcfs.smk"
    include: "workflow/rules/10_pca.smk"


if config["ROH"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"
    include: "workflow/rules/4_genotyping.smk"
    include: "workflow/rules/5_CpG_identification.smk"
    include: "workflow/rules/8.1_vcf_CpG_filtering.smk"
    include: "workflow/rules/8.2_vcf_qual_repeat_filtering.smk"
    include: "workflow/rules/9_merge_vcfs.smk"
    include: "workflow/rules/11_ROH.smk"


if config["snpEff"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"
    include: "workflow/rules/4_genotyping.smk"
    include: "workflow/rules/5_CpG_identification.smk"
    include: "workflow/rules/8.1_vcf_CpG_filtering.smk"
    include: "workflow/rules/8.2_vcf_qual_repeat_filtering.smk"
    include: "workflow/rules/9_merge_vcfs.smk"
    include: "workflow/rules/12_snpEff.smk"


###
if config["gerp"]:
    include: "workflow/rules/0.1_reference_genome_preps.smk"
    include: "workflow/rules/0.2_repeat_identification.smk"
    include: "workflow/rules/1.1_fastq_processing.smk"
    include: "workflow/rules/2_mapping.smk"
    include: "workflow/rules/3.1_bam_rmdup_realign_indels.smk"
    include: "workflow/rules/3.2_historical_bam_mapDamage.smk"
    include: "workflow/rules/3.3_bam_subsampling.smk"
    include: "workflow/rules/4_genotyping.smk"
    include: "workflow/rules/5_CpG_identification.smk"
    include: "workflow/rules/8.1_vcf_CpG_filtering.smk"
    include: "workflow/rules/8.2_vcf_qual_repeat_filtering.smk"
    include: "workflow/rules/9_merge_vcfs.smk"
    include: "workflow/rules/13_GERP.smk"
###


##########################################################################
############################# PSEUDORULES ################################
##########################################################################


rule all:
    input:
        all_outputs,


##########################################################################
################################ REPORT ##################################
##########################################################################


if (config["bam_rmdup_realign_indels"]
    or config["mlRho"]
    or config["pca"]
    or config["ROH"]
    or config["snpEff"]
    or config["gerp"]):

    if os.path.exists(config["historical_samples"]) or os.path.exists(config["modern_samples"]):
        onsuccess:
            shell(
                r"""
                snakemake --unlock --cores 1 &&
                snakemake --report GenErode_pipeline_report.html --cores 1
                """)
            all_outputs.append("GenErode_pipeline_report.html")
