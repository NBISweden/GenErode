##########################################################################
### 8.2 Quality filtering and repeat filtering of VCF files per sample

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    all_outputs.append("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/multiqc/multiqc_report.html")
    all_outputs.append("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/multiqc/multiqc_report.html")

if os.path.exists(config["modern_samples"]):
    all_outputs.append("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/multiqc/multiqc_report.html")
    all_outputs.append("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/multiqc/multiqc_report.html")


# Functions used by rules of this part of the pipeline
def depth_file_vcf(wildcards):
    """Select correct depth stats file for each sample"""
    if wildcards.sample in HIST_NOT_SUBSAMPLED_SAMPLES:
        dpstats = "results/historical/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample)
    elif wildcards.sample in MODERN_NOT_SUBSAMPLED_SAMPLES:
        dpstats = "results/modern/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample)
    elif wildcards.sample in HIST_NOT_RESCALED_SUBSAMPLED_SAMPLES:
        dpstats = "results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample, DP=config["subsampling_depth"])
    elif wildcards.sample in HIST_RESCALED_SUBSAMPLED_SAMPLES:
        dpstats = "results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample, DP=config["subsampling_depth"])
    elif wildcards.sample in MODERN_SUBSAMPLED_SAMPLES:
        dpstats = "results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample, DP=config["subsampling_depth"])
    return dpstats

def historical_quality_filtered_vcf_multiqc_inputs(wildcards):
    """Input for historical_quality_filtered_vcf_multiqc"""
    rescaled_not_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
        sample=HIST_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    not_rescaled_not_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    rescaled_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
        sample=HIST_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    not_rescaled_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
        sample=HIST_NOT_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    outlist = (rescaled_not_subsampled_not_CpG + not_rescaled_not_subsampled_not_CpG + rescaled_subsampled_not_CpG + not_rescaled_subsampled_not_CpG)
    if config["CpG_from_vcf"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.no{CpG_method}.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
            CpG_method="CpG_vcf",)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.Q30.sorted.no{CpG_method}.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
            CpG_method="CpG_vcf",)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.no{CpG_method}.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],
            CpG_method="CpG_vcf",)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.no{CpG_method}.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],
            CpG_method="CpG_vcf",)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG)
    elif config["CpG_from_reference"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG)
    elif config["CpG_from_vcf_and_reference"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG)
    return outlist

def modern_quality_filtered_vcf_multiqc_inputs(wildcards):
    """Input for modern_quality_filtered_vcf_multiqc"""
    not_subsampled_not_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
        sample=MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    subsampled_not_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
        sample=MODERN_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    outlist = (not_subsampled_not_CpG + subsampled_not_CpG)
    if config["CpG_from_vcf"] == True:
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.Q30.sorted.no{CpG_method}.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,
            CpG_method="CpG_vcf",)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.no{CpG_method}.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],
            CpG_method="CpG_vcf",)
        outlist += (not_subsampled_CpG + subsampled_CpG)
    elif config["CpG_from_reference"] == True:
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (not_subsampled_CpG + subsampled_CpG)
    elif config["CpG_from_vcf_and_reference"] == True:
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (not_subsampled_CpG + subsampled_CpG)
    return outlist

def historical_repmasked_vcf_multiqc_inputs(wildcards):
    """Input for historical_repmasked_vcf_multiqc"""
    rescaled_not_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
        sample=HIST_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    not_rescaled_not_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    rescaled_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
        sample=HIST_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    not_rescaled_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
        sample=HIST_NOT_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    outlist = (rescaled_not_subsampled_not_CpG + not_rescaled_not_subsampled_not_CpG + rescaled_subsampled_not_CpG + not_rescaled_subsampled_not_CpG)
    if config["CpG_from_vcf"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG)
    elif config["CpG_from_reference"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG)
    elif config["CpG_from_vcf_and_reference"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG)
    return outlist

def modern_repmasked_vcf_multiqc_inputs(wildcards):
    """Input for modern_repmasked_vcf_multiqc"""
    not_subsampled_not_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
        sample=MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    subsampled_not_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
        sample=MODERN_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    outlist = (not_subsampled_not_CpG + subsampled_not_CpG)
    if config["CpG_from_vcf"] == True:
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (not_subsampled_CpG + subsampled_CpG)
    elif config["CpG_from_reference"] == True:
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (not_subsampled_CpG + subsampled_CpG)
    elif config["CpG_from_vcf_and_reference"] == True:
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (not_subsampled_CpG + subsampled_CpG)
    return outlist


# snakemake rules
rule remove_snps_near_indels:
    """remove SNPs within 5 bp of an indel"""
    input:
        bcf="results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.bcf",
    output:
        snps=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.bcf"),
    threads: 2
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_remove_snps_near_indels.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools filter -g 5 -O b --threads {threads} -o {output.snps} {input.bcf} 2> {log}
        """


rule filter_vcfs_qual_dp:
    """Remove indels, genotypes of genotype quality < 30 and keep only sites within depth thresholds 
    that were determined from bam files earlier in the pipeline"""
    """Note that the depth filter is recalculated for subsampled bam files, according to the target depth for subsampling"""
    input:
        bcf=rules.remove_snps_near_indels.output,
        dp=depth_file_vcf,
    output:
        filtered=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.bcf"),
    threads: 2
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_filter_vcfs_qual_dp.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        minDP=`head -n 1 {input.dp} | cut -d' ' -f 2`
        maxDP=`head -n 1 {input.dp} | cut -d' ' -f 3`

        # check minimum depth threshold
        if awk "BEGIN{{exit ! ($minDP < 3)}}"
        then
          minDP=3
        fi

        bcftools filter -i "(DP4[0]+DP4[1]+DP4[2]+DP4[3])>$minDP & (DP4[0]+DP4[1]+DP4[2]+DP4[3])<$maxDP & QUAL>=30 & INDEL=0" -O b \
        --threads {threads} -o {output.filtered} {input.bcf} 2> {log}
        """


rule filter_vcfs_allelic_balance:
    """Removes heterozygote sites with allelic imbalance from the vcf files that could be due to contamination, sequencing or mapping errors"""
    """For example, sites where 9 reads support the reference allele and one read the alternative allele"""
    input:
        bcf=rules.filter_vcfs_qual_dp.output.filtered,
    output:
        filtered=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.bcf"),
    threads: 2
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_filter_vcfs_allelic_balance.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools view -e 'GT="0/1" & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) < 0.2' {input.bcf} | \
        bcftools view -e 'GT="0/1" & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]) > 0.8' -Ob > {output.filtered} 2> {log}
        """


rule index_filtered_vcfs:
    """Index vcf files before any downstream processing"""
    input:
        bcf=rules.filter_vcfs_allelic_balance.output.filtered,
    output:
        index=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.bcf.csi"),
    group:
        "qual_filtered_vcf_group"
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_index_filtered_vcfs.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools index -o {output.index} {input.bcf} 2> {log}
        """


rule filtered_vcf_stats:
    """Obtain summary stats of vcf files filtered for quality"""
    input:
        bcf=rules.filter_vcfs_allelic_balance.output.filtered,
        index=rules.index_filtered_vcfs.output.index,
    output:
        stats="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
    group:
        "qual_filtered_vcf_group"
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_filtered_vcf_stats.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools stats {input.bcf} > {output.stats} 2> {log}
        """


rule historical_quality_filtered_vcf_multiqc:
    """Collect all stats files from quality filtered historical vcf files"""
    input:
        historical_quality_filtered_vcf_multiqc_inputs,
    output:
        stats="results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/multiqc/multiqc_report.html",
    params:
        indir="results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/",
        outdir="results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/multiqc",
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/historical/" + REF_NAME + "/historical_quality_filtered_vcf_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_quality_filtered_vcf_multiqc:
    """Collect all stats files from quality filtered modern vcf files"""
    input:
        modern_quality_filtered_vcf_multiqc_inputs,
    output:
        stats="results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/multiqc/multiqc_report.html",
    params:
        indir="results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/",
        outdir="results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/multiqc",
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/modern/" + REF_NAME + "/modern_quality_filtered_vcf_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule filtered_bcf2vcf:
    """Convert bcf format to vcf.gz for removal of sites"""
    input:
        bcf=rules.filter_vcfs_allelic_balance.output.filtered,
    output:
        vcf=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.vcf.gz"),
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_filtered_bcf2vcf.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools convert -O z -o {output.vcf} {input.bcf} 2> {log}
        """


rule remove_repeats_vcf:
    """Remove repeats from vcf files"""
    input:
        vcf=rules.filtered_bcf2vcf.output.vcf,
        bed=rules.make_no_repeats_bed.output.no_rep_bed_dir,
        genomefile=rules.genome_file.output.genomefile,
    output:
        filtered=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.vcf.gz"),
    threads: 6
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_remove_repeats_vcf.log",
    singularity:
        "docker://verku/bedtools-2.29.2" # replace with link to NBIS Dockerhub repo
    shell:
        """
        bedtools intersect -a {input.vcf} -b {input.bed} -header -sorted -g {input.genomefile} | bgzip -c > {output.filtered} 2> {log}
        """


rule filtered_vcf2bcf:
    """Convert the repeat masked vcf back to bcf"""
    input:
        filtered=rules.remove_repeats_vcf.output.filtered,
    output:
        bcf="results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
    threads: 2
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_filtered_vcf2bcf.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools convert -O b -o {output.bcf} {input.filtered} 2> {log}
        """


rule index_repmasked_vcfs:
    """Index vcf files before any downstream processing"""
    input:
        bcf=rules.filtered_vcf2bcf.output.bcf,
    output:
        index="results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
    group:
        "repmasked_vcf_group"
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_index_repmasked_vcfs.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools index -o {output.index} {input.bcf} 2> {log}
        """


rule repmasked_vcf_stats:
    """Obtain summary stats of repeat masked vcf files"""
    input:
        bcf=rules.filtered_vcf2bcf.output.bcf,
        index=rules.index_repmasked_vcfs.output.index,
    output:
        stats="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
    group:
        "repmasked_vcf_group"
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{processed}_repmasked_vcf_stats.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools stats {input.bcf} > {output.stats} 2> {log}
        """


rule historical_repmasked_vcf_multiqc:
    """Collect all stats files from repeat masked historical vcf files"""
    input:
        historical_repmasked_vcf_multiqc_inputs,
    output:
        stats="results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/multiqc/multiqc_report.html",
    params:
        indir="results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/",
        outdir="results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/multiqc",
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/historical/" + REF_NAME + "/historical_repmasked_vcf_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_repmasked_vcf_multiqc:
    """Collect all stats files from repeat masked modern vcf files"""
    input:
        modern_repmasked_vcf_multiqc_inputs,
    output:
        stats="results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/multiqc/multiqc_report.html",
    params:
        indir="results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/",
        outdir="results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/multiqc",
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/modern/" + REF_NAME + "/modern_repmasked_vcf_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
