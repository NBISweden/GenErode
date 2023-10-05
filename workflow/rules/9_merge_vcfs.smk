##########################################################################
### 9. Merge VCF files and filter for biallelic sites, missingness and sex-chromosomal contigs/scaffolds

# Code collecting output files from this part of the pipeline
all_outputs.append("results/all/vcf/" + REF_NAME + "/stats/vcf_merged_missing/multiqc/multiqc_report.html")


# Functions used by rules of this part of the pipeline
def merge_all_inputs(wildcards):
    """Input for merge_all_vcfs"""
    rescaled_not_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
        sample=HIST_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    not_rescaled_not_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    rescaled_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
        sample=HIST_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    not_rescaled_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
        sample=HIST_NOT_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    not_subsampled_not_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
        sample=MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    subsampled_not_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
        sample=MODERN_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    outlist = (rescaled_not_subsampled_not_CpG + not_rescaled_not_subsampled_not_CpG + rescaled_subsampled_not_CpG + not_rescaled_subsampled_not_CpG + not_subsampled_not_CpG + subsampled_not_CpG)
    if config["CpG_from_vcf"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG + not_subsampled_CpG + subsampled_CpG)
    elif config["CpG_from_reference"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG + not_subsampled_CpG + subsampled_CpG)
    elif config["CpG_from_vcf_and_reference"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG + not_subsampled_CpG + subsampled_CpG)
    return outlist


def merge_all_index_inputs(wildcards):
    """Input for merge_all_vcfs"""
    rescaled_not_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
        sample=HIST_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    not_rescaled_not_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    rescaled_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
        sample=HIST_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    not_rescaled_subsampled_not_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
        sample=HIST_NOT_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    not_subsampled_not_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
        sample=MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
    subsampled_not_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
        sample=MODERN_SUBSAMPLED_NOT_CpG_SAMPLES,
        DP=config["subsampling_depth"],)
    outlist = (rescaled_not_subsampled_not_CpG + not_rescaled_not_subsampled_not_CpG + rescaled_subsampled_not_CpG + not_rescaled_subsampled_not_CpG + not_subsampled_not_CpG + subsampled_not_CpG)
    if config["CpG_from_vcf"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG + not_subsampled_CpG + subsampled_CpG)
    elif config["CpG_from_reference"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG + not_subsampled_CpG + subsampled_CpG)
    elif config["CpG_from_vcf_and_reference"] == True:
        rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        not_rescaled_not_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
        rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_rescaled_subsampled_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        not_subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
        subsampled_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
            DP=config["subsampling_depth"],)
        outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG + not_subsampled_CpG + subsampled_CpG)
    return outlist

def missingness_filtered_vcf_multiqc_inputs(wildcards):
    """Input for missingness_filtered_vcf_multiqc"""
    if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
        return expand("results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_merged_missing/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            dataset=["all", "historical", "modern"],
            fmiss=config["f_missing"],
            chr=CHR,)
    elif os.path.exists(config["historical_samples"]):
        return expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_merged_missing/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            fmiss=config["f_missing"],
            chr=CHR,)
    elif os.path.exists(config["modern_samples"]):
        return expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_merged_missing/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            fmiss=config["f_missing"],
            chr=CHR,)


# snakemake rules
rule merge_all_vcfs:
    """Merge all samples into one VCF file, containing only SNPs"""
    input:
        bcf=merge_all_inputs,
        index=merge_all_index_inputs,
    output:
        merged=temp("results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf"),
    threads: 6
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + "_merge_all_vcfs.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        files=`echo {input.bcf} | awk '{{print NF}}'`
        if [ $files -gt 1 ] # check if there are at least 2 files for merging. If there is only one file, copy the bcf file.
        then
          bcftools merge -m snps -O b -o {output.merged} {input.bcf} 2> {log}
        else
          cp {input.bcf} {output.merged} && touch {output.merged} 2> {log}
          echo "Only one file present for merging. Copying the input bcf file." >> {log}
        fi
        """


rule index_merged_vcf:
    """Index vcf files"""
    input:
        bcf="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf",
    output:
        index=temp("results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf.csi"),
    group:
        "merged_vcf_group"
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_index_merged_vcfs.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools index -o {output.index} {input.bcf} 2> {log}
        """


rule merged_vcf_stats:
    """Obtain summary stats of merged vcf file before removing sites with missing data"""
    input:
        merged="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf",
        index=rules.index_merged_vcf.output,
    output:
        stats="results/all/vcf/" + REF_NAME + "/stats/vcf_merged/" + REF_NAME + ".all.merged.snps.bcf.stats.txt",
    group:
        "merged_vcf_group"
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_merged_vcf_stats.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools stats {input.merged} > {output.stats} 2> {log}
        """


rule merged_vcf_multiqc:
    """Collect all stats files from merged vcf files filtered for missing data"""
    input:
        rules.merged_vcf_stats.output.stats,
    output:
        stats="results/all/vcf/" + REF_NAME + "/stats/vcf_merged/multiqc/multiqc_report.html",
    params:
        indir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged/",
        outdir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged/multiqc",
    group:
        "merged_vcf_group"
    log:
        "results/logs/9_merge_vcfs/all/" + REF_NAME + "/merged_vcf_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule filter_vcf_biallelic:
    input:
        bcf="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf",
        index="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf.csi",
        stats=rules.merged_vcf_stats.output.stats,
        multiqc=rules.merged_vcf_multiqc.output.stats,
    output:
        bcf=temp("results/all/vcf/" + REF_NAME + ".all.merged.biallelic.bcf"),
        index=temp("results/all/vcf/" + REF_NAME + ".all.merged.biallelic.bcf.csi"),
    threads: 2
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_filter_vcf_biallelic.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools view -m2 -M2 -v snps -Ob -o {output.bcf} {input.bcf} 2> {log} &&
        bcftools index -f {output.bcf} 2>> {log}
        """


rule biallelic_filtered_vcf_stats:
    """Obtain summary stats of merged vcf file"""
    input:
        bcf=rules.filter_vcf_biallelic.output.bcf,
        index=rules.filter_vcf_biallelic.output.index,
    output:
        stats="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_biallelic/" + REF_NAME + ".all.merged.biallelic.vcf.stats.txt",
    group:
        "biallelic_filtered_vcf_group"
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_biallelic_filtered_vcf_stats.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools stats {input.bcf} > {output.stats} 2> {log}
        """


rule biallelic_filtered_vcf_multiqc:
    """Collect all stats files from merged vcf files filtered for biallelic sites"""
    input:
        rules.biallelic_filtered_vcf_stats.output.stats,
    output:
        stats="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_biallelic/multiqc/multiqc_report.html",
    params:
        indir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_biallelic/",
        outdir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_biallelic/multiqc",
    group:
        "biallelic_filtered_vcf_group"
    log:
        "results/logs/9_merge_vcfs/all/" + REF_NAME + "/biallelic_filtered_vcf_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule filter_vcf_missing:
    input:
        bcf=rules.filter_vcf_biallelic.output.bcf,
        index=rules.filter_vcf_biallelic.output.index,
        multiqc=rules.biallelic_filtered_vcf_multiqc.output.stats,
    output:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.genome.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.genome.vcf.gz.csi",
    threads: 2
    params:
        fmiss=config["f_missing"],
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_fmissing{fmiss}_filter_vcf_missing.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        # only include sites with zero missing data
        if [[ `echo 0.0 {params.fmiss} | awk '{{print ($1 == $2)}}'` == 1 ]]
        then
          bcftools view -i 'F_MISSING = {params.fmiss}' -Oz -o {output.vcf} {input.bcf} 2> {log}
        # include all sites
        elif [[ `echo 1.0 {params.fmiss} | awk '{{print ($1 == $2)}}'` == 1 ]]
        then
          bcftools view -i 'F_MISSING <= {params.fmiss}' -Oz -o {output.vcf} {input.bcf} 2> {log}
        # include sites with less than the fraction f_missing of missing data
        elif [[ `echo 0.0 {params.fmiss} 1.0 | awk '{{print ($1 < $2 && $2 < $3)}}'` == 1 ]]
        then 
          bcftools view -i 'F_MISSING < {params.fmiss}' -Oz -o {output.vcf} {input.bcf} 2> {log}
        fi
        
        bcftools index -f {output.vcf} 2>> {log}
        """


rule remove_chromosomes:
    input:
        bcf=rules.filter_vcf_missing.output.vcf,
        index=rules.filter_vcf_missing.output.index,
    output:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.autos.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.autos.vcf.gz.csi",
    threads: 2
    params:
        exclude = ",".join(sexchromosomeList) # parse list with contigs/scaffolds to exclude and convert to format chr1,chr2,chr3 for removal with bcftools view
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_fmissing{fmiss}.autos_remove_chromosomes.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools view {input.bcf} \
        -t ^{params.exclude} \
        -O z -o {output.vcf}

        bcftools index -f {output.vcf} 2>> {log}
        """


rule filtered_vcf2bed:
    """Convert the VCF file after removal of missing data (and optionally sex chromosomes) to BED file containing the remaining sites"""
    input:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    output:
        bed="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.bed",
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_fmissing{fmiss}.{chr}_filtered_vcf2bed.log",
    singularity:
        "docker://nbisweden/generode-bedtools-2.29.2"
    shell:
        """
        gzip -cd {input.vcf} | grep -v "^#" | awk -F'\t' '{{print $1, $2-1, $2}}' OFS='\t' > {output.bed} 2> {log}
        """


rule extract_historical_samples:
    input:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
        bed=rules.filtered_vcf2bed.output.bed,
    output:
        vcf="results/historical/vcf/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/historical/vcf/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".historical_fmissing{fmiss}.{chr}_extract_historical_samples.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    params:
        samples=hist_sm,
        all_samples=ALL_SAMPLES,
    shell:
        """
        samples_edited=`echo {params.samples} | sed 's/ /,/g'`
        samples_len=`echo {params.samples} | wc -w` # count the number of historical samples
        all_samples_len=`echo {params.all_samples} | wc -w` # count the number of all samples

        if [ $samples_len != $all_samples_len ]
        then
          bcftools view -Oz -s $samples_edited -o {output.vcf} {input.vcf} 2> {log} &&
          bcftools index -f {output.vcf} 2>> {log}
        else
          cp {input.vcf} {output.vcf} && touch {output.vcf} 2> {log} &&
          bcftools index -f {output.vcf} 2>> {log}
          echo "Only historical samples present. Copying the input vcf file." >> {log}
        fi
        """


rule extract_modern_samples:
    input:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
        bed=rules.filtered_vcf2bed.output.bed,
    output:
        vcf="results/modern/vcf/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/modern/vcf/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".modern_fmissing{fmiss}.{chr}_extract_modern_samples.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    params:
        samples=mod_sm,
        all_samples=ALL_SAMPLES,
    shell:
        """
        samples_edited=`echo {params.samples} | sed 's/ /,/g'`
        samples_len=`echo {params.samples} | wc -w` # count the number of historical samples
        all_samples_len=`echo {params.all_samples} | wc -w` # count the number of all samples

        if [ $samples_len != $all_samples_len ]
        then
          bcftools view -Oz -s $samples_edited -o {output.vcf} {input.vcf} 2> {log} &&
          bcftools index -f {output.vcf} 2>> {log}
        else
          cp {input.vcf} {output.vcf} && touch {output.vcf} 2> {log} &&
          bcftools index -f {output.vcf} 2>> {log}
          echo "Only modern samples present. Copying the input vcf file." >> {log}
        fi
        """


rule missingness_filtered_vcf_stats:
    """Obtain summary stats of merged vcf file"""
    input:
        merged="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    output:
        stats="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_merged_missing/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}_missingness_filtered_vcf_stats.log",
    singularity:
        "docker://quay.io/biocontainers/bcftools:1.9--h68d8f2e_9"
    shell:
        """
        bcftools stats {input.merged} > {output.stats} 2> {log}
        """


rule missingness_filtered_vcf_multiqc:
    """Collect all stats files from merged vcf files filtered for missing data"""
    input:
        missingness_filtered_vcf_multiqc_inputs,
    output:
        stats="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_missing/multiqc/multiqc_report.html",
    params:
        indir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_missing/",
        outdir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_missing/multiqc",
    log:
        "results/logs/9_merge_vcfs/all/" + REF_NAME + "/missingness_filtered_vcf_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """