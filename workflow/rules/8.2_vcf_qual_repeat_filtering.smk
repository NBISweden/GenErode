##########################################################################
### 8.2 Quality filtering and repeat filtering of VCF files per sample

# Code collecting output files from this part of the pipeline
vcf_proc_outputs=[]

if os.path.exists(config["historical_samples"]):
    vcf_proc_outputs.append("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/multiqc/multiqc_report.html")
    vcf_proc_outputs.append("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/multiqc/multiqc_report.html")

if os.path.exists(config["modern_samples"]):
    vcf_proc_outputs.append("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/multiqc/multiqc_report.html")
    vcf_proc_outputs.append("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/multiqc/multiqc_report.html")


# Functions used by rules of this part of the pipeline
def historical_quality_filtered_vcf_multiqc_inputs(wildcards):
    """Input for historical_quality_filtered_vcf_multiqc"""
    outlist = []
    outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
        sample=HIST_NOT_CpG_SAMPLES,)
    if config["CpG_from_vcf"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_CpG_SAMPLES,)
    elif config["CpG_from_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_CpG_SAMPLES,)
    elif config["CpG_from_vcf_and_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=HIST_CpG_SAMPLES,)
    return outlist

def modern_quality_filtered_vcf_multiqc_inputs(wildcards):
    """Input for modern_quality_filtered_vcf_multiqc"""
    outlist = []
    outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
        sample=MOD_NOT_CpG_SAMPLES,)
    if config["CpG_from_vcf"] == True:
        outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_reference"] == True:
        outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_vcf_and_reference"] == True:
        outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
            sample=MOD_CpG_SAMPLES,)
    return outlist

def historical_repmasked_vcf_multiqc_inputs(wildcards):
    """Input for historical_repmasked_vcf_multiqc"""
    outlist = []
    outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
        sample=HIST_NOT_CpG_SAMPLES,)
    if config["CpG_from_vcf"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_CpG_SAMPLES,)
    elif config["CpG_from_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_CpG_SAMPLES,)
    elif config["CpG_from_vcf_and_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=HIST_CpG_SAMPLES,)
    return outlist

def modern_repmasked_vcf_multiqc_inputs(wildcards):
    """Input for modern_repmasked_vcf_multiqc"""
    outlist = []
    outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
        sample=MOD_NOT_CpG_SAMPLES,)
    if config["CpG_from_vcf"] == True:
        outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_reference"] == True:
        outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_vcf_and_reference"] == True:
        outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
            sample=MOD_CpG_SAMPLES,)
    return outlist


# snakemake rules
rule remove_snps_near_indels:
    """remove SNPs within 5 bp of an indel"""
    input:
        bcf="results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.bcf",
    output:
        snps=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.bcf"),
    threads: 2
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_remove_snps_near_indels.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools filter -g 5 -O b --threads {threads} -o {output.snps} {input.bcf} 2> {log}
        """


rule filter_vcfs_qual_dp:
    """
    Remove indels, genotypes of genotype quality < 30 and keep only sites within depth thresholds 
    that were determined from bam files earlier in the pipeline.
    Note that the depth filter is recalculated for subsampled bam files, according to the target depth for subsampling.
    """
    input:
        bcf=rules.remove_snps_near_indels.output,
        dp=depth_file,
    output:
        filtered=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.bcf"),
    threads: 2
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_filter_vcfs_qual_dp.log",
    singularity:
        bcftools_container
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
    """
    Removes heterozygote sites with allelic imbalance from the vcf files that could be due to contamination, sequencing or mapping errors.
    For example, sites where 9 reads support the reference allele and one read the alternative allele.
    """
    input:
        bcf=rules.filter_vcfs_qual_dp.output.filtered,
    output:
        filtered=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.bcf"),
    threads: 2
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_filter_vcfs_allelic_balance.log",
    singularity:
        bcftools_container
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
        index=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.bcf.csi"),
    group:
        "qual_filtered_vcf_group"
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_index_filtered_vcfs.log",
    singularity:
        bcftools_container
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
        stats="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_qual_filtered/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.bcf.stats.txt",
    group:
        "qual_filtered_vcf_group"
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_filtered_vcf_stats.log",
    singularity:
        bcftools_container
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
        multiqc_container
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
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule filtered_bcf2vcf:
    """Convert bcf format to vcf.gz for removal of sites"""
    input:
        bcf=rules.filter_vcfs_allelic_balance.output.filtered,
    output:
        vcf=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.vcf.gz"),
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_filtered_bcf2vcf.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools convert -O z -o {output.vcf} {input.bcf} 2> {log}
        """


rule remove_repeats_vcf:
    """Remove repeats from vcf files"""
    input:
        vcf=rules.filtered_bcf2vcf.output.vcf,
        bed=rules.make_no_repeats_bed.output.no_rep_bed,
        genomefile=rules.genome_file.output.genomefile,
    output:
        filtered=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.vcf.gz"),
    threads: 6
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_remove_repeats_vcf.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools intersect -a {input.vcf} -b {input.bed} -header -sorted -g {input.genomefile} | bgzip -c > {output.filtered} 2> {log}
        """


rule filtered_vcf2bcf:
    """Convert the repeat masked vcf back to bcf"""
    input:
        filtered=rules.remove_repeats_vcf.output.filtered,
    output:
        bcf="results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
    threads: 2
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_filtered_vcf2bcf.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools convert -O b -o {output.bcf} {input.filtered} 2> {log}
        """


rule index_repmasked_vcfs:
    """Index vcf files before any downstream processing"""
    input:
        bcf=rules.filtered_vcf2bcf.output.bcf,
    output:
        index="results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
    group:
        "repmasked_vcf_group"
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_index_repmasked_vcfs.log",
    singularity:
        bcftools_container
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
        stats="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_repmasked/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.bcf.stats.txt",
    group:
        "repmasked_vcf_group"
    log:
        "results/logs/8.2_vcf_qual_repeat_filtering/{dataset}/" + REF_NAME + "/{sample}.{filtered}_repmasked_vcf_stats.log",
    singularity:
        bcftools_container
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
        multiqc_container
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
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
