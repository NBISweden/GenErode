##########################################################################
### 4. Genotyping

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    all_outputs.append("results/historical/vcf/" + REF_NAME + "/stats/vcf_sorted/multiqc/multiqc_report.html")

if os.path.exists(config["modern_samples"]):
    all_outputs.append("results/modern/vcf/" + REF_NAME + "/stats/vcf_sorted/multiqc/multiqc_report.html")


# snakemake rules
rule variant_calling:
    """Call variants in historical and modern samples (each sample on its own)"""
    """Input bam files can be just realn or additionally resca"""
    """Minimum mapping quality for a read to be considered: 30"""
    """Minimum base quality for a base to be considered: 30"""
    """-B: Disabled probabilistic realignment for the computation of base alignment quality (BAQ). BAQ is the Phred-scaled probability of a read base being misaligned. Applying this option greatly helps to reduce false SNPs caused by misalignments"""
    input:
        ref=config["ref_path"],
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.bam",
    output:
        bcf=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.Q30.bcf"),
    resources:
        cpus_per_task=3,
    log:
        "results/logs/4_genotyping/{dataset}/" + REF_NAME + "/{sample}.{processed}_variant_calling.log",
    singularity:
        "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0"
    shell:
        """
        bcftools mpileup -Ou -Q 30 -q 30 -B -f {input.ref} {input.bam} | bcftools call -c -M -O b --threads {resources.cpus_per_task} -o {output.bcf} 2> {log}
        """


rule sort_vcfs:
    """Sort vcf files before any downstream processing"""
    input:
        bcf=rules.variant_calling.output.bcf,
    output:
        sort="results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.Q30.sorted.bcf",
    resources:
        cpus_per_task=2,
    params:
        tmpdir="results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.Q30.sorted_XXXXXX/",
    log:
        "results/logs/4_genotyping/{dataset}/" + REF_NAME + "/{sample}.{processed}_sort_vcfs.log",
    singularity:
        "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0"
    shell:
        """
        bcftools sort -O b -T {params.tmpdir} -o {output.sort} {input.bcf} 2> {log}
        """


rule index_sorted_vcfs:
    """Index sorted and fixed vcf files before any downstream processing"""
    input:
        sort=rules.sort_vcfs.output.sort,
    output:
        index="results/{dataset}/vcf/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.Q30.sorted.bcf.csi",
    group:
        "sorted_vcf_group"
    log:
        "results/logs/4_genotyping/{dataset}/" + REF_NAME + "/{sample}.{processed}_index_sorted_vcfs.log",
    singularity:
        "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0"
    shell:
        """
        bcftools index -o {output.index} {input.sort} 2> {log}
        """


rule sorted_vcf_stats:
    """Obtain summary stats of vcf files"""
    input:
        sort=rules.sort_vcfs.output.sort,
        index=rules.index_sorted_vcfs.output.index,
    output:
        stats="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_sorted/{sample}.merged.rmdup.merged.{processed}.Q30.sorted.vcf.stats.txt",
    group:
        "sorted_vcf_group"
    log:
        "results/logs/4_genotyping/{dataset}/" + REF_NAME + "/{sample}.{processed}_sorted_vcf_stats.log",
    singularity:
        "https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0"
    shell:
        """
        bcftools stats {input.sort} > {output.stats} 2> {log}
        """


rule historical_sorted_vcf_multiqc:
    """Collect all stats files from historical vcf files (sorted)"""
    input:
        not_rescaled_not_subsampled=expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_sorted/{sample}.merged.rmdup.merged.realn.Q30.sorted.vcf.stats.txt",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_SAMPLES,),
        rescaled_not_subsampled=expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_sorted/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.vcf.stats.txt",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_SAMPLES,),
        not_rescaled_subsampled=expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_sorted/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.vcf.stats.txt",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_SAMPLES,
            DP=config["subsampling_depth"],),
        rescaled_subsampled=expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_sorted/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.vcf.stats.txt",
            sample=HIST_RESCALED_SUBSAMPLED_SAMPLES,
            DP=config["subsampling_depth"],),
    output:
        stats="results/historical/vcf/" + REF_NAME + "/stats/vcf_sorted/multiqc/multiqc_report.html",
    params:
        indir="results/historical/vcf/" + REF_NAME + "/stats/vcf_sorted/",
        outdir="results/historical/vcf/" + REF_NAME + "/stats/vcf_sorted/multiqc/",
    log:
        "results/logs/4_genotyping/historical/" + REF_NAME + "/historical_sorted_vcf_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_sorted_vcf_multiqc:
    """Collect all stats files from modern vcf files (sorted)"""
    input:
        not_subsampled=expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_sorted/{sample}.merged.rmdup.merged.realn.Q30.sorted.vcf.stats.txt",
            sample=MODERN_NOT_SUBSAMPLED_SAMPLES,),
        subsampled=expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_sorted/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.vcf.stats.txt",
            sample=MODERN_SUBSAMPLED_SAMPLES,
            DP=config["subsampling_depth"],),
    output:
        stats="results/modern/vcf/" + REF_NAME + "/stats/vcf_sorted/multiqc/multiqc_report.html",
    params:
        indir="results/modern/vcf/" + REF_NAME + "/stats/vcf_sorted/",
        outdir="results/modern/vcf/" + REF_NAME + "/stats/vcf_sorted/multiqc/",
    log:
        "results/logs/4_genotyping/modern/" + REF_NAME + "/modern_sorted_vcf_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
