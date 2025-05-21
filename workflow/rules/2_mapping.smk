##########################################################################
### 2. Mapping

# Code collecting output files from this part of the pipeline
mapping_outputs=[]
if os.path.exists(config["historical_samples"]):
    mapping_outputs.append("results/historical/mapping/" + REF_NAME + "/stats/bams_sorted/multiqc/multiqc_report.html")

if os.path.exists(config["modern_samples"]):
    mapping_outputs.append("results/modern/mapping/" + REF_NAME + "/stats/bams_sorted/multiqc/multiqc_report.html")


localrules:
    readgroup_ID_historical,
    readgroup_ID_modern,


# snakemake rules
rule map_historical:
    """Map trimmed and merged reads from historical samples to reference. BWA aln for short Illumina reads, parameters according to Palkopoulou et al. 2015"""
    input:
        ref=config["ref_path"],
        index=rules.bwa_index_reference.output,
        fastq_hist=rules.fastp_historical.output.merged,
    output:
        sai=temp("results/historical/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sai"),
    log:
        "results/logs/2_mapping/historical/" + REF_NAME + "/{sample}_{index}_{lane}_map_historical.log",
    threads: 8
    singularity:
        bwa_container
    shell:
        """
        bwa aln -l 16500 -n 0.01 -o 2 -t {threads} {input.ref} {input.fastq_hist} > {output.sai} 2> {log}
        """


rule readgroup_ID_historical:
    """Generate read group file"""
    input:
        sai=rules.map_historical.output.sai,
    output:
        rg=temp("results/historical/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.readgroup.txt"),
    log:
        "results/logs/2_mapping/historical/" + REF_NAME + "/{sample}_{index}_{lane}_readgroup_ID_historical.log",
    run:
        FULLSAMPLENAME = "{}_{}_{}".format(wildcards.sample, wildcards.index, wildcards.lane)
        shell("""
            echo '@RG\\tID:{ID}\\tSM:{SM}\\tPL:{PL}\\tLB:{LB}' > {{output.rg}} 2> {{log}}
            """.format(FULLSAMPLENAME=FULLSAMPLENAME, **hist_rg_dict[FULLSAMPLENAME]))


rule sai2bam:
    """Convert SAI files to BAM files and add read groups"""
    input:
        ref=config["ref_path"],
        index=rules.bwa_index_reference.output,
        fastq_hist=rules.fastp_historical.output.merged,
        sai=rules.map_historical.output.sai,
        rg=rules.readgroup_ID_historical.output,
    output:
        bam="results/historical/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sorted.bam",
    log:
        "results/logs/2_mapping/historical/" + REF_NAME + "/{sample}_{index}_{lane}_sai2bam.log",
    threads: 8
    singularity:
        bwa_samtools_container
    shell:
        """
        bwa samse -r $(cat {input.rg}) {input.ref} {input.sai} {input.fastq_hist} | \
        samtools sort -@ {threads} - > {output.bam} 2> {log}
        """


rule readgroup_ID_modern:
    """Generate read group file"""
    input:
        ref=config["ref_path"],
        index=rules.bwa_index_reference.output,
        fastq_mod_R1=rules.fastp_modern.output.R1_trimmed,
        fastq_mod_R2=rules.fastp_modern.output.R2_trimmed,
    output:
        rg=temp("results/modern/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.readgroup.txt"),
    log:
        "results/logs/2_mapping/modern/" + REF_NAME + "/{sample}_{index}_{lane}_readgroup_ID_modern.log",
    run:
        FULLSAMPLENAME = "{}_{}_{}".format(wildcards.sample, wildcards.index, wildcards.lane)
        shell("""
            echo '@RG\\tID:{ID}\\tSM:{SM}\\tPL:{PL}\\tLB:{LB}' > {{output.rg}} 2> {{log}}
            """.format(FULLSAMPLENAME=FULLSAMPLENAME, **mod_rg_dict[FULLSAMPLENAME]))


rule map_modern:
    """
    Map trimmed reads from modern samples to reference using BWA mem for long Illumina reads.
    Shorter split hits are marked as secondary for Picard.
    """
    input:
        ref=config["ref_path"],
        index=rules.bwa_index_reference.output,
        fastq_mod_R1=rules.fastp_modern.output.R1_trimmed,
        fastq_mod_R2=rules.fastp_modern.output.R2_trimmed,
        rg=rules.readgroup_ID_modern.output,
    output:
        bam="results/modern/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sorted.bam",
    log:
        "results/logs/2_mapping/modern/" + REF_NAME + "/{sample}_{index}_{lane}_map_modern.log",
    threads: 8
    singularity:
        bwa_samtools_container
    shell:
        """
        bwa mem -M -t {threads} -R $(cat {input.rg}) {input.ref} {input.fastq_mod_R1} {input.fastq_mod_R2} | \
        samtools sort -@ {threads} - > {output.bam} 2> {log}
        """


rule index_sorted_bams:
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sorted.bam",
    output:
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sorted.bam.bai",
    log:
        "results/logs/2_mapping/{dataset}/" + REF_NAME + "/{sample}_{index}_{lane}_index_sorted_bams.log",
    singularity:
        bwa_samtools_container
    shell:
        """
        samtools index {input.bam} {output.index} 2> {log}
        """


rule sorted_bam_stats:
    """Basic statistics on mapping output"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sorted.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sorted.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_sorted/{sample}_{index}_{lane}.sorted.bam.stats.txt",
    log:
        "results/logs/2_mapping/{dataset}/" + REF_NAME + "/{sample}_{index}_{lane}_sorted_bam_stats.log",
    singularity:
        bwa_samtools_container
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule sorted_bam_qualimap:
    """More detailed stats"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sorted.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}_{lane}.sorted.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_sorted/{sample}_{index}_{lane}.sorted.bam.qualimap/qualimapReport.html",
        results="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_sorted/{sample}_{index}_{lane}.sorted.bam.qualimap/genome_results.txt",
        outdir=directory("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_sorted/{sample}_{index}_{lane}.sorted.bam.qualimap/"),
    threads: 8
    resources:
        mem_mb=64000,
    params:
        outdir="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_sorted/{sample}_{index}_{lane}.sorted.bam.qualimap",
    log:
        "results/logs/2_mapping/{dataset}/" + REF_NAME + "/{sample}_{index}_{lane}_sorted_bam_qualimap.log",
    singularity:
        qualimap_container
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G -nt {threads} -outdir {params.outdir} -outformat html 2> {log}
        """


rule historical_raw_bam_multiqc:
    """Summarize all stats results from all historical bam files after mapping"""
    input:
        stats=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_sorted/{sampleindexlane}.sorted.bam.stats.txt",
            sampleindexlane=hist_pipeline_bam_sm_idx_ln,),
        qualimap=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_sorted/{sampleindexlane}.sorted.bam.qualimap/qualimapReport.html",
            sampleindexlane=hist_pipeline_bam_sm_idx_ln,),
    output:
        "results/historical/mapping/"+ REF_NAME + "/stats/bams_sorted/multiqc/multiqc_report.html",
    params:
        indir="results/historical/mapping/" + REF_NAME + "/stats/bams_sorted/",
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_sorted/multiqc/",
    log:
        "results/logs/2_mapping/historical/" + REF_NAME + "/historical_raw_bam_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_raw_bam_multiqc:
    """Summarize all stats results from all modern bam files after mapping"""
    input:
        stats=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_sorted/{sampleindexlane}.sorted.bam.stats.txt",
            sampleindexlane=mod_pipeline_bam_sm_idx_ln,),
        qualimap=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_sorted/{sampleindexlane}.sorted.bam.qualimap/qualimapReport.html",
            sampleindexlane=mod_pipeline_bam_sm_idx_ln,),
    output:
        "results/modern/mapping/" + REF_NAME + "/stats/bams_sorted/multiqc/multiqc_report.html",
    params:
        indir="results/modern/mapping/" + REF_NAME + "/stats/bams_sorted/",
        outdir="results/modern/mapping/" + REF_NAME + "/stats/bams_sorted/multiqc/",
    log:
        "results/logs/2_mapping/modern/" + REF_NAME + "/modern_raw_bam_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
