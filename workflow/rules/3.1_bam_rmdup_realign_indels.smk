##########################################################################
### 3.1 BAM file processing: sort, merge samples from different lanes per PCR/index, remove duplicates, merge bam files per sample, realign indels, calculate depth

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    all_outputs.append("results/historical/mapping/" + REF_NAME + "/stats/bams_merged_index/multiqc/multiqc_report.html")
    all_outputs.append("results/historical/mapping/" + REF_NAME + "/stats/bams_rmdup/multiqc/multiqc_report.html")
    all_outputs.append("results/historical/mapping/" + REF_NAME + "/stats/bams_merged_sample/multiqc/multiqc_report.html")
    all_outputs.append("results/historical/mapping/" + REF_NAME + "/stats/bams_indels_realigned/multiqc/multiqc_report.html")

if os.path.exists(config["modern_samples"]):
    all_outputs.append("results/modern/mapping/" + REF_NAME + "/stats/bams_merged_index/multiqc/multiqc_report.html")
    all_outputs.append("results/modern/mapping/" + REF_NAME + "/stats/bams_rmdup/multiqc/multiqc_report.html")
    all_outputs.append("results/modern/mapping/" + REF_NAME + "/stats/bams_merged_sample/multiqc/multiqc_report.html")
    all_outputs.append("results/modern/mapping/" + REF_NAME + "/stats/bams_indels_realigned/multiqc/multiqc_report.html")


# Functions used by rules of this part of the pipeline
def merge_hist_bams_per_index_inputs(wildcards):
    """Input for merging historical bam files per index"""
    SAMPLEIDX = "{}_{}".format(wildcards.sample, wildcards.index)
    SAMPLEIDXLN_LIST = hist_sampleidxln_dict[SAMPLEIDX]
    return expand("results/historical/mapping/" + REF_NAME + "/{sampleindexlane}.sorted.bam",
        sampleindexlane=SAMPLEIDXLN_LIST,)

def merge_mod_bams_per_index_inputs(wildcards):
    """Input for merging modern bam files per index"""
    SAMPLEIDX = "{}_{}".format(wildcards.sample, wildcards.index)
    SAMPLEIDXLN_LIST = mod_sampleidxln_dict[SAMPLEIDX]
    return expand("results/modern/mapping/" + REF_NAME + "/{sampleindexlane}.sorted.bam",
        sampleindexlane=SAMPLEIDXLN_LIST,)

def merge_hist_bams_per_sample_inputs(wildcards):
    """Input for merging historical bam files per sample"""
    SAMPLEIDX_LIST = hist_sampleidx_dict[wildcards.sample]
    return expand("results/historical/mapping/" + REF_NAME + "/{sampleindex}.merged.rmdup.bam",
        sampleindex=SAMPLEIDX_LIST,)

def merge_mod_bams_per_sample_inputs(wildcards):
    """Input for merging modern bam files per sample"""
    SAMPLEIDX_LIST = mod_sampleidx_dict[wildcards.sample]
    return expand("results/modern/mapping/" + REF_NAME + "/{sampleindex}.merged.rmdup.bam",
        sampleindex=SAMPLEIDX_LIST,)


# snakemake rules
rule merge_historical_bams_per_index:
    """Merge files per sample and index in case some were sequenced on different lanes"""
    input:
        merge_hist_bams_per_index_inputs,
    output:
        merged=temp("results/historical/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam"),
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/historical/" + REF_NAME + "/{sample}_{index}_merge_historical_bams_per_index.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        files=`echo {input} | awk '{{print NF}}'`
        if [ $files -gt 1 ] # check if there are at least 2 files for merging. If there is only one file, copy the sorted bam file.
        then
            samtools merge {output.merged} {input} 2> {log}
        else
            cp {input} {output.merged} && touch {output.merged} 2> {log}
            echo "Only one file present for merging. Copying the input bam file." >> {log}
        fi
        """


rule merge_modern_bams_per_index:
    """Merge files per sample and index in case some were sequenced on different lanes"""
    input:
        merge_mod_bams_per_index_inputs,
    output:
        merged=temp("results/modern/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam"),
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/modern/" + REF_NAME + "/{sample}_{index}_merge_modern_bams_per_index.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        files=`echo {input} | awk '{{print NF}}'`
        if [ $files -gt 1 ] # check if there are at least 2 files for merging. If there is only one file, copy the sorted bam file.
        then
            samtools merge {output.merged} {input} 2> {log}
        else
            cp {input} {output.merged} && touch {output.merged} 2> {log}
            echo "Only one file present for merging. Copying the input bam file." >> {log}
        fi
        """


rule index_merged_index_bams:
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam",
    output:
        index=temp("results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam.bai"),
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_{index}_index_merged_index_bams.log",
    group:
        "merged_index_bam_group"
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools index {input.bam} {output.index} 2> {log}
        """


rule merged_index_bam_stats:
    """Basic statistics on mapping output"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_index/{sample}_{index}.merged.bam.stats.txt",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_{index}_merged_index_bam_stats.log",
    group:
        "merged_index_bam_group"
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule merged_index_bam_qualimap:
    """More detailed stats"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_index/{sample}_{index}.merged.bam.qualimap/qualimapReport.html",
        results="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_index/{sample}_{index}.merged.bam.qualimap/genome_results.txt",
        outdir=directory("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_index/{sample}_{index}.merged.bam.qualimap"),
    threads: 8
    resources:
        mem_mb=64000,
    params:
        outdir="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_index/{sample}_{index}.merged.bam.qualimap",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_{index}_merged_index_bam_qualimap.log",
    singularity:
        "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G -nt {threads} -outdir {params.outdir} -outformat html 2> {log}
        """


rule historical_merged_index_bam_multiqc:
    """Summarize all stats from all historical bam files"""
    input:
        merged=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_merged_index/{sampleindex}.merged.bam.stats.txt",
            sampleindex=hist_pipeline_bam_sm_idx,),
        qualimap=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_merged_index/{sampleindex}.merged.bam.qualimap/qualimapReport.html",
            sampleindex=hist_pipeline_bam_sm_idx,),
    output:
        "results/historical/mapping/" + REF_NAME + "/stats/bams_merged_index/multiqc/multiqc_report.html",
    params:
        indir="results/historical/mapping/" + REF_NAME + "/stats/bams_merged_index/",
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_merged_index/multiqc/",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/historical/" + REF_NAME + "/historical_merged_index_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_merged_index_bam_multiqc:
    """Summarize all stats from all modern bam files"""
    input:
        merged=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_merged_index/{sampleindex}.merged.bam.stats.txt",
            sampleindex=mod_pipeline_bam_sm_idx,),
        qualimap=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_merged_index/{sampleindex}.merged.bam.qualimap/qualimapReport.html",
            sampleindex=mod_pipeline_bam_sm_idx,),
    output:
        "results/modern/mapping/" + REF_NAME + "/stats/bams_merged_index/multiqc/multiqc_report.html",
    params:
        indir="results/modern/mapping/" + REF_NAME + "/stats/bams_merged_index/",
        outdir="results/modern/mapping/" + REF_NAME + "/stats/bams_merged_index/multiqc/",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/modern/" + REF_NAME + "/modern_merged_index_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule rmdup_historical_bams:
    """Remove PCR duplicates from historical samples using Pontus Skoglund's custom script for duplicate removal (checks both ends of a read)"""
    """The script was modified so that also unmapped reads are printed to the output bam file so that they are not lost"""
    input:
        merged=rules.merge_historical_bams_per_index.output.merged,
        index="results/historical/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam.bai",
    output:
        rmdup=temp("results/historical/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam"),
    threads: 6
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/historical/" + REF_NAME + "/{sample}_{index}_rmdup_historical_bams.log",
    singularity:
        "oras://community.wave.seqera.io/library/samtools_python:2e56d0f345426c81"
    shell:
        """
        samtools view -@ {threads} -h {input.merged} | python3 workflow/scripts/samremovedup.py | \
        samtools view -b -o {output.rmdup} 2> {log}
        """


rule rmdup_modern_bams:
    """Mark duplicates in modern samples"""
    input:
        merged=rules.merge_modern_bams_per_index.output.merged,
        index="results/modern/mapping/" + REF_NAME + "/{sample}_{index}.merged.bam.bai",
    output:
        rmdup=temp("results/modern/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam"),
        metrix=temp("results/modern/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup_metrics.txt"),
    threads: 2
    resources:
        mem_mb=16000,
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/modern/" + REF_NAME + "/{sample}_{index}_rmdup_modern_bams.log",
    singularity:
        "docker://quay.io/biocontainers/picard:2.26.6--hdfd78af_0"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        picard MarkDuplicates -Xmx${{mem}}g INPUT={input.merged} OUTPUT={output.rmdup} METRICS_FILE={output.metrix} 2> {log}
        """


rule index_rmdup_bams:
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam",
    output:
        index=temp("results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam.bai"),
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_{index}_index_rmdup_bams.log",
    group:
        "rmdup_bam_group"
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools index {input.bam} {output.index} 2> {log}
        """


rule rmdup_bam_stats:
    """Basic statistics on mapping output"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_rmdup/{sample}_{index}.merged.rmdup.bam.stats.txt",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_{index}_rmdup_bam_stats.log",
    group:
        "rmdup_bam_group"
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule rmdup_bam_qualimap:
    """More detailed stats"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}_{index}.merged.rmdup.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_rmdup/{sample}_{index}.merged.rmdup.bam.qualimap/qualimapReport.html",
        results="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_rmdup/{sample}_{index}.merged.rmdup.bam.qualimap/genome_results.txt",
        outdir=directory("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_rmdup/{sample}_{index}.merged.rmdup.bam.qualimap"),
    threads: 8
    resources:
        mem_mb=64000,
    params:
        outdir="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_rmdup/{sample}_{index}.merged.rmdup.bam.qualimap",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_{index}_rmdup_bam_qualimap.log",
    singularity:
        "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G -nt {threads} -outdir {params.outdir} -outformat html 2> {log}
        """


rule historical_rmdup_bam_multiqc:
    """Summarize all stats from all historical bam files"""
    input:
        rmdup=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_rmdup/{sampleindex}.merged.rmdup.bam.stats.txt",
            sampleindex=hist_pipeline_bam_sm_idx,),
        qualimap=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_rmdup/{sampleindex}.merged.rmdup.bam.qualimap/qualimapReport.html",
            sampleindex=hist_pipeline_bam_sm_idx,),
    output:
        "results/historical/mapping/" + REF_NAME + "/stats/bams_rmdup/multiqc/multiqc_report.html",
    params:
        indir="results/historical/mapping/" + REF_NAME + "/stats/bams_rmdup/",
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_rmdup/multiqc/",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/historical/" + REF_NAME + "/historical_rmdup_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_rmdup_bam_multiqc:
    """Summarize all stats from all modern bam files"""
    input:
        rmdup=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_rmdup/{sampleindex}.merged.rmdup.bam.stats.txt",
            sampleindex=mod_pipeline_bam_sm_idx,),
        qualimap=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_rmdup/{sampleindex}.merged.rmdup.bam.qualimap/qualimapReport.html",
            sampleindex=mod_pipeline_bam_sm_idx,),
    output:
        "results/modern/mapping/" + REF_NAME + "/stats/bams_rmdup/multiqc/multiqc_report.html",
    params:
        indir="results/modern/mapping/" + REF_NAME + "/stats/bams_rmdup/",
        outdir="results/modern/mapping/" + REF_NAME + "/stats/bams_rmdup/multiqc/",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/modern/" + REF_NAME + "/modern_rmdup_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule merge_historical_bams_per_sample:
    """Merge files per sample"""
    input:
        merge_hist_bams_per_sample_inputs,
    output:
        merged=temp("results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam"),
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/historical/" + REF_NAME + "/{sample}_merge_historical_bams_per_sample.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        files=`echo {input} | awk '{{print NF}}'`
        if [ $files -gt 1 ] # check if there are at least 2 files for merging. If there is only one file, copy the sorted bam file.
        then
            samtools merge {output.merged} {input} 2> {log}
        else
            cp {input} {output.merged} && touch {output.merged} 2> {log}
            echo "Only one file present for merging. Copying the input bam file." >> {log}
        fi
        """


rule merge_modern_bams_per_sample:
    """Merge files per sample"""
    input:
        merge_mod_bams_per_sample_inputs,
    output:
        merged=temp("results/modern/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam"),
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/modern/" + REF_NAME + "/{sample}_merge_modern_bams_per_sample.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        files=`echo {input} | awk '{{print NF}}'`
        if [ $files -gt 1 ] # check if there are at least 2 files for merging. If there is only one file, copy the sorted bam file.
        then
            samtools merge {output.merged} {input} 2> {log}
        else
            cp {input} {output.merged} && touch {output.merged} 2> {log}
            echo "Only one file present for merging. Copying the input bam file." >> {log}
        fi
        """


rule index_merged_sample_bams:
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam",
    output:
        index=temp("results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam.bai"),
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_index_merged_sample_bams.log",
    group:
        "merged_sample_bam_group"
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools index {input.bam} {output.index} 2> {log}
        """


rule merged_sample_bam_stats:
    """Basic statistics on mapping output"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_sample/{sample}.merged.rmdup.merged.bam.stats.txt",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_merged_sample_bam_stats.log",
    group:
        "merged_sample_bam_group"
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule merged_sample_bam_qualimap:
    """More detailed stats"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_sample/{sample}.merged.rmdup.merged.bam.qualimap/qualimapReport.html",
        results="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_sample/{sample}.merged.rmdup.merged.bam.qualimap/genome_results.txt",
        outdir=directory("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_sample/{sample}.merged.rmdup.merged.bam.qualimap"),
    threads: 8
    resources:
        mem_mb=64000,
    params:
        outdir="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_merged_sample/{sample}.merged.rmdup.merged.bam.qualimap",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_merged_sample_bam_qualimap.log",
    singularity:
        "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G -nt {threads} -outdir {params.outdir} -outformat html 2> {log}
        """


rule historical_merged_sample_bam_multiqc:
    """Summarize all stats from all historical bam files"""
    input:
        merged=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_merged_sample/{sample}.merged.rmdup.merged.bam.stats.txt",
            sample=hist_pipeline_bam_sm,),
        qualimap=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_merged_sample/{sample}.merged.rmdup.merged.bam.qualimap/qualimapReport.html",
            sample=hist_pipeline_bam_sm,),
    output:
        "results/historical/mapping/" + REF_NAME + "/stats/bams_merged_sample/multiqc/multiqc_report.html",
    params:
        indir="results/historical/mapping/" + REF_NAME + "/stats/bams_merged_sample/",
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_merged_sample/multiqc/",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/historical/" + REF_NAME + "/historical_merged_sample_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_merged_sample_bam_multiqc:
    """Summarize all stats from all modern bam files"""
    input:
        merged=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_merged_sample/{sample}.merged.rmdup.merged.bam.stats.txt",
            sample=mod_pipeline_bam_sm,),
        qualimap=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_merged_sample/{sample}.merged.rmdup.merged.bam.qualimap/qualimapReport.html",
            sample=mod_pipeline_bam_sm,),
    output:
        "results/modern/mapping/" + REF_NAME + "/stats/bams_merged_sample/multiqc/multiqc_report.html",
    params:
        indir="results/modern/mapping/" + REF_NAME + "/stats/bams_merged_sample/",
        outdir="results/modern/mapping/" + REF_NAME + "/stats/bams_merged_sample/multiqc/",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/modern/" + REF_NAME + "/modern_merged_sample_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule indel_realigner_targets:
    """Create a realignment target list to realign indels to improve alignments to reference and reduce the amount of false positive SNPs"""
    input:
        ref=config["ref_path"],
        fdict=REF_DIR + "/" + REF_NAME + ".dict",
        fai=config["ref_path"] + ".fai",
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam",
        bai="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam.bai",
    output:
        target_list=temp("results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn_targets.list"),
    threads: 8
    resources:
        mem_mb=64000,
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_indel_realigner_targets.log",
    singularity:
        "docker://broadinstitute/gatk3:3.7-0"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        java -jar -Xmx${{mem}}g /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator -R {input.ref} -I {input.bam} -o {output.target_list} -nt {threads} 2> {log}
        """


rule indel_realigner:
    """Realign indels to improve alignments to reference and reduce the amount of false positive SNPs"""
    input:
        ref=config["ref_path"],
        fdict=REF_DIR + "/" + REF_NAME + ".dict",
        fai=config["ref_path"] + ".fai",
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam",
        bai="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.bam.bai",
        target_list=rules.indel_realigner_targets.output,
    output:
        realigned="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam",
    threads: 8
    resources:
        mem_mb=64000,
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_indel_realigner.log",
    singularity:
        "docker://broadinstitute/gatk3:3.7-0"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        java -jar -Xmx${{mem}}g /usr/GenomeAnalysisTK.jar -T IndelRealigner -R {input.ref} -I {input.bam} -targetIntervals {input.target_list} -o {output.realigned} 2> {log}
        """


rule index_realigned_bams:
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam",
    output:
        index=temp("results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam.bai"),
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_index_realigned_bams.log",
    group:
        "realigned_bam_group"
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools index {input.bam} {output.index} 2> {log}
        """


rule realigned_bam_stats:
    """Basic statistics on mapping output"""
    input:
        bam=rules.indel_realigner.output.realigned,
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.bam.stats.txt",
    group:
        "realigned_bam_group"
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_realigned_bam_stats.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule realigned_bam_fastqc:
    """Run fastqc on all realigned bam files"""
    input:
        bam=rules.indel_realigner.output.realigned,
    output:
        html="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/fastqc/{sample}.merged.rmdup.merged.realn_fastqc.html",
        zip="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/fastqc/{sample}.merged.rmdup.merged.realn_fastqc.zip",
        dir=directory("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/fastqc/{sample}.merged.rmdup.merged.realn_fastqc"),
    params:
        dir="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/fastqc",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_realigned_bam_fastqc.log",
    threads: 2
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    shell:
        """
        fastqc -o {params.dir} -t {threads} --extract {input.bam} 2> {log}
        """


rule realigned_bam_qualimap:
    """More detailed stats"""
    input:
        bam=rules.indel_realigner.output.realigned,
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.bam.qualimap/qualimapReport.html",
        results="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.bam.qualimap/genome_results.txt",
        outdir=directory("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.bam.qualimap"),
    threads: 8
    resources:
        mem_mb=64000,
    params:
        outdir="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.bam.qualimap",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_realigned_bam_qualimap.log",
    singularity:
        "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G -nt {threads} -outdir {params.outdir} -outformat html 2> {log}
        """


rule realigned_bam_depth:
    """Get average genome-wide depth per bam file, excluding repeats and filtering for high quality"""
    """Will be used to filter out sites outside the estimated depth thresholds"""
    """samtools depth to get depth per site, piped into awk to get average across all sites, 
    piped into awk to calculate 1/3*average and 10*average depth (rounded to integer)"""
    input:
        bam=rules.indel_realigner.output.realigned,
        no_rep_bed=REF_DIR + "/" + REF_NAME + ".repma.bed",
    output:
        tmp=temp("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dp"),
        dp="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dpstats.txt",
    group:
        "realigned_bam_group"
    params:
        minDP=config["minDP"],
        maxDP=config["maxDP"],
        cov=config["zerocoverage"],
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_realigned_bam_depth.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        if [ {params.cov} = "True" ] # include sites with missing data / zero coverage
        then
            samtools depth -a -Q 30 -q 30 -b {input.no_rep_bed} {input.bam} > {output.tmp} 2> {log} &&
            awk '{{sum+=$3}} END {{ print sum/NR }}' {output.tmp} | awk -v min={params.minDP} -v max={params.maxDP} \
            '{{ printf "%.0f %.0f %.0f", $1, $1*min, $1*max }}' > {output.dp} 2>> {log}
        elif [ {params.cov} = "False" ] # exclude sites with missing data / zero coverage
        then
            samtools depth -Q 30 -q 30 -b {input.no_rep_bed} {input.bam} > {output.tmp} 2> {log} &&
            awk '{{sum+=$3}} END {{ print sum/NR }}' {output.tmp} | awk -v min={params.minDP} -v max={params.maxDP} \
            '{{ printf "%.0f %.0f %.0f", $1, $1*min, $1*max }}' > {output.dp} 2>> {log}
        fi
        """


rule plot_dp_hist:
    """Plot a depth histogram per sample with depth filtering thresholds and genome-wide average value"""
    input:
        dp=rules.realigned_bam_depth.output.dp,
        tmp=rules.realigned_bam_depth.output.tmp,
    output:
        pdf=report("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dp.hist.pdf",
            caption="../report/depth_plot.rst",
            category="BAM file processing",),
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/{dataset}/" + REF_NAME + "/{sample}_plot_dp_hist.log",
    script:
        "../scripts/dp_hist_plot.py"


rule historical_realigned_bam_multiqc:
    """Summarize all stats and qualimap results from all historical bam files until indel realignment"""
    input:
        expand("results/historical/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30{extension}",
            sample=hist_pipeline_bam_sm,
            extension=[".bam.dp.hist.pdf", 
                ".bam.stats.txt", 
                ".bam.qualimap/qualimapReport.html",
                "_fastqc.html",
                "_fastqc.zip",],),
    output:
        stats=report(
            "results/historical/mapping/" + REF_NAME + "/stats/bams_indels_realigned/multiqc/multiqc_report.html",
            caption="../report/historical_realigned_bam_multiqc.rst",
            category="BAM file processing",),
    params:
        indir="results/historical/mapping/" + REF_NAME + "/stats/bams_indels_realigned/",
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_indels_realigned/multiqc",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/historical/" + REF_NAME + "/historical_realigned_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_realigned_bam_multiqc:
    """Summarize all stats and qualimap results from all modern bam files"""
    input:
        expand("results/modern/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dp.hist.pdf",
            sample=mod_pipeline_bam_sm,
            extension=[".bam.dp.hist.pdf", 
                ".bam.stats.txt", 
                ".bam.qualimap/qualimapReport.html",
                "_fastqc.html",
                "_fastqc.zip",],),),
    output:
        stats=report(
            "results/modern/mapping/" + REF_NAME + "/stats/bams_indels_realigned/multiqc/multiqc_report.html",
            caption="../report/modern_realigned_bam_multiqc.rst",
            category="BAM file processing",),
    params:
        indir="results/modern/mapping/" + REF_NAME + "/stats/bams_indels_realigned/",
        outdir="results/modern/mapping/" + REF_NAME + "/stats/bams_indels_realigned/multiqc",
    log:
        "results/logs/3.1_bam_rmdup_realign_indels/modern/" + REF_NAME + "/modern_realigned_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
