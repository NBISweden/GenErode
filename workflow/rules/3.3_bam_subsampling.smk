##########################################################################
### 3.3 OPTIONAL BAM file subsampling to same depth

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    if len(HIST_SUBSAMPLED_SAMPLES) > 0:
        all_outputs.append("results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/multiqc/multiqc_report.html")

if os.path.exists(config["modern_samples"]):
    if len(MODERN_SUBSAMPLED_SAMPLES) > 0:
        all_outputs.append("results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/multiqc/multiqc_report.html")


# Functions for this step of the pipeline
def historical_subsampled_bam_multiqc_inputs(wildcards):
    """Collect all inputs for multiqc"""
    not_rescaled_subsampled_pipeline_bams=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}{extension}",
        sample=HIST_PIPELINE_NOT_RESCALED_SUBSAMPLED_SAMPLES,
        DP=config["subsampling_depth"],
        extension=[".bam.stats.txt", 
            ".bam.qualimap/qualimapReport.html", 
            ".bam.qualimap/genome_results.txt",
            ".repma.Q30.bam.dpstats.txt"],)
    rescaled_subsampled_pipeline_bams=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}{extension}",
        sample=HIST_PIPELINE_RESCALED_SUBSAMPLED_SAMPLES,
        DP=config["subsampling_depth"],
        extension=[".bam.stats.txt", 
            ".bam.qualimap/qualimapReport.html", 
            ".bam.qualimap/genome_results.txt",
            ".repma.Q30.bam.dpstats.txt"],)
    not_rescaled_subsampled_user_bams=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.userprovided.mapped_q30.subs_dp{DP}{extension}",
        sample=HIST_USER_NOT_RESCALED_SUBSAMPLED_SAMPLES,
        DP=config["subsampling_depth"],
        extension=[".bam.stats.txt", 
            ".bam.qualimap/qualimapReport.html", 
            ".bam.qualimap/genome_results.txt",
            ".repma.Q30.bam.dpstats.txt"],)
    rescaled_subsampled_user_bams=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.userprovided.rescaled.mapped_q30.subs_dp{DP}{extension}",
        sample=HIST_USER_RESCALED_SUBSAMPLED_SAMPLES,
        DP=config["subsampling_depth"],
        extension=[".bam.stats.txt", 
            ".bam.qualimap/qualimapReport.html", 
            ".bam.qualimap/genome_results.txt",
            ".repma.Q30.bam.dpstats.txt"],)
    return not_rescaled_subsampled_pipeline_bams + rescaled_subsampled_pipeline_bams + not_rescaled_subsampled_user_bams + rescaled_subsampled_user_bams

def modern_subsampled_bam_multiqc_inputs(wildcards):
    """Collect all inputs for multiqc"""
    subsampled_pipeline_bams=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}{extension}",
        sample=MODERN_PIPELINE_SUBSAMPLED_SAMPLES,
        DP=config["subsampling_depth"],
        extension=[".bam.stats.txt", 
            ".bam.qualimap/qualimapReport.html", 
            ".bam.qualimap/genome_results.txt",
            ".repma.Q30.bam.dpstats.txt"],)
    subsampled_user_bams=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.userprovided.mapped_q30.subs_dp{DP}{extension}",
        sample=MODERN_USER_SUBSAMPLED_SAMPLES,
        DP=config["subsampling_depth"],
        extension=[".bam.stats.txt", 
            ".bam.qualimap/qualimapReport.html", 
            ".bam.qualimap/genome_results.txt",
            ".repma.Q30.bam.dpstats.txt"],)
    return subsampled_pipeline_bams + subsampled_user_bams


# snakemake rules
rule filter_bam_mapped_mq:
    """Remove unmapped reads before subsampling"""
    """Keep only reads with mapping quality > 30"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.{processed}.bam",
    output:
        filtered=temp("results/{dataset}/mapping/" + REF_NAME + "/{sample}.{processed}.mapped_q30.bam"),
    threads: 2
    log:
        "results/logs/3.3_bam_subsampling/{dataset}/" + REF_NAME + "/{sample}.{processed}_filter_bam_mapped_mq.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools view -h -b -F 4 -q 30 -@ {threads} -o {output.filtered} {input.bam} 2> {log}
        """


rule subsample_bams:
    """Subsample bam files to target depth per sample"""
    input:
        bam=rules.filter_bam_mapped_mq.output.filtered,
        dp="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.{processed}.repma.Q30.bam.dpstats.txt",
    output:
        subsam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam",
    threads: 2
    params:
        DP=config["subsampling_depth"],
    log:
        "results/logs/3.3_bam_subsampling/{dataset}/" + REF_NAME + "/{sample}.{processed}.subs_dp{DP}_subsample_bams.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        depth=`head -n 1 {input.dp} | cut -d' ' -f 1`
        frac=`awk -v s={params.DP} -v d=$depth "BEGIN {{print s/d}}"`
        if [ `awk 'BEGIN {{print ('$frac' <= 1.0)}}'` = 1 ] # awk will return 1 if the statement is true, and 0 if it is false
        then
            samtools view -h -b -s $frac -@ {threads} -o {output.subsam} {input.bam} 2> {log}
        else
            echo "!!!\nWarning [genome erosion workflow]: The sample {input.bam} has a lower average depth than the target depth for subsampling. \
            Remove the sample from the subsampling list in the config file or choose a lower target depth.\n!!!" >> {log}
        fi
        """


rule index_subsampled_bams:
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam",
    output:
        index=temp("results/{dataset}/mapping/" + REF_NAME + "/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam.bai"),
    log:
        "results/logs/3.3_bam_subsampling/{dataset}/" + REF_NAME + "/{sample}.{processed}.subs_dp{DP}_index_subsampled_bams.log",
    group:
        "subsampled_bam_group"
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools index {input.bam} {output.index} 2> {log}
        """


rule subsampled_bam_stats:
    """Basic statistics on mapping output"""
    input:
        bam=rules.subsample_bams.output.subsam,
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam.stats.txt",
    group:
        "subsampled_bam_group"
    log:
        "results/logs/3.3_bam_subsampling/{dataset}/" + REF_NAME + "/{sample}.{processed}.subs_dp{DP}_subsampled_bam_stats.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule subsampled_bam_depth:
    """Get average genome-wide depth per bam file, excluding repeats and filtering for high quality"""
    """Will be used to filter out sites outside the estimated depth thresholds"""
    """samtools depth to get depth per site, piped into awk to get average across all sites, 
    piped into awk to calculate 1/3*average and 2*average depth (rounded to integer)"""
    input:
        bam=rules.subsample_bams.output.subsam,
        no_rep_bed=REF_DIR + "/" + REF_NAME + ".repma.bed",
    output:
        tmp=temp("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.{processed}.mapped_q30.subs_dp{DP}.repma.Q30.bam.dp"),
        dp="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.{processed}.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt",
    group:
        "subsampled_bam_group"
    params:
        minDP=config["minDP"],
        maxDP=config["maxDP"],
        cov=config["zerocoverage"],
    log:
        "results/logs/3.3_bam_subsampling/{dataset}/" + REF_NAME + "/{sample}.{processed}.subs_dp{DP}_subsampled_bam_depth.log",
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


rule subsampled_bam_qualimap:
    """More detailed stats"""
    input:
        bam=rules.subsample_bams.output.subsam,
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam.qualimap/qualimapReport.html",
        results="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam.qualimap/genome_results.txt",
        outdir=directory("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam.qualimap"),
    threads: 8
    resources:
        mem_mb=64000,
    params:
        outdir="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.{processed}.mapped_q30.subs_dp{DP}.bam.qualimap",
    log:
        "results/logs/3.3_bam_subsampling/{dataset}/" + REF_NAME + "/{sample}.{processed}.subs_dp{DP}_subsampled_bam_qualimap.log",
    singularity:
        "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G -nt {threads} -outdir {params.outdir} -outformat html 2> {log}
        """


rule historical_subsampled_bam_multiqc:
    """Summarize all stats and qualimap results from subsampled historical bam files"""
    input:
        historical_subsampled_bam_multiqc_inputs,
    output:
        stats="results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/multiqc/multiqc_report.html",
    params:
        indir="results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/",
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/multiqc",
    log:
        "results/logs/3.3_bam_subsampling/historical/historical_subsampled_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_subsampled_bam_multiqc:
    """Summarize all stats and qualimap results from subsampled modern bam files"""
    input:
        modern_subsampled_bam_multiqc_inputs,
    output:
        stats="results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/multiqc/multiqc_report.html",
    params:
        indir="results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/",
        outdir="results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/multiqc",
    log:
        "results/logs/3.3_bam_subsampling/modern/modern_subsampled_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
