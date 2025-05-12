##########################################################################
### 3.2 OPTIONAL Base quality rescaling in historical BAM files with ancient DNA damage patterns

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    if len(HIST_RESCALED_SAMPLES) > 0:
        all_outputs.append("results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/multiqc/multiqc_report.html")


# Functions for this step of the pipeline
# return correct bam file path depending on bam file type
def rescale_historical_bam_inputs(wildcards):
    if wildcards.sample in hist_pipeline_bam_sm:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam"
    elif wildcards.sample in hist_user_bam_sm:
        return "results/historical/mapping/" + REF_NAME + "/{sample}.userprovided.bam"


def rescaled_bam_multiqc_inputs(wildcards):
    """Collect all inputs for multiqc"""
    # Get the list of all rescaled bam files
    pipeline_bams_rescaled = expand("results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.merged.rmdup.merged.realn.rescaled{extension}",
        sample=HIST_PIPELINE_RESCALED_SAMPLES,
        extension=[".bam.stats.txt", 
        ".bam.qualimap/qualimapReport.html", 
        ".bam.qualimap/genome_results.txt",],)
    pipeline_bams_rescaled_fastqc = expand("results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/fastqc/{sample}.merged.rmdup.merged.realn.rescaled{extension}",
        sample=HIST_PIPELINE_RESCALED_SAMPLES,
        extension=["_fastqc.html",
        "_fastqc.zip",],)
    userprocessed_bams_rescaled = expand("results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.userprovided.rescaled{extension}",
        sample=HIST_USER_RESCALED_SAMPLES,
        extension=[".bam.stats.txt", 
        ".bam.qualimap/qualimapReport.html", 
        ".bam.qualimap/genome_results.txt",],)
    userprocessed_bams_rescaled_fastqc = expand("results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/fastqc/{sample}.userprovided.rescaled{extension}",
        sample=HIST_USER_RESCALED_SAMPLES,
        extension=["_fastqc.html",
        "_fastqc.zip",],)
    return pipeline_bams_rescaled + pipeline_bams_rescaled_fastqc + userprocessed_bams_rescaled + userprocessed_bams_rescaled_fastqc


# snakemake rules
rule rescale_historical:
    """Rescale base quality scores of likely damaged positions"""
    input:
        bam=rescale_historical_bam_inputs,
        ref=config["ref_path"],
    output:
        fragmis="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/Fragmisincorporation_plot.pdf",
        length="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/Length_plot.pdf",
        dna="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/dnacomp.txt",
        dna_csv="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/dnacomp_genome.csv",
        misin="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/misincorporation.txt",
        mcmcstat="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/Stats_out_MCMC_iter_summ_stat.csv",
        mcmcprob="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/Stats_out_MCMC_correct_prob.csv",
        mcmchist="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/Stats_out_MCMC_hist.pdf",
        mcmciter="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/Stats_out_MCMC_iter.csv",
        mcmcpost="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/Stats_out_MCMC_post_pred.pdf",
        mcmctrace="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/Stats_out_MCMC_trace.pdf",
        GtoA_freq="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/3pGtoA_freq.txt",
        CtoT_freq="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/5pCtoT_freq.txt",
        lgdist="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/lgdistribution.txt",
        log="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/Runtime_log.txt",
        tmp=temp("results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/{sample}.{processed}.rescaled.bam"),
        rescaled="results/historical/mapping/" + REF_NAME + "/{sample}.{processed}.rescaled.bam",
    threads: 4
    params:
        dir="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.bam.mapDamage/",
    log:
        "results/logs/3.2_historical_bam_mapDamage/" + REF_NAME + "/{sample}.{processed}_rescale_historical.log",
    singularity:
        "docker://biocontainers/mapdamage:v2.0.9dfsg-1-deb_cv1"
    shell:
        """
        mapDamage -i {input.bam} -r {input.ref} -d {params.dir} --merge-reference-sequences --rescale 2> {log} &&
        cp {output.tmp} {output.rescaled} 2>> {log}
        """


rule index_rescaled_bams:
    input:
        bam="results/historical/mapping/" + REF_NAME + "/{sample}.{processed}.rescaled.bam",
    output:
        index="results/historical/mapping/" + REF_NAME + "/{sample}.{processed}.rescaled.bam.bai",
    log:
        "results/logs/3.2_historical_bam_mapDamage/" + REF_NAME + "/historical/{sample}.{processed}_index_rescaled_bams.log",
    group:
        "rescaled_bam_group"
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools index {input.bam} {output.index} 2> {log}
        """


rule rescaled_bam_stats:
    """Basic statistics on mapping output"""
    input:
        bam=rules.rescale_historical.output.rescaled,
        bai="results/historical/mapping/" + REF_NAME + "/{sample}.{processed}.rescaled.bam.bai",
    output:
        stats="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.rescaled.bam.stats.txt",
    group:
        "rescaled_bam_group"
    params:
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/",
    log:
        "results/logs/3.2_historical_bam_mapDamage/" + REF_NAME + "/{sample}.{processed}_rescaled_bam_stats.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule rescaled_bam_fastqc:
    """Run fastqc on rescaled bam files"""
    input:
        bam=rules.rescale_historical.output.rescaled,
        bai="results/historical/mapping/" + REF_NAME + "/{sample}.{processed}.rescaled.bam.bai",
    output:
        html="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/fastqc/{sample}.{processed}.rescaled_fastqc.html",
        zip="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/fastqc/{sample}.{processed}.rescaled_fastqc.zip",
        dir=directory("results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/fastqc/{sample}.{processed}.rescaled_fastqc"),
    params:
        dir="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/fastqc",
    log:
        "results/logs/3.2_historical_bam_mapDamage/" + REF_NAME + "/{sample}.{processed}_rescaled_bam_fastqc.log",
    threads: 2
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    shell:
        """
        fastqc -o {params.dir} -t {threads} --extract {input.bam} 2> {log}
        """


rule rescaled_bam_qualimap:
    """More detailed stats"""
    input:
        bam=rules.rescale_historical.output.rescaled,
        bai="results/historical/mapping/" + REF_NAME + "/{sample}.{processed}.rescaled.bam.bai",
    output:
        stats="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.rescaled.bam.qualimap/qualimapReport.html",
        results="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.rescaled.bam.qualimap/genome_results.txt",
        outdir=directory("results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.rescaled.bam.qualimap"),
    threads: 8
    resources:
        mem_mb=64000,
    params:
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/{sample}.{processed}.rescaled.bam.qualimap",
    log:
        "results/logs/3.2_historical_bam_mapDamage/" + REF_NAME + "/{sample}.{processed}_rescaled_bam_qualimap.log",
    singularity:
        "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G  -nt {threads} -outdir {params.outdir} -outformat html 2> {log}
        """


rule rescaled_bam_multiqc:
    """Summarize all stats results from all historical bam files after rescaling"""
    input:
        rescaled_bam_multiqc_inputs
    output:
        stats="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/multiqc/multiqc_report.html",
    params:
        indir="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/",
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_rescaled/multiqc",
    log:
        "results/logs/3.2_historical_bam_mapDamage/" + REF_NAME + "/rescaled_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
