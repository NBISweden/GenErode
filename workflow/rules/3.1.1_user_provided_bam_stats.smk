##########################################################################
### 3.1.1 Statistics for user-provided BAM files

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    if len(hist_user_bam_sm) > 0:
        all_outputs.append("results/historical/mapping/" + REF_NAME + "/stats/bams_user_provided/multiqc/multiqc_report.html")

if os.path.exists(config["modern_samples"]):
    if len(mod_user_bam_sm) > 0:
        all_outputs.append("results/modern/mapping/" + REF_NAME + "/stats/bams_user_provided/multiqc/multiqc_report.html")


# snakemake rules
localrules: userprovided_bam_historical_symbolic_links, userprovided_bam_modern_symbolic_links

rule userprovided_bam_historical_symbolic_links:
    """Make symbolic links to the user-provided bam files based on a metadata table"""
    input:
        ancient(config["historical_samples"]),
    output:
        bam="results/historical/mapping/" + REF_NAME + "/{sample}.userprovided.bam",
    params:
        abs_bam=lambda wildcards, output: os.path.abspath(output.bam),
    log:
        "data/logs/3.1.1_user_provided_bams/historical/" + REF_NAME + "/{sample}_userprovided_bam_historical_symbolic_links.log",
    run:
        SAMPLENAME = "{}".format(wildcards.sample)
        shell("""
        ln -s {bam} {{params.abs_bam}} 2> {{log}}
        """.format(SAMPLENAME=SAMPLENAME, **hist_user_bam_symlinks_dict[SAMPLENAME]))


rule userprovided_bam_modern_symbolic_links:
    """Make symbolic links to the user-provided bam files based on a metadata table"""
    input:
        ancient(config["modern_samples"]),
    output:
        bam="results/modern/mapping/" + REF_NAME + "/{sample}.userprovided.bam",
    params:
        abs_bam=lambda wildcards, output: os.path.abspath(output.bam),
    log:
        "data/logs/3.1.1_user_provided_bams/modern/" + REF_NAME + "/{sample}_userprovided_bam_modern_symbolic_links.log",
    run:
        SAMPLENAME = "{}".format(wildcards.sample)
        shell("""
        ln -s {bam} {{params.abs_bam}} 2> {{log}}
        """.format(SAMPLENAME=SAMPLENAME, **mod_user_bam_symlinks_dict[SAMPLENAME]))


rule index_userprovided_bams:
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.userprovided.bam",
    output:
        index=temp("results/{dataset}/mapping/" + REF_NAME + "/{sample}.userprovided.bam.bai"),
    log:
        "results/logs/3.1.1_user_provided_bams/{dataset}/" + REF_NAME + "/{sample}_index_userprovided_bams.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools index {input.bam} {output.index} 2> {log}
        """


rule userprovided_bam_stats:
    """Basic statistics on mapping output"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.userprovided.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}.userprovided.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.bam.stats.txt",
    group:
        "userprovided_bam_group"
    log:
        "results/logs/3.1.1_user_provided_bams/{dataset}/" + REF_NAME + "/{sample}_userprovided_bam_stats.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule userprovided_bam_fastqc:
    """Run fastqc on all userprovided bam files"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.userprovided.bam",
    output:
        html="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/fastqc/{sample}.userprovided_fastqc.html",
        zip="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/fastqc/{sample}.userprovided_fastqc.zip",
        dir=directory("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/fastqc/{sample}.userprovided_fastqc"),
    params:
        dir="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/fastqc",
    log:
        "results/logs/3.1.1_user_provided_bams/{dataset}/" + REF_NAME + "/{sample}_userprovided_bam_fastqc.log",
    threads: 2
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    shell:
        """
        fastqc -o {params.dir} -t {threads} --extract {input.bam} 2> {log}
        """


rule userprovided_bam_qualimap:
    """More detailed stats"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.userprovided.bam",
        index="results/{dataset}/mapping/" + REF_NAME + "/{sample}.userprovided.bam.bai",
    output:
        stats="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.bam.qualimap/qualimapReport.html",
        results="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.bam.qualimap/genome_results.txt",
        outdir=directory("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.bam.qualimap"),
    threads: 8
    resources:
        mem_mb=64000,
    params:
        outdir="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.bam.qualimap",
    log:
        "results/logs/3.1.1_user_provided_bams/{dataset}/" + REF_NAME + "/{sample}_userprovided_bam_qualimap.log",
    singularity:
        "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G -nt {threads} -outdir {params.outdir} -outformat html 2> {log}
        """


rule userprovided_bam_depth:
    """Get average genome-wide depth per bam file, excluding repeats and filtering for high quality"""
    """Will be used to filter out sites outside the estimated depth thresholds"""
    """samtools depth to get depth per site, piped into awk to get average across all sites, 
    piped into awk to calculate 1/3*average and 10*average depth (rounded to integer)"""
    input:
        bam="results/{dataset}/mapping/" + REF_NAME + "/{sample}.userprovided.bam",
        no_rep_bed=REF_DIR + "/" + REF_NAME + ".repma.bed",
    output:
        tmp=temp("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.repma.Q30.bam.dp"),
        dp="results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.repma.Q30.bam.dpstats.txt",
    params:
        minDP=config["minDP"],
        maxDP=config["maxDP"],
        cov=config["zerocoverage"],
    log:
        "results/logs/3.1.1_user_provided_bams/{dataset}/" + REF_NAME + "/{sample}_userprovided_bam_depth.log",
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


rule plot_userprovided_bam_dp_hist:
    """Plot a depth histogram per sample with depth filtering thresholds and genome-wide average value"""
    input:
        dp=rules.userprovided_bam_depth.output.dp,
        tmp=rules.userprovided_bam_depth.output.tmp,
    output:
        pdf=report("results/{dataset}/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.repma.Q30.bam.dp.hist.pdf",
            caption="../report/depth_plot.rst",
            category="BAM file processing",),
    log:
        "results/logs/3.1.1_user_provided_bams/{dataset}/" + REF_NAME + "/{sample}_plot_userprovided_bam_dp_hist.log",
    script:
        "../scripts/dp_hist_plot.py"


rule historical_userprovided_bam_multiqc:
    """Summarize all stats and qualimap results from all historical bam files until indel realignment"""
    input:
        dp=expand("results/historical/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.repma.Q30.{extension}",
            sample=hist_user_bam_sm,
            extension=[".bam.dp.hist.pdf", 
                ".bam.stats.txt", 
                ".bam.qualimap/qualimapReport.html",
                "_fastqc.html",
                "_fastqc.zip",],),
    output:
        stats=report(
            "results/historical/mapping/" + REF_NAME + "/stats/bams_user_provided/multiqc/multiqc_report.html",
            caption="../report/historical_userprovided_bam_multiqc.rst",
            category="BAM file processing",),
    params:
        indir="results/historical/mapping/" + REF_NAME + "/stats/bams_user_provided/",
        outdir="results/historical/mapping/" + REF_NAME + "/stats/bams_user_provided/multiqc",
    log:
        "results/logs/3.1.1_user_provided_bams/historical/" + REF_NAME + "/historical_userprovided_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_userprovided_bam_multiqc:
    """Summarize all stats and qualimap results from all modern bam files"""
    input:
        dp=expand("results/modern/mapping/" + REF_NAME + "/stats/bams_user_provided/{sample}.userprovided.repma.Q30.{extension}",
            sample=mod_user_bam_sm,
            extension=[".bam.dp.hist.pdf", 
                ".bam.stats.txt", 
                ".bam.qualimap/qualimapReport.html",
                "_fastqc.html",
                "_fastqc.zip",],),
    output:
        stats=report(
            "results/modern/mapping/" + REF_NAME + "/stats/bams_user_provided/multiqc/multiqc_report.html",
            caption="../report/modern_userprovided_bam_multiqc.rst",
            category="BAM file processing",),
    params:
        indir="results/modern/mapping/" + REF_NAME + "/stats/bams_user_provided/",
        outdir="results/modern/mapping/" + REF_NAME + "/stats/bams_user_provided/multiqc",
    log:
        "results/logs/3.1.1_user_provided_bams/modern/" + REF_NAME + "/modern_userprovided_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """