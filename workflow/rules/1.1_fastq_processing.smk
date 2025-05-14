##########################################################################
### 1. Adapter removal and read merging

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    all_outputs.append("data/raw_reads_symlinks/historical/stats/multiqc/multiqc_report.html")
    all_outputs.append("results/historical/trimming/stats/multiqc/multiqc_report.html")

if os.path.exists(config["modern_samples"]):
    all_outputs.append("data/raw_reads_symlinks/modern/stats/multiqc/multiqc_report.html")
    all_outputs.append("results/modern/trimming/stats/multiqc/multiqc_report.html")


# snakemake rules
localrules: fastq_historical_symbolic_links, fastq_modern_symbolic_links

rule fastq_historical_symbolic_links:
    """Make symbolic links to the original fastq-files based on a metadata table"""
    input:
        ancient(config["historical_samples"]),
    output:
        fastq_r1="data/raw_reads_symlinks/historical/{sample}_{index}_{lane}_R1.fastq.gz",
        fastq_r2="data/raw_reads_symlinks/historical/{sample}_{index}_{lane}_R2.fastq.gz",
    params:
        abs_fastq_r1=lambda wildcards, output: os.path.abspath(output.fastq_r1),
        abs_fastq_r2=lambda wildcards, output: os.path.abspath(output.fastq_r2),
    log:
        "data/logs/1.1_fastq_processing/historical/{sample}_{index}_{lane}_fastq_historical_symbolic_links.log",
    run:
        FULLSAMPLENAME = "{}_{}_{}".format(wildcards.sample, wildcards.index, wildcards.lane)
        shell("""
        ln -s {R1} {{params.abs_fastq_r1}} 2> {{log}}
        ln -s {R2} {{params.abs_fastq_r2}} 2>> {{log}}
        """.format(FULLSAMPLENAME=FULLSAMPLENAME, **hist_fastq_symlinks_dict[FULLSAMPLENAME]))


rule fastqc_historical_raw:
    """Run fastqc on all samples"""
    input:
        fastq="data/raw_reads_symlinks/historical/{sample}_{index}_{lane}_R{nr}.fastq.gz",
    output:
        html="data/raw_reads_symlinks/historical/stats/{sample}_{index}_{lane}_R{nr}_fastqc.html",
        zip="data/raw_reads_symlinks/historical/stats/{sample}_{index}_{lane}_R{nr}_fastqc.zip",
        dir=directory("data/raw_reads_symlinks/historical/stats/{sample}_{index}_{lane}_R{nr}_fastqc"),
    params:
        dir="data/raw_reads_symlinks/historical/stats",
    log:
        "data/logs/1.1_fastq_processing/historical/{sample}_{index}_{lane}_R{nr}_fastqc_historical_raw.log",
    threads: 2
    singularity:
        fastqc_container
    shell:
        """
        fastqc -o {params.dir} -t {threads} --extract {input.fastq} 2> {log}
        """


rule multiqc_historical_raw:
    """Summarize all fastqc results for raw historical samples"""
    input:
        expand("data/raw_reads_symlinks/historical/stats/{sampleindexlane}_R{nr}_fastqc.html",
            sampleindexlane=hist_pipeline_bam_sm_idx_ln,
            nr=[1, 2],),
        expand("data/raw_reads_symlinks/historical/stats/{sampleindexlane}_R{nr}_fastqc.zip",
            sampleindexlane=hist_pipeline_bam_sm_idx_ln,
            nr=[1, 2],),
    output:
        stats="data/raw_reads_symlinks/historical/stats/multiqc/multiqc_report.html",
    params:
        indir="data/raw_reads_symlinks/historical/stats/",
        outdir="data/raw_reads_symlinks/historical/stats/multiqc/",
    log:
        "data/logs/1.1_fastq_processing/historical/multiqc_historical_raw.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule fastq_modern_symbolic_links:
    """Make symbolic links to the original fastq-files based on a metadata table"""
    input:
        ancient(config["modern_samples"]),
    output:
        fastq_r1="data/raw_reads_symlinks/modern/{sample}_{index}_{lane}_R1.fastq.gz",
        fastq_r2="data/raw_reads_symlinks/modern/{sample}_{index}_{lane}_R2.fastq.gz",
    params:
        abs_fastq_r1=lambda wildcards, output: os.path.abspath(output.fastq_r1),
        abs_fastq_r2=lambda wildcards, output: os.path.abspath(output.fastq_r2),
    log:
        "data/logs/1.1_fastq_processing/modern/{sample}_{index}_{lane}_fastq_modern_symbolic_links.log",
    run:
        FULLSAMPLENAME = "{}_{}_{}".format(wildcards.sample, wildcards.index, wildcards.lane)
        shell("""
        ln -s {R1} {{params.abs_fastq_r1}} 2> {{log}}
        ln -s {R2} {{params.abs_fastq_r2}} 2>> {{log}}
            """.format(FULLSAMPLENAME=FULLSAMPLENAME, **mod_fastq_symlinks_dict[FULLSAMPLENAME]))


rule fastqc_modern_raw:
    """Run fastqc on all samples"""
    input:
        fastq="data/raw_reads_symlinks/modern/{sample}_{index}_{lane}_R{nr}.fastq.gz",
    output:
        html="data/raw_reads_symlinks/modern/stats/{sample}_{index}_{lane}_R{nr}_fastqc.html",
        zip="data/raw_reads_symlinks/modern/stats/{sample}_{index}_{lane}_R{nr}_fastqc.zip",
        dir=directory("data/raw_reads_symlinks/modern/stats/{sample}_{index}_{lane}_R{nr}_fastqc"),
    params:
        dir="data/raw_reads_symlinks/modern/stats",
    log:
        "data/logs/1.1_fastq_processing/modern/{sample}_{index}_{lane}_R{nr}_fastqc_modern_raw.log",
    threads: 2
    singularity:
        fastqc_container
    shell:
        """
        fastqc -o {params.dir} -t {threads} --extract {input.fastq} 2> {log}
        """


rule multiqc_modern_raw:
    """Summarize all fastqc results for all raw modern samples"""
    input:
        expand("data/raw_reads_symlinks/modern/stats/{sampleindexlane}_R{nr}_fastqc.html",
            sampleindexlane=mod_pipeline_bam_sm_idx_ln,
            nr=[1, 2],),
        expand("data/raw_reads_symlinks/modern/stats/{sampleindexlane}_R{nr}_fastqc.zip",
            sampleindexlane=mod_pipeline_bam_sm_idx_ln,
            nr=[1, 2],),
    output:
        stats="data/raw_reads_symlinks/modern/stats/multiqc/multiqc_report.html",
    params:
        indir="data/raw_reads_symlinks/modern/stats/",
        outdir="data/raw_reads_symlinks/modern/stats/multiqc/",
    log:
        "data/logs/1.1_fastq_processing/modern/multiqc_modern_raw.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule fastp_historical:
    """Remove adapters, quality trim (phred 15) and merge overlapping paired-end reads in historical samples.
    fastp automatically detects adapter sequences for removal.
    NextSeq and NovaSeq samples are automatically detected and poly-G tails are removed.
    Minimum read length specified in config file.
    """
    input:
        R1=rules.fastq_historical_symbolic_links.output.fastq_r1,
        R2=rules.fastq_historical_symbolic_links.output.fastq_r2,
    output:
        R1_un=temp("results/historical/trimming/{sample}_{index}_{lane}_R1_unmerged.fastq.gz"),
        R2_un=temp("results/historical/trimming/{sample}_{index}_{lane}_R2_unmerged.fastq.gz"),
        merged="results/historical/trimming/{sample}_{index}_{lane}_trimmed_merged.fastq.gz",
        html="results/historical/trimming/stats/{sample}_{index}_{lane}_fastp_report.html",
        json=temp("results/historical/trimming/stats/{sample}_{index}_{lane}_fastp_report.json"),
    params:
        readlength=config["hist_readlength"],
        report="fastp report for {sample}_{index}_{lane}",
    log:
        "results/logs/1.1_fastq_processing/historical/{sample}_{index}_{lane}_fastp_historical.log",
    threads: 4
    singularity:
        fastp_container
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -p -c --merge --overlap_len_require 15 --overlap_diff_limit 1 \
        --merged_out={output.merged} -o {output.R1_un} -O {output.R2_un} -h {output.html} -j {output.json} \
        -R '{params.report}' -w {threads} -l {params.readlength} 2> {log}
        """


rule fastp_modern:
    """Remove adapters from modern samples and quality trim (phred 15).
    fastp automatically detects adapter sequences for removal.
    NextSeq and NovaSeq samples are automatically detected and poly-G tails are removed.
    Minimum read length specified in config file.
    """
    input:
        R1=rules.fastq_modern_symbolic_links.output.fastq_r1,
        R2=rules.fastq_modern_symbolic_links.output.fastq_r2,
    output:
        R1_trimmed="results/modern/trimming/{sample}_{index}_{lane}_R1_trimmed.fastq.gz",
        R2_trimmed="results/modern/trimming/{sample}_{index}_{lane}_R2_trimmed.fastq.gz",
        html="results/modern/trimming/stats/{sample}_{index}_{lane}_fastp_report.html",
        json=temp("results/modern/trimming/stats/{sample}_{index}_{lane}_fastp_report.json"),
    params:
        readlength=config["mod_readlength"],
        report="fastp report for {sample}_{index}_{lane}",
    log:
        "results/logs/1.1_fastq_processing/modern/{sample}_{index}_{lane}_fastp_modern.log",
    threads: 4
    singularity:
        fastp_container
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -p -c -o {output.R1_trimmed} -O {output.R2_trimmed} \
        -h {output.html} -j {output.json} -R '{params.report}' -w {threads} -l {params.readlength} 2> {log}
        """


rule fastqc_historical_merged:
    """Run fastqc on all samples after adapter trimming and merging"""
    input:
        rules.fastp_historical.output.merged,
    output:
        html="results/historical/trimming/stats/{sample}_{index}_{lane}_trimmed_merged_fastqc.html",
        zip="results/historical/trimming/stats/{sample}_{index}_{lane}_trimmed_merged_fastqc.zip",
        dir=directory("results/historical/trimming/stats/{sample}_{index}_{lane}_trimmed_merged_fastqc"),
    params:
        dir="results/historical/trimming/stats",
    log:
        "results/logs/1.1_fastq_processing/historical/{sample}_{index}_{lane}_fastqc_historical_merged.log",
    threads: 2
    singularity:
        fastqc_container
    shell:
        """
        fastqc -o {params.dir} -t {threads} --extract {input} 2> {log}
        """


rule fastqc_historical_unmerged:
    """Run fastqc on all unmerged reads after adapter trimming and merging"""
    input:
        "results/historical/trimming/{sample}_{index}_{lane}_R{nr}_unmerged.fastq.gz",
    output:
        html="results/historical/trimming/stats/{sample}_{index}_{lane}_R{nr}_unmerged_fastqc.html",
        zip="results/historical/trimming/stats/{sample}_{index}_{lane}_R{nr}_unmerged_fastqc.zip",
        dir=directory("results/historical/trimming/stats/{sample}_{index}_{lane}_R{nr}_unmerged_fastqc"),
    params:
        dir="results/historical/trimming/stats",
    log:
        "results/logs/1.1_fastq_processing/historical/{sample}_{index}_{lane}_R{nr}_unmerged_fastqc_historical_unmerged.log",
    threads: 2
    singularity:
        fastqc_container
    shell:
        """
        if [ -s {input} ]
        then
            fastqc -o {params.dir} -t {threads} --extract {input} 2> {log}
        else
            mkdir -p {output.dir} &&
            touch {output.html} && touch {output.zip} &&
            echo "No reads in {input} >> {log}"
        fi
        """


rule fastqc_modern_trimmed:
    """Run fastqc on all samples after adapter trimming and merging"""
    input:
        "results/modern/trimming/{sample}_{index}_{lane}_R{nr}_trimmed.fastq.gz",
    output:
        html="results/modern/trimming/stats/{sample}_{index}_{lane}_R{nr}_trimmed_fastqc.html",
        zip="results/modern/trimming/stats/{sample}_{index}_{lane}_R{nr}_trimmed_fastqc.zip",
        dir=directory("results/modern/trimming/stats/{sample}_{index}_{lane}_R{nr}_trimmed_fastqc"),
    params:
        dir="results/modern/trimming/stats/",
    log:
        "results/logs/1.1_fastq_processing/modern/{sample}_{index}_{lane}_R{nr}_trimmed_fastqc_modern_trimmed.log",
    threads: 2
    singularity:
        fastqc_container
    shell:
        """
        fastqc -o {params.dir} -t {threads} --extract {input} 2> {log}
        """


rule multiqc_historical_trimmed:
    """Summarize all fastqc results for all trimmed and merged historical samples"""
    input:
        expand("results/historical/trimming/stats/{sampleindexlane}_trimmed_merged_fastqc.{ext}",
            sampleindexlane=hist_pipeline_bam_sm_idx_ln,
            ext=["html", "zip"],),
        expand("results/historical/trimming/stats/{sampleindexlane}_fastp_report.html",
            sampleindexlane=hist_pipeline_bam_sm_idx_ln,),
        expand("results/historical/trimming/stats/{sampleindexlane}_R{nr}_unmerged_fastqc.{ext}",
            sampleindexlane=hist_pipeline_bam_sm_idx_ln,
            nr=[1, 2],
            ext=["html", "zip"],),
    output:
        stats="results/historical/trimming/stats/multiqc/multiqc_report.html",
    params:
        indir="results/historical/trimming/stats/",
        outdir="results/historical/trimming/stats/multiqc/",
    log:
        "results/logs/1.1_fastq_processing/historical/multiqc_historical_trimmed.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule multiqc_modern_trimmed:
    """Summarize all fastqc results for all trimmed modern samples"""
    input:
        expand("results/modern/trimming/stats/{sampleindexlane}_R{nr}_trimmed_fastqc.{ext}",
            sampleindexlane=mod_pipeline_bam_sm_idx_ln,
            nr=[1, 2],
            ext=["html", "zip"]),
        expand("results/modern/trimming/stats/{sampleindexlane}_fastp_report.html",
            sampleindexlane=mod_pipeline_bam_sm_idx_ln,)
    output:
        stats="results/modern/trimming/stats/multiqc/multiqc_report.html",
    params:
        indir="results/modern/trimming/stats/",
        outdir="results/modern/trimming/stats/multiqc/",
    log:
        "results/logs/1.1_fastq_processing/modern/multiqc_modern_trimmed.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
