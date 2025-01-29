##########################################################################
### 13. Compute GERP++ scores
# Authors: Marcin Kierczak, Tom van der Valk, Verena Kutschera

# Code collecting output files from this part of the pipeline
all_outputs.append(expand("results/gerp/" + REF_NAME + ".{chr}.ancestral.rates.gerp.hist.pdf",
    chr=CHR,))

if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
    all_outputs.append(expand("results/gerp/all/" + REF_NAME + ".all.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_plot.pdf",
            fmiss=config["f_missing"],
            chr=CHR,
            minGERP=config["min_gerp"],
            maxGERP=config["max_gerp"],))

elif os.path.exists(config["historical_samples"]):
    all_outputs.append(expand("results/gerp/historical/" + REF_NAME + ".historical.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_plot.pdf",
            fmiss=config["f_missing"],
            chr=CHR,
            minGERP=config["min_gerp"],
            maxGERP=config["max_gerp"],))

elif os.path.exists(config["modern_samples"]):
    all_outputs.append(expand("results/gerp/modern/" + REF_NAME + ".modern.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_plot.pdf",
            fmiss=config["f_missing"],
            chr=CHR,
            minGERP=config["min_gerp"],
            maxGERP=config["max_gerp"],))


# Functions used by rules of this part of the pipeline
def rel_load_table_inputs(wildcards):
    """Collect output files for pipeline report"""
    outlist = []
    if os.path.exists(config["historical_samples"]):
        rescaled_not_subsampled_not_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
            sample=HIST_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,
            minGERP=config["min_gerp"],
            maxGERP=config["max_gerp"],)
        not_rescaled_not_subsampled_not_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
            sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,
            minGERP=config["min_gerp"],
            maxGERP=config["max_gerp"],)
        rescaled_subsampled_not_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
            sample=HIST_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
            DP=config["subsampling_depth"],
            fmiss=config["f_missing"],
            chr=CHR,
            minGERP=config["min_gerp"],
            maxGERP=config["max_gerp"],)
        not_rescaled_subsampled_not_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
            sample=HIST_NOT_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
            DP=config["subsampling_depth"],
            fmiss=config["f_missing"],
            chr=CHR,
            minGERP=config["min_gerp"],
            maxGERP=config["max_gerp"],)
        outlist += (rescaled_not_subsampled_not_CpG + not_rescaled_not_subsampled_not_CpG + rescaled_subsampled_not_CpG + not_rescaled_subsampled_not_CpG)
        if config["CpG_from_vcf"] == True:
            rescaled_not_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            not_rescaled_not_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            rescaled_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                DP=config["subsampling_depth"],
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            not_rescaled_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                DP=config["subsampling_depth"],
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG)
        elif config["CpG_from_reference"] == True:
            rescaled_not_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            not_rescaled_not_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            rescaled_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                DP=config["subsampling_depth"],
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            not_rescaled_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                DP=config["subsampling_depth"],
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG)
        elif config["CpG_from_vcf_and_reference"] == True:
            rescaled_not_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            not_rescaled_not_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            rescaled_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                DP=config["subsampling_depth"],
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            not_rescaled_subsampled_CpG = expand("results/gerp/historical/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                DP=config["subsampling_depth"],
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            outlist += (rescaled_not_subsampled_CpG + not_rescaled_not_subsampled_CpG + rescaled_subsampled_CpG + not_rescaled_subsampled_CpG)
    if os.path.exists(config["modern_samples"]):
        not_subsampled_not_CpG = expand("results/gerp/modern/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
            sample=MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,
            minGERP=config["min_gerp"],
            maxGERP=config["max_gerp"],)
        subsampled_not_CpG = expand("results/gerp/modern/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
            sample=MODERN_SUBSAMPLED_NOT_CpG_SAMPLES,
            DP=config["subsampling_depth"],
            fmiss=config["f_missing"],
            chr=CHR,
            minGERP=config["min_gerp"],
            maxGERP=config["max_gerp"],)
        outlist += (not_subsampled_not_CpG + subsampled_not_CpG)
        if config["CpG_from_vcf"] == True:
            not_subsampled_CpG = expand("results/gerp/modern/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            subsampled_CpG = expand("results/gerp/modern/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                DP=config["subsampling_depth"],
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            outlist += (not_subsampled_CpG + subsampled_CpG)
        elif config["CpG_from_reference"] == True:
            not_subsampled_CpG = expand("results/gerp/modern/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            subsampled_CpG = expand("results/gerp/modern/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                DP=config["subsampling_depth"],
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            outlist += (not_subsampled_CpG + subsampled_CpG)
        elif config["CpG_from_vcf_and_reference"] == True:
            not_subsampled_CpG = expand("results/gerp/modern/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            subsampled_CpG = expand("results/gerp/modern/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.Q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
                sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                DP=config["subsampling_depth"],
                fmiss=config["f_missing"],
                chr=CHR,
                minGERP=config["min_gerp"],
                maxGERP=config["max_gerp"],)
            outlist += (not_subsampled_CpG + subsampled_CpG)
    return outlist

def all_GERP_outputs(wildcards):
    """Collect output files for report"""
    if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
        return "results/gerp/all/" + REF_NAME + ".all.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt"
    elif os.path.exists(config["historical_samples"]):
        return "results/gerp/historical/" + REF_NAME + ".historical.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt"
    elif os.path.exists(config["modern_samples"]):
        return "results/gerp/modern/" + REF_NAME + ".modern.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt"


# snakemake rules
localrules:
    fasta_to_fa,
    fna_to_fa,
    split_chunk_bed_files,
    relative_mutational_load_plot,


rule fasta_to_fa:
    """Rename fasta files for outgroup species."""
    input:
        ref=GERP_REF_PATH + "/{gerpref}.fasta.gz",
    output:
        ref=GERP_REF_PATH + "/{gerpref}.fa.gz",
    log:
        "results/logs/13_GERP/fastq/{gerpref}_fasta_to_fa.log",
    shell:
        """
        mv {input} {output} 2> {log}
        """


rule fna_to_fa:
    """Rename fasta files for outgroup species."""
    input:
        ref=GERP_REF_PATH + "/{gerpref}.fna.gz",
    output:
        ref=GERP_REF_PATH + "/{gerpref}.fa.gz",
    log:
        "results/logs/13_GERP/fastq/{gerpref}_fna_to_fa.log",
    shell:
        """
        mv {input} {output} 2> {log}
        """


rule outgroups2fastq:
    """Parse the re-formatted reference into fastq-format with 35 bp read length."""
    input:
        fasta=GERP_REF_PATH + "/{gerpref}.fa.gz",
    output:
        fastq=temp("results/gerp/fastq_files/{gerpref}.fq.gz"),
    log:
        "results/logs/13_GERP/fastq/{gerpref}_outgroups2fastq.log",
    shell:
        """
        python3 workflow/scripts/fa2fq.py {input.fasta} {output.fastq} 2> {log}
        """


rule outgroup_fastqc:
    """Run FastQC on fastq-files from re-formatted outgroup reference genomes."""
    input:
        fastq=rules.outgroups2fastq.output,
    output:
        html="results/gerp/fastq_files/stats/{gerpref}_fastqc.html",
        zip="results/gerp/fastq_files/stats/{gerpref}_fastqc.zip",
        dir=directory("results/gerp/fastq_files/stats/{gerpref}_fastqc/"),
    params:
        dir="results/gerp/fastq_files/stats/",
    log:
        "results/logs/13_GERP/fastq/{gerpref}_outgroup_fastqc.log",
    threads: 2
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    shell:
        """
        fastqc -o {params.dir} -t {threads} --extract {input.fastq} 2> {log}
        """


rule outgroup_fastqc_multiqc:
    """Summarize all fastqc results for re-formatted outgroup reference genomes."""
    input:
        expand("results/gerp/fastq_files/stats/{gerpref}_fastqc.html",
            gerpref=GERP_REF_NAMES,),
        expand("results/gerp/fastq_files/stats/{gerpref}_fastqc.zip",
            gerpref=GERP_REF_NAMES,),
    output:
        stats="results/gerp/fastq_files/stats/multiqc/multiqc_report.html",
    params:
        indir="results/gerp/fastq_files/stats/",
        outdir="results/gerp/fastq_files/stats/multiqc/",
    log:
        "results/logs/13_GERP/fastq/outgroup_fastqc_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule align2target:
    """
    Map fastq reads to the reference genome of the target species.
    Parameters:
        -F 256 -- read maps only to one location
        -F 4 -- keep only mapped reads
        -q 1 -- is it redundant with -F 4? 
        add multi-threading to samtools using -@
    """
    input:
        target=config["ref_path"],
        index=rules.bwa_index_reference.output,
        fastq=rules.outgroups2fastq.output,
        stats=rules.outgroup_fastqc_multiqc.output.stats,
    output:
        bam=temp("results/gerp/alignment/" + REF_NAME + "/{gerpref}.bam"),
    threads: 8
    params:
        extra=r"-R '@RG\tID:{input.fastq}\tSM:{input.fastq}\tPL:ILLUMINA\tPI:330'",
    log:
        "results/logs/13_GERP/alignment/" + REF_NAME + "/{gerpref}_align2target.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        bwa mem {params.extra} -t {threads} {input.target} {input.fastq} | \
            samtools view -@ {threads} -h -q 1 -F 4 -F 256 | grep -v XA:Z | grep -v SA:Z | \
            samtools view -@ {threads} -b - | samtools sort -@ {threads} - > {output.bam} 2> {log}
        """


rule index_gerp_bams:
    """Index bam files."""
    input:
        bam=rules.align2target.output.bam,
    output:
        index=temp("results/gerp/alignment/" + REF_NAME + "/{gerpref}.bam.bai"),
    log:
        "results/logs/13_GERP/alignment/" + REF_NAME + "/{gerpref}_index_gerp_bams.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools index {input.bam} {output.index} 2> {log}
        """


rule gerp_bam_stats:
    """Get some simple statistics on the bam files."""
    input:
        bam=rules.align2target.output.bam,
        index=rules.index_gerp_bams.output.index,
    output:
        stats="results/gerp/alignment/" + REF_NAME + "/stats/{gerpref}.bam.stats.txt",
    log:
        "results/logs/13_GERP/alignment/" + REF_NAME + "/{gerpref}_gerp_bam_stats.log",
    singularity:
        "oras://community.wave.seqera.io/library/bwa_samtools:58df1856e12c14b9"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule gerp_bam_multiqc:
    """Summarize all flagstat reports for gerp bam files."""
    input:
        expand("results/gerp/alignment/" + REF_NAME + "/stats/{gerpref}.bam.stats.txt",
            gerpref=GERP_REF_NAMES,),
    output:
        stats="results/gerp/alignment/" + REF_NAME + "/stats/multiqc/multiqc_report.html",
    params:
        indir="results/gerp/alignment/" + REF_NAME + "/stats/",
        outdir="results/gerp/alignment/" + REF_NAME + "/stats/multiqc/",
    log:
        "results/logs/13_GERP/alignment/" + REF_NAME + "/gerp_bam_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule bam2fasta:
    """
    Parse bam files to fasta format.
    This analysis is run as one job per genome chunk, but is internally run per contig.
    """
    input:
        bam=rules.align2target.output.bam,
        index=rules.index_gerp_bams.output.index,
        stats=rules.gerp_bam_multiqc.output.stats,
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
    output:
        fasta_dir=temp(directory("results/gerp/{chr}_chunks/" + REF_NAME + "/fasta/{gerpref}_{chunk}/")),
    params:
        gerpref="{gerpref}",
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/fasta/{gerpref}_{chunk}_bam2fasta.log",
    threads: 2
    singularity:
        "oras://community.wave.seqera.io/library/samtools_python:2e56d0f345426c81"
    shell:
        """
        if [ ! -d {output.fasta_dir} ]; then
          mkdir -p {output.fasta_dir}
        else
          touch -m {output.fasta_dir}
        fi

        for contig in $(awk -F'\t' '{{print $1}}' {input.chunk_bed}) # run the analysis per contig
        do
          samtools mpileup -aa -r $contig --no-output-ends {input.bam} | python3 workflow/scripts/filter_mpile.py > {output.fasta_dir}/{params.gerpref}_${{contig}}.mpile 2> {log} &&
          python3 workflow/scripts/sequence_to_fastafile.py {output.fasta_dir}/{params.gerpref}_${{contig}}.mpile $contig {params.gerpref} 2>> {log}
        done

        # Check if the fasta files have been created
        for contig in $(awk -F'\t' '{{print $1}}' {input.chunk_bed}) # check each contig
        do
          if [ -f {output.fasta_dir}/{params.gerpref}_${{contig}}.fasta ]; then
            echo "BAM file converted to fasta for" $contig >> {log}
          else
            echo "BAM file conversion to fasta failed for" $contig >> {log} &&
            rm -r {output.fasta_dir} && # Remove the output directory so that Snakemake knows the rule failed
            echo "Removed {output.fasta_dir}" >> {log} &&
            exit 1 # Break the loop so that the Snakemake rule fails 
          fi
        done
        """


rule split_ref_contigs:
    """
    Split the reference genome of the target species into contigs for concatenation with the outgroups.
    Replace contig in output fasta file with reference genome name.
    """
    input:
        ref=config["ref_path"],
        fai=config["ref_path"] + ".fai",
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
        stats=rules.gerp_bam_multiqc.output.stats,
    output:
        fasta_dir=temp(directory("results/gerp/{chr}_chunks/" + REF_NAME + "/fasta/" + REF_NAME + "_{chunk}/")),
    params:
        gerpref=REF_NAME,
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/fasta/" + REF_NAME + "_{chunk}_split_ref_contigs.log",
    threads: 1
    singularity:
        "oras://community.wave.seqera.io/library/seqtk:1.4--e75a8dec899d1be8"
    shell:
        """
        if [ ! -d {output.fasta_dir} ]; then
          mkdir -p {output.fasta_dir};
        fi

        for contig in $(awk -F'\t' '{{print $1}}' {input.chunk_bed}) # run the analysis per contig
        do
          echo $contig > {output.fasta_dir}/${{contig}}.lst &&
          seqtk subseq {input.ref} {output.fasta_dir}/${{contig}}.lst | sed "s/$contig/{params.gerpref}/g" | \
          seqtk seq > {output.fasta_dir}/{params.gerpref}_${{contig}}.fasta 2> {log} &&
          echo $contig "extracted from reference" >> {log}
        done
        """


rule concatenate_fasta_per_contig:
    """
    Concatenate fasta files from all outgroup genomes per contig.
    This analysis is run as one job per genome chunk, but is internally run per contig.
    """
    input:
        gerpref_fasta=expand("results/gerp/{chr}_chunks/" + REF_NAME + "/fasta/{gerpref}_{{chunk}}/", chr=CHR, gerpref=GERP_REF_NAMES),
        ref_fasta=rules.split_ref_contigs.output,
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
    output:
        concatenated_fasta_dir=temp(directory("results/gerp/{chr}_chunks/" + REF_NAME + "/fasta/concatenated_{chunk}/")),
    params:
        fasta_dir="results/gerp/{chr}_chunks/" + REF_NAME + "/fasta",
        chunk=lambda wildcards: "{wildcards.chunk}",
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/fasta/{chunk}_concatenate_fasta_per_contig.log",
    threads: 2
    run:
        if not os.path.exists(output.concatenated_fasta_dir):
            os.makedirs(output.concatenated_fasta_dir)
        chunk_contigs = []
        with open(input.chunk_bed, "r") as bed:
            for line in bed:
                contigname = line.strip().split("\t")[0]
                chunk_contigs.append(contigname)
        for contig in chunk_contigs:  # run the analysis per contig
            ref_fastas = []
            for (ref) in (ALL_GERP_REF_NAMES):  # this list includes the reference genome of the target species
                ref_fastas.append("{fasta_dir}/{ref}_{chunk}/{ref}_{contig}.fasta".format(fasta_dir=params.fasta_dir, ref=ref, chunk=params.chunk, contig=contig,))
            cat_fastas = " ".join(ref_fastas)
            shell("""
                if [ ! -d {{output.concatenated_fasta_dir}} ]; then
                    mkdir -p {{output.concatenated_fasta_dir}};
                fi
                cat {cat_fastas} > {{output.concatenated_fasta_dir}}/{contig}.fasta 2>> {{log}} &&
                echo "Fasta files concatenated for {contig}" >> {{log}}
                """.format(cat_fastas=cat_fastas, contig=contig))


rule compute_gerp:
    """
    Compute GERP++ scores.
    '-s 0.001': tree scaling factor to re-scale tree. Set to 0.001 to scale
    a tree from millions of years (such as trees from www.timetree.org) to 
    billions of years to ensure consistently scaled GERP scores across GenErode
    runs. GERP++ gerpcol default: 1.0
    '-v': verbose mode
    '-a': alignment in mfa format
    Output only includes positions, no contig names.
    This analysis is run as one job per genome chunk, but is internally run per contig.
    """
    input:
        concatenated_fasta_dir=rules.concatenate_fasta_per_contig.output,
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
        tree=config["tree"],
    output:
        gerp_dir=temp(directory("results/gerp/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}_gerp_raw/")),
    params:
        name=REF_NAME,
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}_compute_gerp.log",
    threads: 4
    singularity:
        "https://depot.galaxyproject.org/singularity/gerp:2.1--h1b792b2_2"
    shell:
        """
        if [ ! -d {output.gerp_dir} ]; then 
            mkdir -p {output.gerp_dir}; 
        fi
        for contig in $(awk -F'\t' '{{print $1}}' {input.chunk_bed}) # run the analysis per contig
        do
          gerpcol -v -f {input.concatenated_fasta_dir}/${{contig}}.fasta -t {input.tree} -a -s 0.001 -e {params.name} 2> {log} &&
          mv {input.concatenated_fasta_dir}/${{contig}}.fasta.rates {output.gerp_dir} 2>> {log} &&
          echo "Computed GERP++ scores for" $contig >> {log} 
        done
        """


rule gerp2coords:
    """
    Convert GERP-scores to the correct genomic coordinates. 
    Script currently written to output positions without contig names.
    This analysis is run as one job per genome chunk, but is internally run per contig.
    """
    input:
        concatenated_fasta_dir=rules.concatenate_fasta_per_contig.output,
        gerp_dir=rules.compute_gerp.output,
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
    output:
        gerp_coords_dir=temp(directory("results/gerp/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}_gerp_coords/")),
    params:
        name=REF_NAME,
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}_gerp2coords.log",
    threads: 2
    run:
        chunk_contigs = []
        with open(input.chunk_bed, "r") as file:
            for line in file:
                contigname = line.strip().split("\t")[0]
                chunk_contigs.append(contigname)
        for contig in chunk_contigs:  # run the analysis per contig
            shell("""
                if [ ! -d {{output.gerp_coords_dir}} ]; then 
                    mkdir -p {{output.gerp_coords_dir}}; 
                fi
                python3 workflow/scripts/gerp_to_position.py {{input.concatenated_fasta_dir}}/{contig}.fasta \
                {{input.gerp_dir}}/{contig}.fasta.rates {{params.name}} 2>> {{log}} &&
                mv {{input.gerp_dir}}/{contig}.fasta.rates.parsed {{output.gerp_coords_dir}} 2>> {{log}} &&
                echo "GERP-score coordinates converted for {contig}" >> {{log}}
                """.format(contig=contig))


rule get_ancestral_state:
    """Get the ancestral state of each position in the focal reference genome."""
    input:
        concatenated_fasta_dir=rules.concatenate_fasta_per_contig.output,
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
    output:
        fasta_ancestral_dir=temp(directory("results/gerp/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}_fasta_ancestral/")),
    params:
        name=REF_NAME,
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}_get_ancestral_state.log",
    threads: 2
    run:
        chunk_contigs = []
        with open(input.chunk_bed, "r") as file:
            for line in file:
                contigname = line.strip().split("\t")[0]
                chunk_contigs.append(contigname)
        for contig in chunk_contigs:  # run the analysis per contig
            shell("""
                if [ ! -d {{output.fasta_ancestral_dir}} ]; then 
                    mkdir -p {{output.fasta_ancestral_dir}}; 
                fi
                python3 workflow/scripts/fasta_to_major_allele.py {{input.concatenated_fasta_dir}}/{contig}.fasta {{params.name}} 2>> {{log}} &&
                mv {{input.concatenated_fasta_dir}}/{contig}.fasta.parsed {{output.fasta_ancestral_dir}} 2>> {{log}} &&
                echo "Ancestral states obtained for {contig}" >> {{log}}
                """.format(contig=contig))


rule produce_contig_out:
    """Merge the ancestral allele and gerp-scores into one file per contig."""
    input:
        fasta_ancestral_dir=rules.get_ancestral_state.output.fasta_ancestral_dir,
        gerp_coords_dir=rules.gerp2coords.output.gerp_coords_dir,
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
    output:
        gerp_merged_dir=temp(directory("results/gerp/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}_gerp_merged/")),
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}_produce_contig_out.log",
    threads: 2
    run:
        chunk_contigs = []
        with open(input.chunk_bed, "r") as file:
            for line in file:
                contigname = line.strip().split("\t")[0]
                chunk_contigs.append(contigname)
        for contig in chunk_contigs:  # run the analysis per contig
            shell("""
                if [ ! -d {{output.gerp_merged_dir}} ]; then 
                    mkdir -p {{output.gerp_merged_dir}}; 
                fi
                paste {{input.fasta_ancestral_dir}}/{contig}.fasta.parsed {{input.gerp_coords_dir}}/{contig}.fasta.rates.parsed | \
                sed "s/^/{contig}\t/g" > {{output.gerp_merged_dir}}/{contig}.fasta.parsed.rates 2>> {{log}} &&
                echo "GERP-scores and ancestral states merged for {contig}" >> {{log}}
                """.format(contig=contig))


rule merge_gerp_per_chunk:
    """Merge results per genome chunk into one file."""
    input:
        gerp_merged_dir=rules.produce_contig_out.output.gerp_merged_dir,
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
    output:
        gerp_chunks_merged=temp("results/gerp/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}.fasta.parsed.rates"),
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}_merge_per_chunk.log",
    threads: 2
    run:
        chunk_contigs = []
        with open(input.chunk_bed, "r") as file:
            for line in file:
                contigname = line.strip().split("\t")[0]
                chunk_contigs.append(contigname)
        contig_rates = []  # collect all file paths for concatenation in a list
        for contig in chunk_contigs:
            contig_rates.append("{gerp_merged_dir}/{contig}.fasta.parsed.rates".format(gerp_merged_dir=input.gerp_merged_dir, contig=contig))
        cat_rates = " ".join(contig_rates)  # concatenate list to a string with whitespace as separator, as input for cat command
        shell("""cat {cat_rates} > {{output.gerp_chunks_merged}} 2>> {{log}}""".format(cat_rates=cat_rates))


rule merge_gerp_gz:
    """Merge results per contig into one file."""
    input:
        gerp_chunks_merged=expand("results/gerp/{chr}_chunks/" + REF_NAME + "/gerp/{chunk}.fasta.parsed.rates", chr=CHR, chunk=CHUNKS),
    output:
        gerp_out="results/gerp/" + REF_NAME + ".{chr}.ancestral.rates.gz",
    log:
        "results/logs/13_GERP/" + REF_NAME + ".{chr}.merge_gerp_gz.log",
    threads: 2
    shell:
        """
        cat {input.gerp_chunks_merged} | gzip - > {output.gerp_out} 2> {log}
        """


rule plot_gerp_hist:
    """Plot the GERP scores as histogram"""
    input:
        gerp_out=rules.merge_gerp_gz.output.gerp_out,
    output:
        pdf=report("results/gerp/" + REF_NAME + ".{chr}.ancestral.rates.gerp.hist.pdf",
            caption="../report/gerp_plot.rst",
            category="GERP",),
    threads: 2
    log:
        "results/logs/13_GERP/" + REF_NAME + ".{chr}.plot_gerp_hist.log",
    script:
        "../scripts/gerp_hist_plot.py"


rule split_vcf_files:
    """Split individual VCF files into chunks for more resource-efficient merging with GERP results"""
    input:
        vcf=rules.filter_biallelic_missing_vcf.output.filtered,
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
        genomefile=REF_DIR + "/" + REF_NAME + ".genome",
    output:
        vcf_chunk=temp("results/gerp/{chr}_chunks/" + REF_NAME + "/{dataset}/vcf/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.{chunk}.vcf.gz"),
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/{dataset}/vcf/{sample}.{processed}_fmissing{fmiss}.{chr}.{chunk}_split_vcf_chunks.log",
    singularity:
        "oras://community.wave.seqera.io/library/bedtools_htslib:06ed4722f423d939"
    shell:
        """
        bedtools intersect -a {input.vcf} -b {input.chunk_bed} -g {input.genomefile} -header | gzip - > {output.vcf_chunk} 2> {log}
        """


rule split_chunk_bed_files:
    """Split the chunk bed files into 10 million basepair windows"""
    input:
        chunk_bed=REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/{chunk}.bed",
    output:
        chunk_win_bed=temp(REF_DIR + "/gerp/" + REF_NAME + "/split_bed_files_{chr}/windows/{chunk}_10Mwindows.bed"),
    log:
        "results/logs/13_GERP/" + REF_NAME + ".{chunk}_{chr}_split_chunk_bed_files.log",
    singularity:
        "oras://community.wave.seqera.io/library/bedtools_htslib:06ed4722f423d939"
    shell:
        """
        bedtools makewindows -b {input.chunk_bed} -w 10000000 > {output.chunk_win_bed} 2> {log}
        """


rule gerp_derived_alleles:
    """
    Add the number of derived alleles per sample and position to the gerp output for each window, run separately for each sample.
    """
    input:
        gerp_out=rules.merge_gerp_per_chunk.output.gerp_chunks_merged,
        gerp_merged_dir=rules.produce_contig_out.output.gerp_merged_dir,
        vcf=rules.split_vcf_files.output.vcf_chunk,
        chunk_win_bed=rules.split_chunk_bed_files.output.chunk_win_bed,
    output:
        gerp_alleles_dir=temp(directory("results/gerp/{chr}_chunks/" + REF_NAME + "/{dataset}/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.{chunk}_gerp_derived_alleles/")),
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/{dataset}/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.{chunk}_gerp_derived_alleles.log",
    threads: 4
    shell:
        """
        if [ ! -d {output.gerp_alleles_dir} ]; then 
            mkdir -p {output.gerp_alleles_dir}; 
        fi
        while IFS="\t" read -r contig start end; do 
            python3 workflow/scripts/gerp_derived_alleles.py {input.gerp_merged_dir}/${{contig}}.fasta.parsed.rates \
            {input.vcf} ${{contig}} ${{start}} ${{end}} {output.gerp_alleles_dir}/${{contig}}_${{start}}_${{end}}.fasta.parsed.rates.derived_alleles;
        done < {input.chunk_win_bed} 2>> {log}
        """


rule merge_gerp_alleles_per_chunk:
    """Merge results per genome chunk into one file per sample."""
    input:
        gerp_alleles_dir=rules.gerp_derived_alleles.output.gerp_alleles_dir,
        chunk_win_bed=rules.split_chunk_bed_files.output.chunk_win_bed,
    output:
        gerp_chunks_merged=temp("results/gerp/{chr}_chunks/" + REF_NAME + "/{dataset}/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.{chunk}.fasta.parsed.rates.derived_alleles"),
    log:
        "results/logs/13_GERP/{chr}_chunks/" + REF_NAME + "/{dataset}/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.{chunk}_merge_gerp_alleles_per_chunk.log",
    threads: 4
    run:
        chunk_windows = []
        with open(input.chunk_win_bed, "r") as file:
            for line in file:
                line_tup = ()
                contigname = line.strip().split("\t")[0]
                start = line.strip().split("\t")[1]
                end = line.strip().split("\t")[2]
                line_tup = (contigname, start, end)
                chunk_windows.append(line_tup)
        window_rates_alleles = []  # collect all file paths for concatenation in a list
        for window in chunk_windows:
            window_rates_alleles.append("{gerp_alleles_dir}/{contig}_{start}_{end}.fasta.parsed.rates.derived_alleles".format(
                gerp_alleles_dir=input.gerp_alleles_dir, contig=window[0], start=window[1], end=window[2]))
        cat_rates_alleles = " ".join(window_rates_alleles)  # concatenate list to a string with whitespace as separator, as input for awk command
        shell("""awk 'FNR>1 || NR==1' {cat_rates_alleles} > {{output.gerp_chunks_merged}} 2>> {{log}}""".format(cat_rates_alleles=cat_rates_alleles))  # awk command concatenates all files without their header line, except the first file


rule merge_gerp_alleles_gz:
    """Merge results into one file per sample."""
    input:
        gerp_chunks_merged=expand("results/gerp/{chr}_chunks/" + REF_NAME + "/{{dataset}}/{{sample}}.merged.rmdup.merged.{{processed}}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.{chunk}.fasta.parsed.rates.derived_alleles",
            chunk=CHUNKS,
            fmiss=config["f_missing"],
            chr=CHR,
            min_gerp=config["min_gerp"],
            max_gerp=config["max_gerp"],),
    output:
        gerp_out="results/gerp/{dataset}/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.ancestral.rates.derived.alleles.gz",
    log:
        "results/logs/13_GERP/{dataset}/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.merge_gerp_alleles_gz.log",
    threads: 4
    shell:
        """
        awk 'FNR>1 || NR==1' {input.gerp_chunks_merged} | gzip - > {output.gerp_out} 2> {log}
        """


rule relative_mutational_load_per_sample:
    """Calculate the relative mutational load per sample."""
    input:
        gerp_out="results/gerp/{dataset}/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.ancestral.rates.derived.alleles.gz",
    output:
        mut_load=temp("results/gerp/{dataset}/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt"),
    params:
        min_gerp=config["min_gerp"],
        max_gerp=config["max_gerp"],
    log:
        "results/logs/13_GERP/{dataset}/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.relative_mutational_load_table.gerp_{minGERP}_{maxGERP}.log",
    threads: 2
    shell:
        """
        python3 workflow/scripts/gerp_rel_mut_load_sample.py {input.gerp_out} {params.min_gerp} {params.max_gerp} {output.mut_load} 2> {log}
        """


rule relative_mutational_load_table:
    """Combine the relative mutational load across all samples."""
    input:
        rel_load_table_inputs,
    output:
        mut_load=report("results/gerp/{dataset}/" + REF_NAME + ".{dataset}.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_table.txt",
            caption="../report/relative_mutational_load_table.rst",
            category="GERP",),
    log:
        "results/logs/13_GERP/{dataset}/" + REF_NAME + ".{dataset}.fmissing{fmiss}.{chr}.relative_mutational_load_table.gerp_{minGERP}_{maxGERP}.log",
    threads: 4
    script:
        "../scripts/gerp_rel_mut_load_table.py"


rule relative_mutational_load_plot:
    """Plot the relative mutational load for modern and historical samples"""
    input:
        all_GERP_outputs,
    output:
        plot=report("results/gerp/{dataset}/" + REF_NAME + ".{dataset}.fmissing{fmiss}.{chr}.relative_mutational_load.gerp_{minGERP}_{maxGERP}_plot.pdf",
            caption="../report/relative_mutational_load_plot.rst",
            category="GERP",),
    log:
        "results/logs/13_GERP/{dataset}/" + REF_NAME + ".{dataset}.fmissing{fmiss}.{chr}.relative_mutational_load_plot.gerp_{minGERP}_{maxGERP}.log",
    script:
        "../scripts/gerp_rel_mut_load_plot.py"
