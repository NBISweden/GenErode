##########################################################################
### 1.2 OPTIONAL Mapping of historical reads to human and animal mitochondrial genomes

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    all_outputs.append("results/historical/mitogenomes_mapping/stats/multiqc/multiqc_report.html")


# Functions used by rules of this part of the pipeline
def get_mapped_reads(inlist):
    """
    Parse the samtools flagstat files from a python list to return a 
    dictionary containing the filename and the number of mapped reads
    """
    readsDict = {}
    for i in inlist:
        bn = os.path.basename(i)
        name = "_".join(bn.split("_", 3)[:3])
        with open(i, "r") as f:
            for line in f:
                if " mapped (" in line:
                    reads = int(line.strip().split()[0])
                    if reads == 0:
                        readsDict[name] = reads + 1
                    else:
                        readsDict[name] = reads
    return readsDict

def ratio_of_mapped_reads(readsDictOtherSpecies, readsDictTargetSpecies):
    """
    Take two dictionaries containing the filename and the number of mapped reads for the same set of reads 
    mapped to a different species and to the target species and return the ratio of mapped reads
    """
    ratio = { k: round(readsDictOtherSpecies[k] / readsDictTargetSpecies[k], 2) for k in readsDictTargetSpecies if k in readsDictOtherSpecies }
    return ratio

def historical_mapped_reads_ratios_species_inputs(wildcards):
    """Input for mapped_reads_ratios"""
    if config["map_unmerged_reads"] == True:
        species = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_" + MITO_NAME + ".sorted.bam.stats.txt",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged", "unmerged"],)
    elif config["map_unmerged_reads"] == False:
        species = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_" + MITO_NAME + ".sorted.bam.stats.txt",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged"],)
    return species

def historical_mapped_reads_ratios_other_inputs(wildcards):
    """Input for mapped_reads_ratios"""
    if config["map_unmerged_reads"] == True:
        other = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_{{mitoref}}.sorted.bam.stats.txt",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged", "unmerged"],)
    elif config["map_unmerged_reads"] == False:
        other = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_{{mitoref}}.sorted.bam.stats.txt",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged"],)
    return other

def historical_mito_bams_merge_files_inputs(wildcards):
    """
    Input for merge_read_ratio_files depending if 
    only merged or also unmerged reads should be mapped
    """
    if config["map_unmerged_reads"] == True:
        ratiotables = expand("results/historical/mitogenomes_mapping/stats/{mitoref}_vs_" + MITO_NAME + "_{reads}.txt",
            reads=["merged", "unmerged"],
            mitoref=[HUMAN_NAME, CHICK_NAME, COW_NAME, PIG_NAME, MOUSE_NAME],)
    elif config["map_unmerged_reads"] == False:
        ratiotables = expand("results/historical/mitogenomes_mapping/stats/{mitoref}_vs_" + MITO_NAME + "_{reads}.txt",
            reads=["merged"],
            mitoref=[HUMAN_NAME, CHICK_NAME, COW_NAME, PIG_NAME, MOUSE_NAME],)
    return ratiotables

def merge_hist_mito_bams_per_sample_inputs(wildcards):
    """Input for merge_historical_mitogenome_bams_per_sample"""
    SAMPLE = "{}".format(wildcards.sample)
    SAMPLEIDXLN_LIST = hist_sample_dict[SAMPLE]
    return expand("results/historical/mitogenomes_mapping/{sampleindexlane}_merged_{{mitoref}}.sorted.bam",
        sampleindexlane=SAMPLEIDXLN_LIST,)

def historical_mito_bams_multiqc_inputs(wildcards):
    """Input for historical_mito_bams_multiqc and all"""
    if config["map_unmerged_reads"] == True:
        quali_report = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_{mitoref}.sorted.bam.qualimap/qualimapReport.html",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged", "unmerged"],
            mitoref=[MITO_NAME, HUMAN_NAME, CHICK_NAME, COW_NAME, PIG_NAME, MOUSE_NAME],)
        summary = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_{mitoref}.sorted.bam.qualimap/genome_results.txt",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged", "unmerged"],
            mitoref=[MITO_NAME, HUMAN_NAME, CHICK_NAME, COW_NAME, PIG_NAME, MOUSE_NAME],)
        stats = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_{mitoref}.sorted.bam.stats.txt",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged", "unmerged"],
            mitoref=[MITO_NAME, HUMAN_NAME, CHICK_NAME, COW_NAME, PIG_NAME, MOUSE_NAME],)
        ratios = ["results/historical/mitogenomes_mapping/stats/ratios_of_merged_and_unmerged_mapped_to_various_mitochondrial_genomes_vs_" + MITO_NAME + ".txt"]
        merged_bam_stats = expand("results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.stats.txt",
            sample=hist_sm,
            mitoref=[MITO_NAME],)
        merged_bam_quali_report = expand("results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.qualimap/qualimapReport.html",
            sample=hist_sm,
            mitoref=[MITO_NAME],)
        merged_bam_summary = expand("results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.qualimap/genome_results.txt",
            sample=hist_sm,
            mitoref=[MITO_NAME],)
        outlist = (quali_report + summary + stats + ratios + merged_bam_stats + merged_bam_quali_report + merged_bam_summary)
    elif config["map_unmerged_reads"] == False:
        quali_report = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_{mitoref}.sorted.bam.qualimap/qualimapReport.html",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged"],
            mitoref=[MITO_NAME, HUMAN_NAME, CHICK_NAME, COW_NAME, PIG_NAME, MOUSE_NAME],)
        summary = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_{mitoref}.sorted.bam.qualimap/genome_results.txt",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged"],
            mitoref=[MITO_NAME, HUMAN_NAME, CHICK_NAME, COW_NAME, PIG_NAME, MOUSE_NAME],)
        stats = expand("results/historical/mitogenomes_mapping/stats/{sampleindexlane}_{reads}_{mitoref}.sorted.bam.stats.txt",
            sampleindexlane=hist_sm_idx_ln,
            reads=["merged"],
            mitoref=[MITO_NAME, HUMAN_NAME, CHICK_NAME, COW_NAME, PIG_NAME, MOUSE_NAME],)
        ratios = ["results/historical/mitogenomes_mapping/stats/ratios_of_merged_mapped_to_various_mitochondrial_genomes_vs_" + MITO_NAME + ".txt"]
        merged_bam_stats = expand("results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.stats.txt",
            sample=hist_sm,
            mitoref=[MITO_NAME],)
        merged_bam_quali_report = expand("results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.qualimap/qualimapReport.html",
            sample=hist_sm,
            mitoref=[MITO_NAME],)
        merged_bam_summary = expand("results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.qualimap/genome_results.txt",
            sample=hist_sm,
            mitoref=[MITO_NAME],)
        outlist = (quali_report + summary + stats + ratios + merged_bam_stats + merged_bam_quali_report + merged_bam_summary)
    return outlist


# snakemake rules
localrules: symlink_mito_genome, merge_read_ratio_files

rule symlink_mito_genome:
    """Create a symbolic link of the mitochondrial genome fasta file of the target species"""
    input:
        ref=config["species_mt_path"],
    output:
        ref=MITO_DIR + MITO_NAME + ".fasta",
    params:
        abs_in=lambda wildcards, input: os.path.abspath(input.ref),
        abs_out=lambda wildcards, output: os.path.abspath(output.ref),
    log:
        "results/logs/1.2_map_to_mitogenomes/symlink_mito_genome.log",
    shell:
        """
        ln -s {params.abs_in} {params.abs_out} 2> {log}
        """


rule bwa_index_mito_ref:
    """Index the mitogenomes using bwa"""
    input:
        ref=MITO_DIR + "{mitoref}.fasta",
    output:
        amb=MITO_DIR + "{mitoref}.fasta.amb",
        ann=MITO_DIR + "{mitoref}.fasta.ann",
        bwt=MITO_DIR + "{mitoref}.fasta.bwt",
        pac=MITO_DIR + "{mitoref}.fasta.pac",
        sa=MITO_DIR + "{mitoref}.fasta.sa",
    log:
        "results/logs/1.2_map_to_mitogenomes/{mitoref}_bwa_index_mito_ref.log",
    singularity:
        "docker://biocontainers/bwa:v0.7.17-3-deb_cv1"
    shell:
        """
        bwa index {input.ref} 2> {log}
        """


rule map_historical_merged_to_mito:
    """Map merged reads from historical samples to mitochondrial genomes."""
    """BWA aln for short Illumina reads, parameters according to Palkopoulou et al. 2015"""
    input:
        ref=MITO_DIR + "{mitoref}.fasta",
        index=rules.bwa_index_mito_ref.output,
        merged="results/historical/trimming/{sample}_{index}_{lane}_trimmed_merged.fastq.gz",
    output:
        bam=temp("results/historical/mitogenomes_mapping/{sample}_{index}_{lane}_merged_{mitoref}.sorted.bam"),
    log:
        "results/logs/1.2_map_to_mitogenomes/{sample}_{index}_{lane}_{mitoref}_map_historical_merged_to_mito.log",
    singularity:
        "docker://nbisweden/generode-bwa:latest"
    shell:
        """
        bwa aln -l 16500 -n 0.01 -o 2 -t {threads} {input.ref} {input.merged} | \
        bwa samse {input.ref} - {input.merged} | samtools sort - > {output.bam} 2> {log}
        """


rule map_historical_unmerged_to_mito:
    """Map unmerged reads from historical samples to mitochondrial genomes."""
    """BWA aln for short Illumina reads, parameters according to Palkopoulou et al. 2015"""
    input:
        ref=MITO_DIR + "{mitoref}.fasta",
        index=rules.bwa_index_mito_ref.output,
        R1_un="results/historical/trimming/{sample}_{index}_{lane}_R1_unmerged.fastq.gz",
        R2_un="results/historical/trimming/{sample}_{index}_{lane}_R2_unmerged.fastq.gz",
    output:
        R1_sai=temp("results/historical/mitogenomes_mapping/{sample}_{index}_{lane}_R1_{mitoref}.sai"),
        R2_sai=temp("results/historical/mitogenomes_mapping/{sample}_{index}_{lane}_R2_{mitoref}.sai"),
        bam=temp("results/historical/mitogenomes_mapping/{sample}_{index}_{lane}_unmerged_{mitoref}.sorted.bam"),
    log:
        "results/logs/1.2_map_to_mitogenomes/{sample}_{index}_{lane}_{mitoref}_map_historical_unmerged_to_mito.log",
    singularity:
        "docker://nbisweden/generode-bwa:latest"
    shell:
        """
        bwa aln -l 16500 -n 0.01 -o 2 -t {threads} {input.ref} {input.R1_un} > {output.R1_sai} 2> {log} &&
        bwa aln -l 16500 -n 0.01 -o 2 -t {threads} {input.ref} {input.R2_un} > {output.R2_sai} 2>> {log} &&
        bwa samse {input.ref} {output.R1_sai} {output.R2_sai} {input.R1_un} {input.R2_un} | samtools sort - > {output.bam} 2>> {log}
        """


rule mitogenome_bam_stats:
    """Basic stats of bam files"""
    input:
        bam="results/historical/mitogenomes_mapping/{sample}_{index}_{lane}_{reads}_{mitoref}.sorted.bam",
    output:
        stats="results/historical/mitogenomes_mapping/stats/{sample}_{index}_{lane}_{reads}_{mitoref}.sorted.bam.stats.txt",
    group:
        "historical_mito_bams_group"
    log:
        "results/logs/1.2_map_to_mitogenomes/{sample}_{index}_{lane}_{reads}_{mitoref}_mitogenome_bam_stats.log",
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule historical_mito_bams_qualimap:
    """Get stats of mapping success to mitogenomes"""
    input:
        bam="results/historical/mitogenomes_mapping/{sample}_{index}_{lane}_{reads}_{mitoref}.sorted.bam",
        stats="results/historical/mitogenomes_mapping/stats/{sample}_{index}_{lane}_{reads}_{mitoref}.sorted.bam.stats.txt",
    output:
        report="results/historical/mitogenomes_mapping/stats/{sample}_{index}_{lane}_{reads}_{mitoref}.sorted.bam.qualimap/qualimapReport.html",
        summary="results/historical/mitogenomes_mapping/stats/{sample}_{index}_{lane}_{reads}_{mitoref}.sorted.bam.qualimap/genome_results.txt",
    params:
        outdir="results/historical/mitogenomes_mapping/stats/{sample}_{index}_{lane}_{reads}_{mitoref}.sorted.bam.qualimap/",
    log:
        "results/logs/1.2_map_to_mitogenomes/{sample}_{index}_{lane}_{reads}_{mitoref}_historical_mito_bams_qualimap.log",
    group:
        "historical_mito_bams_group"
    threads: 1
    singularity:
        "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
    shell:
        """
        reads=`head -n1 {input.stats} | cut -d' ' -f 1`
        mem=$(((6 * {threads}) - 2))
        if [ "$reads" -gt 100 ] # check if bam file contains enough reads
        then
          unset DISPLAY
          qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G -nt {threads} -nr 100 -outdir {params.outdir} -outformat html 2> {log}
        else
          mkdir -p {params.outdir} 2> {log}
          touch {output.report} 2>> {log}
          touch {output.summary} 2>> {log}
          echo "Not enough reads to run QualiMap" >> {log}
        fi
        """


rule mapped_reads_ratios:
    """Get the ratios of reads mapped to other species and mapped to target species to identify outliers that might indicate contamination"""
    input:
        species=historical_mapped_reads_ratios_species_inputs,
        other=historical_mapped_reads_ratios_other_inputs,
    output:
        ratiotables=temp("results/historical/mitogenomes_mapping/stats/{mitoref}_vs_" + MITO_NAME + "_{reads}.txt"),
    log:
        "results/logs/1.2_map_to_mitogenomes/{mitoref}_vs_" + MITO_NAME + "_{reads}_mapped_reads_ratios.log",
    run:
        with open(output.ratiotables, "w") as out:
            # Parse the samtools flagstat files from the different input lists to return a dictionary containing the filename and the number of mapped reads for each of them
            speciesDict = get_mapped_reads(input.species)
            otherDict = get_mapped_reads(input.other)

            # Calculate the ratios of mapped reads for each species to the target species reference
            other_species = ratio_of_mapped_reads(otherDict, speciesDict)

            # Write to output file
            for k, v in other_species.items():
                out.write("{} {}\n".format(k, v))


rule merge_read_ratio_files:
    """Merge all read ratio files into one, taking information from file name as first column"""
    input:
        historical_mito_bams_merge_files_inputs,
    output:
        finaltable="results/historical/mitogenomes_mapping/stats/ratios_of_{reads}_mapped_to_various_mitochondrial_genomes_vs_" + MITO_NAME + ".txt",
    log:
        "results/logs/1.2_map_to_mitogenomes/{reads}_merge_read_ratio_files.log",
    shell:
        """
        echo "reference_comparison sample ratio_of_mapped_reads" > {output.finaltable} && 
        awk '{{print FILENAME, $0}}' {input} | sed 's/_*..merged.txt//g' | awk -F'/' '{{print $NF}}' >> {output.finaltable} 2> {log}
        """


rule merge_historical_mitogenome_bams_per_sample:
    """Merge files per sample"""
    input:
        merge_hist_mito_bams_per_sample_inputs,
    output:
        merged="results/historical/mitogenomes_mapping/{sample}_merged_{mitoref}.sorted.merged.bam",
    message:
        "the input files are: {input}"
    log:
        "results/logs/1.2_map_to_mitogenomes/{sample}_{mitoref}_merge_historical_mitogenome_bams_per_sample.log",
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        """
        files=`echo {input} | awk '{{print NF}}'`
        if [ $files -gt 1 ] # check if there are at least 2 files for merging. If there is only one file, copy the sorted bam file.
        then
          samtools merge {output.merged} {input} 2> {log}
        else
          cp {input} {output.merged} && touch {output.merged} 2> {log}
          echo "Only one file present for merging. Copying the sorted bam file." >> {log}
        fi
        """


rule merged_mitogenome_bam_stats:
    """Basic stats of bam files"""
    input:
        bam=rules.merge_historical_mitogenome_bams_per_sample.output.merged,
    output:
        stats="results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.stats.txt",
    log:
        "results/logs/1.2_map_to_mitogenomes/{sample}_{mitoref}_merged_mitogenome_bam_stats.log",
    group:
        "historical_merged_mito_bams_group"
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        """
        samtools flagstat {input.bam} > {output.stats} 2> {log}
        """


rule historical_merged_mito_bams_qualimap:
    """Get stats of mapping success to mitogenomes"""
    input:
        bam=rules.merge_historical_mitogenome_bams_per_sample.output.merged,
        stats="results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.stats.txt",
    output:
        report="results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.qualimap/qualimapReport.html",
        summary="results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.qualimap/genome_results.txt",
        outdir=directory("results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.qualimap/"),
    params:
        outdir="results/historical/mitogenomes_mapping/stats/{sample}_merged_{mitoref}.sorted.merged.bam.qualimap/",
    log:
        "results/logs/1.2_map_to_mitogenomes/{sample}_{mitoref}_historical_merged_mito_bams_qualimap.log",
    group:
        "historical_merged_mito_bams_group"
    threads: 1
    singularity:
        "oras://community.wave.seqera.io/library/qualimap:2.3--95d781b369b835f2"
    shell:
        """
        reads=`head -n1 {input.stats} | cut -d' ' -f 1`
        mem=$(((6 * {threads}) - 2))
        if [ "$reads" -gt 100 ] # check if bam file contains enough reads
        then
          unset DISPLAY
          qualimap bamqc -bam {input.bam} --java-mem-size=${{mem}}G -nt {threads} -nr 100 -outdir {params.outdir} -outformat html 2> {log}
        else
          mkdir -p {params.outdir} 2> {log}
          touch {output.report} 2>> {log}
          touch {output.summary} 2>> {log}
          echo "Not enough reads to run QualiMap" >> {log}
        fi
        """


rule historical_mito_bams_multiqc:
    """Summarize qualimap results for all samples"""
    input:
        historical_mito_bams_multiqc_inputs,
    output:
        "results/historical/mitogenomes_mapping/stats/multiqc/multiqc_report.html",
    params:
        indir="results/historical/mitogenomes_mapping/stats/",
        outdir="results/historical/mitogenomes_mapping/stats/multiqc/",
    log:
        "results/logs/1.2_map_to_mitogenomes/historical_mito_bams_multiqc.log",
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """
