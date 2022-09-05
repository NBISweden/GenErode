##########################################################################
### 0.1 Preparation of the reference genome for downstream analyses

# Code collecting output files from this part of the pipeline
all_outputs.append(expand("{ref_path}.{ext}", 
    ref_path=config["ref_path"], 
    ext=["amb", "ann", "bwt", "pac", "sa", "fai"]))
all_outputs.append(expand(REF_DIR + "/" + REF_NAME + ".{ext}", 
    ext=["dict", "genome", "bed"]))

# snakemake rules
localrules: make_reference_bed

rule bwa_index_reference:
    """Index the reference genome using bwa"""
    input:
        ref=config["ref_path"],
    output:
        amb=config["ref_path"] + ".amb",
        ann=config["ref_path"] + ".ann",
        bwt=config["ref_path"] + ".bwt",
        pac=config["ref_path"] + ".pac",
        sa=config["ref_path"] + ".sa",
    params:
        dir=REF_DIR,
    log:
        "results/logs/0.1_reference_genome_preps/" + REF_NAME + "_bwa_index.log",
    group:
        "reference_prep_group"
    singularity:
        "docker://biocontainers/bwa:v0.7.17-3-deb_cv1"
    shell:
        """
        bwa index -a bwtsw {input.ref} 2> {log}
        """


rule samtools_fasta_index:
    """Index the reference genome using samtools"""
    input:
        ref=config["ref_path"],
    output:
        fai=config["ref_path"] + ".fai",
    log:
        "results/logs/0.1_reference_genome_preps/"
        + REF_NAME
        + "_samtools_fasta_index.log",
    group:
        "reference_prep_group"
    singularity:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        """
        samtools faidx {input.ref} 2> {log}
        """


rule picard_fasta_dict:
    """Index the reference genome with picard for GATK tools"""
    input:
        ref=config["ref_path"],
    output:
        fdict=REF_DIR + "/" + REF_NAME + ".dict",
    params:
        mem="4g",
    log:
        "results/logs/0.1_reference_genome_preps/" + REF_NAME + "_picard_fasta_dict.log",
    group:
        "reference_prep_group"
    singularity:
        "docker://quay.io/biocontainers/picard:2.26.6--hdfd78af_0"
    shell:
        """
        picard CreateSequenceDictionary -Xmx{params.mem} R={input.ref} O={output.fdict} 2> {log}
        """


rule genome_file:
    """Create a genome file for filtering of VCF and BAM files"""
    input:
        fai=rules.samtools_fasta_index.output,
    output:
        genomefile=REF_DIR + "/" + REF_NAME + ".genome",
    log:
        "results/logs/0.1_reference_genome_preps/" + REF_NAME + "_genome_file.log",
    group:
        "reference_prep_group"
    shell:
        """
        awk -v OFS='\t' '{{print $1, $2}}' {input.fai} > {output.genomefile} 2> {log}
        """


rule make_reference_bed:
    """Generate a BED file with all genome regions"""
    input:
        fai=rules.samtools_fasta_index.output,
    output:
        ref_bed=REF_DIR + "/" + REF_NAME + ".bed",
    log:
        "results/logs/0.1_reference_genome_preps/"
        + REF_NAME
        + "_make_reference_bed.log",
    shell:
        """
        awk -v OFS='\t' '{{print $1, "0", $2}}' {input.fai} > {output.ref_bed} 2> {log}
        """
