##########################################################################
### 0.1 Preparation of the reference genome for downstream analyses

# Code collecting output files from this part of the pipeline
ref_prep_outputs=[]
ref_prep_outputs.append(expand("{ref_path}.{ext}", 
    ref_path=config["ref_path"], 
    ext=["amb", "ann", "bwt", "pac", "sa", "fai"]))
ref_prep_outputs.append(expand(REF_DIR + "/" + REF_NAME + ".{ext}", 
    ext=["dict", "bed"]))

# snakemake rules

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
    container:
        bwa_container
    shell:
        """
        bwa index -a bwtsw {input.ref} 2> {log}
        """


rule index_reference:
    """Index the reference genome using samtools and picard"""
    input:
        ref=config["ref_path"],
    output:
        fai=config["ref_path"] + ".fai",
        ref_bed=REF_DIR + "/" + REF_NAME + ".bed",
        fdict=REF_DIR + "/" + REF_NAME + ".dict",
    params:
        mem="4g",
    log:
        "results/logs/0.1_reference_genome_preps/" + REF_NAME + "_index_reference.log",
    container:
        bwa_picard_samtools_container
    shell:
        """
        samtools faidx {input.ref} 2> {log} &&
        awk -v OFS='\t' '{{print $1, "0", $2}}' {output.fai} > {output.ref_bed} 2>> {log} &&
        picard CreateSequenceDictionary -Xmx{params.mem} R={input.ref} O={output.fdict} 2>> {log}
        """

