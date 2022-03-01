##########################################################################
### 6.3 Generate BED files with sex-chromosomal and autosomal scaffolds (e.g. for mlRho or other downstream analyses)

# Code collecting output files from this part of the pipeline
if len(sexchromosomeList) > 0:
    all_outputs.append(expand("results/" + REF_NAME + ".repma.{chr}.bed", chr=["autos", "sexchr"],))
    if config["CpG_from_vcf"] == True:
        all_outputs.append(expand("results/" + REF_NAME + ".noCpG_vcf.repma.{chr}.bed",
            chr=["autos", "sexchr"],))
    elif config["CpG_from_reference"] == True:
        all_outputs.append(expand("results/" + REF_NAME + ".noCpG_ref.repma.{chr}.bed",
            chr=["autos", "sexchr"],))
    elif config["CpG_from_vcf_and_reference"] == True:
        all_outputs.append(expand("results/" + REF_NAME + ".noCpG_vcfref.repma.{chr}.bed",
            chr=["autos", "sexchr"],))


# snakemake rules
rule make_sexchr_bed:
    """Generate a bed file of sex chromosome-linked contigs/scaffolds"""
    input:
        fai=config["ref_path"] + ".fai",
    output:
        sexchr_bed=REF_DIR + "/" + REF_NAME + ".sexchr.bed",
    params:
        sexchr=sexchromosomeList,  # list of sex chromosomal scaffold names created from file provided by user
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + "_make_sexchr_bed.log",
    shell:
        """
        for i in {params.sexchr}; do grep -w $i {input.fai} | awk 'BEGIN {{FS="\t"}}; {{print $1 FS "0" FS $2}}' >> {output.sexchr_bed} 2>> {log}; done
        """


rule make_autosomes_bed:
    """Generate a bed file of autosomal contigs/scaffolds"""
    input:
        ref_bed=rules.make_reference_bed.output,
        sexchr_bed=rules.make_sexchr_bed.output,
    output:
        autosome_bed=REF_DIR + "/" + REF_NAME + ".autos.bed",
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + "_make_autosomes_bed.log",
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    shell:
        """
        bedtools subtract -a {input.ref_bed} -b {input.sexchr_bed} > {output.autosome_bed} 2> {log}
        """


rule intersect_sexchr_repma_beds:
    input:
        no_rep_bed_dir=rules.make_no_repeats_bed.output.no_rep_bed_dir,
        sexchr_bed=rules.make_sexchr_bed.output,
    output:
        repma_sex_chr="results/" + REF_NAME + ".repma.sexchr.bed",
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + "_intersect_sexchr_repma_beds.log",
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    shell:
        """
        bedtools intersect -a {input.no_rep_bed_dir} -b {input.sexchr_bed} > {output.repma_sex_chr} 2> {log}
        """


rule intersect_autos_repma_beds:
    input:
        no_rep_bed_dir=rules.make_no_repeats_bed.output.no_rep_bed_dir,
        autosome_bed=rules.make_autosomes_bed.output,
    output:
        repma_autos="results/" + REF_NAME + ".repma.autos.bed",
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + "_intersect_autos_repma_beds.log",
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    shell:
        """
        bedtools intersect -a {input.no_rep_bed_dir} -b {input.autosome_bed} > {output.repma_autos} 2> {log}
        """


rule intersect_sexchr_noCpG_repma_beds:
    input:
        no_CpG_repma_bed="results/" + REF_NAME + ".no{CpG_method}.repma.bed",
        sexchr_bed=rules.make_sexchr_bed.output,
    output:
        no_CpG_repma_sexchr="results/" + REF_NAME + ".no{CpG_method}.repma.sexchr.bed",
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + ".no{CpG_method}_intersect_sexchr_noCpG_repma_beds.log",
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    shell:
        """
        bedtools intersect -a {input.no_CpG_repma_bed} -b {input.sexchr_bed} > {output.no_CpG_repma_sexchr} 2> {log}
        """


rule intersect_autos_noCpG_repma_beds:
    input:
        no_CpG_repma_bed="results/" + REF_NAME + ".no{CpG_method}.repma.bed",
        autosome_bed=rules.make_autosomes_bed.output,
    output:
        no_CpG_repma_autos="results/" + REF_NAME + ".no{CpG_method}.repma.autos.bed",
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + ".no{CpG_method}_intersect_autos_noCpG_repma_beds.log",
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    shell:
        """
        bedtools intersect -a {input.no_CpG_repma_bed} -b {input.autosome_bed} > {output.no_CpG_repma_autos} 2> {log}
        """
