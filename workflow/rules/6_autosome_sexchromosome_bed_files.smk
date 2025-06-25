##########################################################################
### 6 Generate BED files with sex-chromosomal and autosomal scaffolds (e.g. for mlRho or other downstream analyses)

# Code collecting output files from this part of the pipeline
autos_sexchr_bed_outputs=[]

if len(sexchromosomeList) > 0:
    autos_sexchr_bed_outputs.append(expand("results/" + REF_NAME + ".repma.{chr}.bed", chr=["autos", "sexchr", "genome"],))
    if config["CpG_from_vcf"] == True:
        autos_sexchr_bed_outputs.append(expand("results/" + REF_NAME + ".noCpG_vcf.repma.{chr}.bed",
            chr=["autos", "sexchr", "genome"],))
    elif config["CpG_from_reference"] == True:
        autos_sexchr_bed_outputs.append(expand("results/" + REF_NAME + ".noCpG_ref.repma.{chr}.bed",
            chr=["autos", "sexchr", "genome"],))
    elif config["CpG_from_vcf_and_reference"] == True:
        autos_sexchr_bed_outputs.append(expand("results/" + REF_NAME + ".noCpG_vcfref.repma.{chr}.bed",
            chr=["autos", "sexchr", "genome"],))


# snakemake rules
localrules: rename_genome_bed, rename_noCpG_genome_beds

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
        genomefile=rules.genome_file.output.genomefile,
    output:
        autosome_bed=REF_DIR + "/" + REF_NAME + ".autos.bed",
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + "_make_autosomes_bed.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools subtract -a {input.ref_bed} -b {input.sexchr_bed} \
        -sorted -g {input.genomefile} > {output.autosome_bed} 2> {log}
        """


rule intersect_sexchr_repma_beds:
    input:
        no_rep_bed=rules.make_no_repeats_bed.output.no_rep_bed,
        sexchr_bed=rules.make_sexchr_bed.output,
        genomefile=rules.genome_file.output.genomefile,
    output:
        repma_sex_chr="results/" + REF_NAME + ".repma.sexchr.bed",
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + "_intersect_sexchr_repma_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools intersect -a {input.no_rep_bed} -b {input.sexchr_bed} \
        -sorted -g {input.genomefile} > {output.repma_sex_chr} 2> {log}
        """


rule intersect_autos_repma_beds:
    input:
        no_rep_bed=rules.make_no_repeats_bed.output.no_rep_bed,
        autosome_bed=rules.make_autosomes_bed.output,
        genomefile=rules.genome_file.output.genomefile,
    output:
        repma_autos="results/" + REF_NAME + ".repma.autos.bed",
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + "_intersect_autos_repma_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools intersect -a {input.no_rep_bed} -b {input.autosome_bed} \
        -sorted -g {input.genomefile} > {output.repma_autos} 2> {log}
        """


rule intersect_sexchr_noCpG_repma_beds:
    input:
        no_CpG_repma_bed="results/" + REF_NAME + ".no{CpG_method}.repma.bed",
        sexchr_bed=rules.make_sexchr_bed.output,
        genomefile=rules.genome_file.output.genomefile,
    output:
        no_CpG_repma_sexchr="results/" + REF_NAME + ".no{CpG_method}.repma.sexchr.bed",
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + ".no{CpG_method}_intersect_sexchr_noCpG_repma_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools intersect -a {input.no_CpG_repma_bed} -b {input.sexchr_bed} \
        -sorted -g {input.genomefile} > {output.no_CpG_repma_sexchr} 2> {log}
        """


rule intersect_autos_noCpG_repma_beds:
    input:
        no_CpG_repma_bed="results/" + REF_NAME + ".no{CpG_method}.repma.bed",
        autosome_bed=rules.make_autosomes_bed.output,
        genomefile=rules.genome_file.output.genomefile,
    output:
        no_CpG_repma_autos="results/" + REF_NAME + ".no{CpG_method}.repma.autos.bed",
    threads: 2
    message:
        "the input files are: {input}"
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + ".no{CpG_method}_intersect_autos_noCpG_repma_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools intersect -a {input.no_CpG_repma_bed} -b {input.autosome_bed} \
        -sorted -g {input.genomefile} > {output.no_CpG_repma_autos} 2> {log}
        """


rule rename_genome_bed:
    input:
        REF_DIR + "/" + REF_NAME + ".repma.bed",
    output:
        "results/" + REF_NAME + ".repma.genome.bed",
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + ".rename_genome_bed.log",
    shell:
        """
        cp {input} {output} 2> {log}
        """

rule rename_noCpG_genome_beds:
    input:
        "results/" + REF_NAME + ".no{CpG_method}.repma.bed",
    output:
        "results/" + REF_NAME + ".no{CpG_method}.repma.genome.bed",
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + ".no{CpG_method}_rename_noCpG_genome_beds.log",
    shell:
        """
        cp {input} {output} 2> {log}
        """