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
rule make_sexchr_autosomes_bed:
    """Generate a bed file of autosomal scaffolds/chromosomes and
    a bed file of sex chromosome-linked scaffolds/chromosomes"""
    input:
        fai=config["ref_path"] + ".fai",
        ref_bed=rules.index_reference.output.ref_bed,
        genomefile=rules.index_reference.output.fai,
        no_rep_bed=rules.make_no_repeats_bed.output.no_rep_bed,
    output:
        sexchr_bed=REF_DIR + "/" + REF_NAME + ".sexchr.bed",
        autosome_bed=REF_DIR + "/" + REF_NAME + ".autos.bed",
        repma_sex_chr="results/" + REF_NAME + ".repma.sexchr.bed",
        repma_autos="results/" + REF_NAME + ".repma.autos.bed",
    params:
        sexchr=sexchromosomeList,  # list of sex chromosomal scaffold names created from file provided by user
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + "_make_sexchr_autosomes_bed.log",
    container:
        bedtools_htslib_container
    shell:
        """
        # sex chromosome bed file
        for i in {params.sexchr}; do grep -w $i {input.fai} | \
        awk 'BEGIN {{FS="\t"}}; {{print $1 FS "0" FS $2}}' >> \
        {output.sexchr_bed} 2>> {log}; done &&
        
        # autosomal bed file
        bedtools subtract -a {input.ref_bed} -b {output.sexchr_bed} \
        -sorted -g {input.genomefile} > {output.autosome_bed} 2>> {log} &&

        # intersect sexchr repma beds
        bedtools intersect -a {input.no_rep_bed} -b {output.sexchr_bed} \
        -sorted -g {input.genomefile} > {output.repma_sex_chr} 2>> {log} &&

        # intersect autos repma beds
        bedtools intersect -a {input.no_rep_bed} -b {output.autosome_bed} \
        -sorted -g {input.genomefile} > {output.repma_autos} 2> {log}
        """


rule intersect_sexchr_autos_noCpG_repma_beds:
    input:
        no_CpG_repma_bed="results/" + REF_NAME + ".no{CpG_method}.repma.bed",
        sexchr_bed=rules.make_sexchr_autosomes_bed.output.sexchr_bed,
        autosome_bed=rules.make_sexchr_autosomes_bed.output.autosome_bed,
        genomefile=rules.index_reference.output.fai,
    output:
        no_CpG_repma_sexchr="results/" + REF_NAME + ".no{CpG_method}.repma.sexchr.bed",
        no_CpG_repma_autos="results/" + REF_NAME + ".no{CpG_method}.repma.autos.bed",
    threads: 2
    log:
        "results/logs/6_autosome_sexchromosome_bed_files/" + REF_NAME + ".no{CpG_method}_intersect_sexchr_autos_noCpG_repma_beds.log",
    container:
        bedtools_htslib_container
    shell:
        """
        # no CpG repma sex chromosome bed
        bedtools intersect -a {input.no_CpG_repma_bed} -b {input.sexchr_bed} \
        -sorted -g {input.genomefile} > {output.no_CpG_repma_sexchr} 2> {log} &&

        # no CpG repma autosomes bed
        bedtools intersect -a {input.no_CpG_repma_bed} -b {input.autosome_bed} \
        -sorted -g {input.genomefile} > {output.no_CpG_repma_autos} 2>> {log}
        """
