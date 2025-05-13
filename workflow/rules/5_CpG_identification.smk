##########################################################################
### 5. Identification of CpG sites in VCF files and the reference genome, merging of CpG bed file with repeats bed file for downstream filtering

# Code collecting output files from this part of the pipeline
if config["CpG_from_vcf"] == True:
    all_outputs.append("results/" + REF_NAME + ".CpG_vcf.bed")
    all_outputs.append("results/" + REF_NAME + ".noCpG_vcf.bed")
    all_outputs.append("results/" + REF_NAME + ".noCpG_vcf.repma.bed")
    all_outputs.append("results/" + REF_NAME + ".CpG_vcf.repeats.bed")

elif config["CpG_from_reference"] == True:
    all_outputs.append("results/" + REF_NAME + ".CpG_ref.bed")
    all_outputs.append("results/" + REF_NAME + ".noCpG_ref.bed")
    all_outputs.append("results/" + REF_NAME + ".noCpG_ref.repma.bed")
    all_outputs.append("results/" + REF_NAME + ".CpG_ref.repeats.bed")

elif config["CpG_from_vcf_and_reference"] == True:
    all_outputs.append("results/" + REF_NAME + ".CpG_vcfref.bed")
    all_outputs.append("results/" + REF_NAME + ".noCpG_vcfref.bed")
    all_outputs.append("results/" + REF_NAME + ".noCpG_vcfref.repma.bed")
    all_outputs.append("results/" + REF_NAME + ".CpG_vcfref.repeats.bed")


# Functions used by rules of this part of the pipeline
def CpG_genotype_bed_files_to_merge(wildcards):
    """Collect bed files with CpG sites found in individual samples as input for merge_CpG_genotype_beds"""
    if config["CpG_from_vcf"] == True:
        hist_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.CpG.bed",
            sample=HIST_CpG_SAMPLES,)
    return hist_CpG

def all_CpG_bed_files_to_merge(wildcards):
    """Collect bed files with CpG sites found in individual samples and in the reference genome as input for merge_all_CpG_beds"""
    if config["CpG_from_vcf_and_reference"] == True:
        hist_CpG = expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.CpG.bed",
            sample=HIST_CpG_SAMPLES,)
        mod_CpG = expand("results/modern/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.CpG.bed",
            sample=MODERN_CpG_SAMPLES,)
        ref = expand("results/" + REF_NAME + ".CpG_ref.bed")
    return hist_CpG + mod_CpG + ref


# snakemake rules
rule sorted_bcf2vcf_CpG_id:
    """Convert bcf format to vcf.gz for removal of sites"""
    input:
        bcf="results/{dataset}/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.bcf",
    output:
        vcf=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.CpG_id.vcf.gz"),
    log:
        "results/logs/5_CpG_identification/{dataset}/" + REF_NAME + "/{sample}_sorted_bcf2vcf_CpG_id.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools convert -O z -o {output.vcf} {input.bcf} 2> {log}
        """


rule make_CpG_genotype_bed:
    """Make a bed file of CpG sites for each single individual vcf file"""
    """CpG sites are only identified in variant calls, CpG only found in the reference genome are ignored"""
    input:
        vcf=rules.sorted_bcf2vcf_CpG_id.output.vcf,
    output:
        bed=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.CpG.bed"),
    log:
        "results/logs/5_CpG_identification/{dataset}/" + REF_NAME + "/{sample}_make_CpG_genotype_bed.log",
    shell:
        """
        python workflow/scripts/find_CpG_genotypes.py {input.vcf} {output.bed} 2> {log}
        """


rule make_CpG_reference_bed:
    """Make a bed file of CpG sites found in the reference genome"""
    input:
        ref=rules.ref_upper.output,
    output:
        bed="results/" + REF_NAME + ".CpG_ref.bed",
    params:
        refdirbed=REF_DIR + "/" + REF_NAME + ".CpG_ref",
    log:
        "results/logs/5_CpG_identification/" + REF_NAME + "_make_CpG_reference_bed.log",
    shell:
        """
        python workflow/scripts/find_CpG_ref_sites.py {input.ref} {params.refdirbed} 2> {log} &&
        cp {params.refdirbed}.bed {output.bed} 2>> {log}
        """


rule merge_CpG_genotype_beds:
    input:
        CpG_genotype_bed_files_to_merge,
    output:
        tmp=temp("results/" + REF_NAME + ".concatenated.CpG_vcf.bed"),
        merged=temp("results/" + REF_NAME + ".merged.CpG_vcf.bed"),
    message:
        "the input files are: {input}"
    group:
        "CpG_genotype_bed_formatting_group"
    log:
        "results/logs/5_CpG_identification/" + REF_NAME + "_merge_CpG_genotype_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        files=`echo {input} | awk '{{print NF}}'`
        if [ $files -gt 1 ] # check if there are at least 2 files for merging. If there is only one file, copy the sorted bam file.
        then
            cat {input} | sort -k1,1 -k2,2n > {output.tmp} 2> {log} &&
            bedtools merge -i {output.tmp} > {output.merged} 2>> {log}
        else
            touch {output.tmp} && cp {input} {output.merged} 2> {log} &&
            echo "Only one file present for merging. Copying the input bed file." >> {log}
        fi
        """


rule sort_CpG_genotype_beds:
    input:
        merged_bed=rules.merge_CpG_genotype_beds.output.merged,
        genomefile=rules.genome_file.output.genomefile,
    output:
        sorted_bed="results/" + REF_NAME + ".CpG_vcf.bed",
    group:
        "CpG_genotype_bed_formatting_group"
    log:
        "results/logs/5_CpG_identification/" + REF_NAME + "_sort_CpG_genotype_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools sort -g {input.genomefile} -i {input.merged_bed} > {output.sorted_bed} 2> {log}
        """


rule merge_all_CpG_beds:
    input:
        all_CpG_bed_files_to_merge,
    output:
        tmp=temp("results/" + REF_NAME + ".concatenated.CpG_vcfref.bed"),
        merged=temp("results/" + REF_NAME + ".merged.CpG_vcfref.bed"),
    message:
        "the input files are: {input}"
    group:
        "all_CpG_bed_formatting_group"
    log:
        "results/logs/5_CpG_identification/" + REF_NAME + "_merge_all_CpG_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        files=`echo {input} | awk '{{print NF}}'`
        if [ $files -gt 1 ] # check if there are at least 2 files for merging. If there is only one file, copy the sorted bam file.
        then
            cat {input} | sort -k1,1 -k2,2n > {output.tmp} 2> {log} &&
            bedtools merge -i {output.tmp} > {output.merged} 2>> {log}
        else
            touch {output.tmp} && cp {input} {output.merged} 2> {log} &&
            echo "Only one file present for merging. Copying the input bed file." >> {log}
        fi
        """


rule sort_all_CpG_beds:
    input:
        merged_bed=rules.merge_all_CpG_beds.output.merged,
        genomefile=rules.genome_file.output.genomefile,
    output:
        sorted_bed="results/" + REF_NAME + ".CpG_vcfref.bed",
    group:
        "all_CpG_bed_formatting_group"
    log:
        "results/logs/5_CpG_identification/" + REF_NAME + "_sort_all_CpG_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools sort -g {input.genomefile} -i {input.merged_bed} > {output.sorted_bed} 2> {log}
        """


rule make_noCpG_bed:
    input:
        ref_bed=rules.make_reference_bed.output,
        CpG_bed="results/" + REF_NAME + ".{CpG_method}.bed",
    output:
        no_CpG_bed="results/" + REF_NAME + ".no{CpG_method}.bed",
    threads: 2
    log:
        "results/logs/5_CpG_identification/" + REF_NAME + ".no{CpG_method}_make_no_CpG_bed.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools subtract -a {input.ref_bed} -b {input.CpG_bed} > {output.no_CpG_bed} 2> {log}
        """


rule merge_CpG_repeats_beds:
    input:
        CpG_bed="results/" + REF_NAME + ".{CpG_method}.bed",
        sorted_rep_bed=rules.sort_repeats_bed.output,
    output:
        tmp=temp("results/" + REF_NAME + ".concatenated.{CpG_method}.repeats.bed"),
        merged=temp("results/" + REF_NAME + ".merged.{CpG_method}.repeats.bed"),
    message:
        "the input files are: {input}"
    group:
        "CpG_repeats_bed_formatting_group"
    log:
        "results/logs/5_CpG_identification/" + REF_NAME + ".{CpG_method}_merge_CpG_repeats_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        cat {input[0]} {input[1]} | sort -k1,1 -k2,2n > {output.tmp} 2> {log} &&
        bedtools merge -i {output.tmp} > {output.merged} 2>> {log}
        """


rule sort_CpG_repeats_beds:
    input:
        merged_bed=rules.merge_CpG_repeats_beds.output.merged,
        genomefile=rules.genome_file.output.genomefile,
    output:
        sorted_bed="results/" + REF_NAME + ".{CpG_method}.repeats.bed",
    group:
        "CpG_repeats_bed_formatting_group"
    log:
        "results/logs/5_CpG_identification/" + REF_NAME + ".{CpG_method}_sort_CpG_repeats_beds.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools sort -g {input.genomefile} -i {input.merged_bed} > {output.sorted_bed} 2> {log}
        """


rule make_noCpG_repma_bed:
    input:
        ref_bed=rules.make_reference_bed.output,
        merged_bed=rules.merge_CpG_repeats_beds.output.merged,
    output:
        no_CpG_repma_bed="results/" + REF_NAME + ".no{CpG_method}.repma.bed",
    threads: 2
    log:
        "results/logs/5_CpG_identification/" + REF_NAME + ".no{CpG_method}_make_noCpG_repma_bed.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools subtract -a {input.ref_bed} -b {input.merged_bed} > {output.no_CpG_repma_bed} 2> {log}
        """
