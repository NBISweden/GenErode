##########################################################################
### 8.1 Removal of CpG-prone sites from VCF files

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    if len(HIST_CpG_SAMPLES) > 0:
        all_outputs.append("results/historical/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/multiqc/multiqc_report.html")

if os.path.exists(config["modern_samples"]):
    if len(MOD_CpG_SAMPLES) > 0:
        all_outputs.append("results/modern/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/multiqc/multiqc_report.html")


# Functions used by rules of this part of the pipeline
def historical_CpG_filtered_multiqc_inputs(wildcards):
    """Input for historical_CpG_filtered_multiqc"""
    if config["CpG_from_vcf"] == True:
        return expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/{sample}.Q30.q30.sorted.noCpG_vcf.bcf.stats.txt",
            sample=HIST_CpG_SAMPLES,)
    elif config["CpG_from_reference"] == True:
        return expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/{sample}.Q30.q30.sorted.noCpG_ref.bcf.stats.txt",
            sample=HIST_CpG_SAMPLES,)
    elif config["CpG_from_vcf_and_reference"] == True:
        return expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/{sample}.Q30.q30.sorted.noCpG_vcfref.bcf.stats.txt",
            sample=HIST_CpG_SAMPLES,)

def modern_CpG_filtered_multiqc_inputs(wildcards):
    """Input for modern_CpG_filtered_multiqc"""
    if config["CpG_from_vcf"] == True:
        return expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/{sample}.Q30.q30.sorted.noCpG_vcf.bcf.stats.txt",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_reference"] == True:
        return expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/{sample}.Q30.q30.sorted.noCpG_ref.bcf.stats.txt",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_vcf_and_reference"] == True:
        return expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/{sample}.Q30.q30.sorted.noCpG_vcfref.bcf.stats.txt",
            sample=MOD_CpG_SAMPLES,)


# snakemake rules
rule sorted_bcf2vcf_CpG_removal:
    """Convert bcf format to vcf.gz for removal of sites"""
    input:
        bcf=rules.sort_vcfs.output.sort,
    output:
        vcf=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.CpG_rm.vcf.gz"),
    log:
        "results/logs/8.1_vcf_CpG_filtering/{dataset}/" + REF_NAME + "/{sample}_sorted_bcf2vcf_CpG_removal.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools convert -O z -o {output.vcf} {input.bcf} 2> {log}
        """


rule remove_CpG_vcf:
    """Remove CpG-prone sites from vcf file"""
    input:
        vcf=rules.sorted_bcf2vcf_CpG_removal.output.vcf,
        bed=rules.make_noCpG_bed.output.no_CpG_bed,
        genomefile=rules.genome_file.output.genomefile,
    output:
        filtered=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.no{CpG_method}.vcf.gz"),
    threads: 6
    log:
        "results/logs/8.1_vcf_CpG_filtering/{dataset}/" + REF_NAME + "/{sample}.no{CpG_method}_remove_CpG_vcf.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools intersect -a {input.vcf} -b {input.bed} -header -sorted -g {input.genomefile} | bgzip -c > {output.filtered} 2> {log}
        """


rule CpG_vcf2bcf:
    """Convert vcf back to bcf"""
    input:
        filtered=rules.remove_CpG_vcf.output.filtered,
    output:
        bcf="results/{dataset}/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.no{CpG_method}.bcf",
    threads: 2
    log:
        "results/logs/8.1_vcf_CpG_filtering/{dataset}/" + REF_NAME + "/{sample}.no{CpG_method}_CpG_vcf2bcf.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools convert -O b -o {output.bcf} {input.filtered} 2> {log}
        """


rule index_CpG_bcf:
    """Index the bcf file"""
    input:
        bcf=rules.CpG_vcf2bcf.output.bcf,
    output:
        index="results/{dataset}/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.no{CpG_method}.bcf.csi",
    group:
        "CpG_bcf_group"
    log:
        "results/logs/8.1_vcf_CpG_filtering/{dataset}/" + REF_NAME + "/{sample}.no{CpG_method}_index_CpG_bcf.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools index -o {output.index} {input.bcf} 2> {log}
        """


rule CpG_filtered_vcf_stats:
    """Obtain summary stats of vcf files filtered for CpG sites"""
    input:
        bcf=rules.CpG_vcf2bcf.output.bcf,
        index=rules.index_CpG_bcf.output.index,
    output:
        stats="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/{sample}.Q30.q30.sorted.no{CpG_method}.bcf.stats.txt",
    group:
        "CpG_bcf_group"
    log:
        "results/logs/8.1_vcf_CpG_filtering/{dataset}/" + REF_NAME + "/{sample}.no{CpG_method}_CpG_filtered_vcf_stats.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools stats {input.bcf} > {output.stats} 2> {log}
        """


rule historical_CpG_filtered_vcf_multiqc:
    """Collect all stats files from CpG filtered vcf files from historical samples"""
    input:
        historical_CpG_filtered_multiqc_inputs,
    output:
        stats="results/historical/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/multiqc/multiqc_report.html",
    params:
        indir="results/historical/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/",
        outdir="results/historical/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/multiqc",
    log:
        "results/logs/8.1_vcf_CpG_filtering/historical/" + REF_NAME + "/historical_CpG_filtered_vcf_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_CpG_filtered_vcf_multiqc:
    """Collect all stats files from CpG filtered vcf files from modern samples"""
    input:
        modern_CpG_filtered_multiqc_inputs,
    output:
        stats="results/modern/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/multiqc/multiqc_report.html",
    params:
        indir="results/modern/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/",
        outdir="results/modern/vcf/" + REF_NAME + "/stats/vcf_CpG_filtered/multiqc",
    log:
        "results/logs/8.1_vcf_CpG_filtering/modern/" + REF_NAME + "/modern_CpG_filtered_vcf_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """