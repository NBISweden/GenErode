##########################################################################
### 12. Annotate SNPs in merged BCF files using snpEff
# Authors: Verena Kutschera, Marcin Kierczak

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
    all_outputs.append(expand("results/all/snpEff/" + REF_NAME + ".all.fmissing{fmiss}.{chr}.snpEff_variant_impact_plot.pdf", 
        fmiss=config["f_missing"],
        chr=CHR,))
    all_outputs.append("results/historical/snpEff/" + REF_NAME + "/multiqc/multiqc_report.html")
    all_outputs.append("results/modern/snpEff/" + REF_NAME + "/multiqc/multiqc_report.html")

elif os.path.exists(config["historical_samples"]):
    all_outputs.append(expand("results/historical/snpEff/" + REF_NAME + ".historical.fmissing{fmiss}.{chr}.snpEff_variant_impact_plot.pdf",
        fmiss=config["f_missing"],
        chr=CHR,))
    all_outputs.append("results/historical/snpEff/" + REF_NAME + "/multiqc/multiqc_report.html")

elif os.path.exists(config["modern_samples"]):
    all_outputs.append(expand("results/modern/snpEff/" + REF_NAME + ".modern.fmissing{fmiss}.{chr}.snpEff_variant_impact_plot.pdf",
        fmiss=config["f_missing"],
        chr=CHR,))
    all_outputs.append("results/modern/snpEff/" + REF_NAME + "/multiqc/multiqc_report.html")


# Functions used by rules of this part of the pipeline
def historical_snpEff_multiqc_inputs(wildcards):
    """Input for historical_snpEff_multiqc"""
    outlist = []
    hist_not_CpG = expand("results/historical/snpEff/" + REF_NAME + "/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.csv",
        sample=HIST_NOT_CpG_SAMPLES,
        fmiss=config["f_missing"],
        chr=CHR,)
    outlist += (hist_not_CpG)
    if config["CpG_from_vcf"] == True:
        hist_CpG = expand("results/historical/snpEff/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.csv",
            sample=HIST_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
        outlist += (hist_CpG)
    elif config["CpG_from_reference"] == True:
        hist_CpG = expand("results/historical/snpEff/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.csv",
            sample=HIST_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
        outlist += (hist_CpG)
    elif config["CpG_from_vcf_and_reference"] == True:
        hist_CpG = expand("results/historical/snpEff/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.csv",
            sample=HIST_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
        outlist += (hist_CpG)
    return outlist

def modern_snpEff_multiqc_inputs(wildcards):
    """Input for modern_snpEff_multiqc"""
    outlist = []
    mod_not_CpG = expand("results/modern/snpEff/" + REF_NAME + "/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.csv",
        sample=MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,
        fmiss=config["f_missing"],
        chr=CHR,)
    outlist += (mod_not_CpG)
    if config["CpG_from_vcf"] == True:
        mod_CpG = expand("results/modern/snpEff/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.csv",
            sample=MODERN_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
        outlist += (mod_CpG)
    elif config["CpG_from_reference"] == True:
        mod_CpG = expand("results/modern/snpEff/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.csv",
            sample=MODERN_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
        outlist += (mod_CpG)
    elif config["CpG_from_vcf_and_reference"] == True:
        mod_CpG = expand("results/modern/snpEff/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.csv",
            sample=MODERN_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
        outlist += (mod_CpG)
    return outlist


def all_snpEff_outputs(wildcards):
    """Collect output files for pipeline report"""
    if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
        return historical_snpEff_multiqc_inputs(wildcards) + modern_snpEff_multiqc_inputs(wildcards)
    elif os.path.exists(config["modern_samples"]):
        return modern_snpEff_multiqc_inputs(wildcards)
    elif os.path.exists(config["historical_samples"]):
        return historical_snpEff_multiqc_inputs(wildcards)


# snakemake rules
localrules:
    prepare_db_build,
    snpEff_variant_impact_plot,


rule prepare_db_build:
    """Prepare file structure and files to build database for snpEff annotation"""
    input:
        gtf=config["gtf_path"],
        ref=config["ref_path"],
    output:
        gtf=GTF_DIR + "/snpEff/data/" + REF_NAME + "/genes.gtf",
        ref=GTF_DIR + "/snpEff/data/genomes/" + REF_NAME + ".fa",
    params:
        abs_gtf_in=lambda wildcards, input: os.path.abspath(input.gtf),
        abs_gtf_out=lambda wildcards, output: os.path.abspath(output.gtf),
        abs_ref_in=lambda wildcards, input: os.path.abspath(input.ref),
        abs_ref_out=lambda wildcards, output: os.path.abspath(output.ref),
    log:
        "results/logs/12_snpEff/" + REF_NAME + "_prepare_db_build.log",
    shell:
        """
        ln -s {params.abs_gtf_in} {params.abs_gtf_out} 2> {log} &&
        ln -s {params.abs_ref_in} {params.abs_ref_out} 2>> {log}
        """


rule update_snpEff_config:
    """Edit snpEff config file to build database"""
    input:
        gtf=rules.prepare_db_build.output.gtf,
        ref=rules.prepare_db_build.output.ref,
    output:
        config=GTF_DIR + "/snpEff/data/" + REF_NAME + "/snpEff.config",
    params:
        config="/usr/local/share/snpeff-4.3.1t-3/snpEff.config",
        ref_name=REF_NAME,
        abs_gtf=lambda wildcards, input: os.path.abspath(input.gtf),
        abs_ref=lambda wildcards, input: os.path.abspath(input.ref),
        abs_config=lambda wildcards, output: os.path.abspath(output.config),
    log:
        "results/logs/12_snpEff/" + REF_NAME + "_update_snpEff_config.log",
    singularity:
        snpeff_container
    shell:
        """
        cp {params.config} {params.abs_config} 2> {log} &&
        echo '#{params.ref_name} genome, version {params.ref_name}' >> {params.abs_config} 2>> {log} &&
        echo '{params.ref_name}.genome : {params.ref_name}' >> {params.abs_config} 2>> {log}
        """


rule build_snpEff_db:
    """Build snpEff database for the reference genome and annotation file"""
    input:
        gtf=rules.prepare_db_build.output.gtf,
        ref=rules.prepare_db_build.output.ref,
        config=rules.update_snpEff_config.output.config,
    output:
        db=GTF_DIR + "/snpEff/data/" + REF_NAME + "/snpEffectPredictor.bin",
    threads: 1
    resources:
        mem_mb=8000,
    params:
        ref_name=REF_NAME,
        abs_gtf=lambda wildcards, input: os.path.abspath(input.gtf),
        abs_ref=lambda wildcards, input: os.path.abspath(input.ref),
        abs_config=lambda wildcards, input: os.path.abspath(input.config),
        abs_db_dir=os.path.abspath(GTF_DIR + "/snpEff/"),
        abs_data_dir=os.path.abspath(GTF_DIR + "/snpEff/data/"),
    log:
        os.path.abspath("results/logs/12_snpEff/" + REF_NAME + "_build_snpEff_db.log"),
    singularity:
        snpeff_container
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        cd {params.abs_db_dir}
        java -jar -Xmx${{mem}}g /usr/local/share/snpeff-4.3.1t-3/snpEff.jar build -gtf22 -c {params.abs_config} \
        -dataDir {params.abs_data_dir} -treatAllAsProteinCoding -v {params.ref_name} 2> {log}
        """


rule annotate_vcf:
    """Annotate the VCF files of each individual"""
    input:
        vcf=rules.filter_biallelic_missing_vcf.output.filtered,
        db=rules.build_snpEff_db.output.db,
        config=rules.update_snpEff_config.output.config,
    output:
        ann="results/{dataset}/snpEff/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.ann.vcf",
        csv="results/{dataset}/snpEff/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.csv",
        html="results/{dataset}/snpEff/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}_stats.html",
    threads: 1
    resources:
        mem_mb=8000,
    params:
        ref_name=REF_NAME,
        abs_config=lambda wildcards, input: os.path.abspath(input.config),
        abs_data_dir=os.path.abspath(GTF_DIR + "/snpEff/data/"),
    log:
        "results/logs/12_snpEff/{dataset}/" + REF_NAME + "/{sample}.{filtered}_fmissing{fmiss}.{chr}_annotate_vcf.log",
    singularity:
        snpeff_container
    shell:
        """
        mem=$((({resources.mem_mb} - 2000)/1000))
        java -jar -Xmx${{mem}}g /usr/local/share/snpeff-4.3.1t-3/snpEff.jar  -c {params.abs_config} -dataDir {params.abs_data_dir} -s {output.html} -csvStats {output.csv} \
        -treatAllAsProteinCoding -v -d -lof {params.ref_name} {input.vcf} > {output.ann} 2> {log}
        """


rule historical_snpEff_multiqc:
    """Run multiQC on snpEff results to get a summary for all historical samples"""
    input:
        historical_snpEff_multiqc_inputs,
    output:
        stats=report("results/historical/snpEff/" + REF_NAME + "/multiqc/multiqc_report.html",
            caption="../report/historical_snpEff_multiqc.rst",
            category="snpEff",),
    params:
        indir="results/historical/snpEff/" + REF_NAME + "/",
        outdir="results/historical/snpEff/" + REF_NAME + "/multiqc/",
    log:
        "results/logs/12_snpEff/historical/" + REF_NAME + "/historical_snpEff_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_snpEff_multiqc:
    """Run multiQC on snpEff results to get a summary for all modern samples"""
    input:
        modern_snpEff_multiqc_inputs,
    output:
        stats=report("results/modern/snpEff/" + REF_NAME + "/multiqc/multiqc_report.html",
            caption="../report/modern_snpEff_multiqc.rst",
            category="snpEff",),
    params:
        indir="results/modern/snpEff/" + REF_NAME + "/",
        outdir="results/modern/snpEff/" + REF_NAME + "/multiqc/",
    log:
        "results/logs/12_snpEff/modern/" + REF_NAME + "/modern_snpEff_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule snpEff_variant_impact_table:
    """Extract numbers of SNPs with high, moderate, low, and modifier impact, for modern and historical samples."""
    input:
        all_snpEff_outputs,
    output:
        table=report("results/{dataset}/snpEff/" + REF_NAME + ".{dataset}.fmissing{fmiss}.{chr}.snpEff_variant_impact_table.txt",
            caption="../report/snpEff_variant_impact_table.rst",
            category="snpEff",),
    log:
        "results/logs/12_snpEff/{dataset}/" + REF_NAME + ".{dataset}.fmissing{fmiss}.{chr}.snpEff_variant_impact_table.log",
    script:
        "../scripts/snpEff_variant_impact_table.py"


rule snpEff_variant_impact_plot:
    """Plot numbers of SNPs with high, moderate, low, and modifier impact, for modern and historical samples."""
    input:
        table= "results/{dataset}/snpEff/" + REF_NAME + ".{dataset}.fmissing{fmiss}.{chr}.snpEff_variant_impact_table.txt",
    output:
        plot=report("results/{dataset}/snpEff/" + REF_NAME + ".{dataset}.fmissing{fmiss}.{chr}.snpEff_variant_impact_plot.pdf",
            caption="../report/snpEff_variant_impact_plot.rst",
            category="snpEff",),
    log:
        "results/logs/12_snpEff/{dataset}/" + REF_NAME + ".{dataset}.fmissing{fmiss}.{chr}.snpEff_variant_impact_plot.log",
    script:
        "../scripts/snpEff_variant_impact_plot.py"
