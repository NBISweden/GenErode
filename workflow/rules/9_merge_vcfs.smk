##########################################################################
### 9. Merge VCF files and filter for biallelic sites, missingness and sex-chromosomal contigs/scaffolds

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
    all_outputs.append(expand("results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_merged_missing/multiqc/multiqc_report.html",
        dataset=["all", "historical", "modern"],))
elif os.path.exists(config["historical_samples"]):
    all_outputs.append("results/historical/vcf/" + REF_NAME + "/stats/vcf_merged_missing/multiqc/multiqc_report.html")
elif os.path.exists(config["modern_samples"]):
    all_outputs.append("results/modern/vcf/" + REF_NAME + "/stats/vcf_merged_missing/multiqc/multiqc_report.html")


# Functions used by rules of this part of the pipeline
def merge_all_inputs(wildcards):
    """Input for merge_all_vcfs"""
    outlist = []
    outlist += expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
        sample=HIST_NOT_CpG_SAMPLES,)
    outlist += expand("results/modern/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
        sample=MOD_NOT_CpG_SAMPLES,)
    if config["CpG_from_vcf"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_CpG_SAMPLES,)
        outlist += expand("results/modern/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_CpG_SAMPLES,)
        outlist += expand("results/modern/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_vcf_and_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=HIST_CpG_SAMPLES,)
        outlist += expand("results/modern/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
            sample=MOD_CpG_SAMPLES,)
    return outlist


def merge_all_index_inputs(wildcards):
    """Input for merge_all_vcfs"""
    outlist = []
    outlist += expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
        sample=HIST_NOT_CpG_SAMPLES,)
    outlist += expand("results/modern/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
        sample=MOD_NOT_CpG_SAMPLES,)
    if config["CpG_from_vcf"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_CpG_SAMPLES,)
        outlist += expand("results/modern/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_CpG_SAMPLES,)
        outlist += expand("results/modern/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=MOD_CpG_SAMPLES,)
    elif config["CpG_from_vcf_and_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=HIST_CpG_SAMPLES,)
        outlist += expand("results/modern/vcf/" + REF_NAME + "/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
            sample=MOD_CpG_SAMPLES,)
    return outlist

def missingness_filtered_vcf_multiqc_inputs(wildcards):
    """Input for missingness_filtered_vcf_multiqc"""
    if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
        return expand("results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_merged_missing/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            dataset=["all", "historical", "modern"],
            fmiss=config["f_missing"],
            chr=CHR,)
    elif os.path.exists(config["historical_samples"]):
        return expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_merged_missing/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            fmiss=config["f_missing"],
            chr=CHR,)
    elif os.path.exists(config["modern_samples"]):
        return expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_merged_missing/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            fmiss=config["f_missing"],
            chr=CHR,)


def historical_biallelic_missing_filtered_vcf_multiqc_inputs(wildcards):
    """Input for historical_biallelic_missing_filtered_vcf_multiqc_inputs"""
    outlist = []
    outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
        sample=HIST_NOT_CpG_SAMPLES,
        fmiss=config["f_missing"],
        chr=CHR,)
    if config["CpG_from_vcf"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            sample=HIST_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
    elif config["CpG_from_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            sample=HIST_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
    elif config["CpG_from_vcf_and_reference"] == True:
        outlist += expand("results/historical/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            sample=HIST_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
    return outlist

def modern_biallelic_missing_filtered_vcf_multiqc_inputs(wildcards):
    """Input for modern_biallelic_missing_filtered_vcf_multiqc_inputs"""
    outlist = []
    outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/{sample}.Q30.q30.sorted.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
        sample=MOD_NOT_CpG_SAMPLES,
        fmiss=config["f_missing"],
        chr=CHR,)
    if config["CpG_from_vcf"] == True:
        outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/{sample}.Q30.q30.sorted.noCpG_vcf.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            sample=MOD_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
    elif config["CpG_from_reference"] == True:
        outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/{sample}.Q30.q30.sorted.noCpG_ref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            sample=MOD_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
    elif config["CpG_from_vcf_and_reference"] == True:
        outlist += expand("results/modern/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/{sample}.Q30.q30.sorted.noCpG_vcfref.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
            sample=MOD_CpG_SAMPLES,
            fmiss=config["f_missing"],
            chr=CHR,)
    return outlist


# snakemake rules
rule merge_all_vcfs:
    """Merge all samples into one VCF file, containing only SNPs"""
    input:
        bcf=merge_all_inputs,
        index=merge_all_index_inputs,
    output:
        merged="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf",
    threads: 6
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + "_merge_all_vcfs.log",
    singularity:
        bcftools_container
    shell:
        """
        files=`echo {input.bcf} | awk '{{print NF}}'`
        if [ $files -gt 1 ] # check if there are at least 2 files for merging. If there is only one file, copy the bcf file.
        then
            bcftools merge -m snps -O b -o {output.merged} {input.bcf} 2> {log}
        else
            cp {input.bcf} {output.merged} && touch {output.merged} 2> {log}
            echo "Only one file present for merging. Copying the input bcf file." >> {log}
        fi
        """


rule index_merged_vcf:
    """Index vcf files"""
    input:
        bcf="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf",
    output:
        index="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf.csi",
    group:
        "merged_vcf_group"
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_index_merged_vcfs.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools index -o {output.index} {input.bcf} 2> {log}
        """


rule merged_vcf_stats:
    """Obtain summary stats of merged vcf file before removing sites with missing data"""
    input:
        merged="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf",
        index=rules.index_merged_vcf.output,
    output:
        stats="results/all/vcf/" + REF_NAME + "/stats/vcf_merged/" + REF_NAME + ".all.merged.snps.bcf.stats.txt",
    group:
        "merged_vcf_group"
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_merged_vcf_stats.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools stats {input.merged} > {output.stats} 2> {log}
        """


rule merged_vcf_multiqc:
    """Collect all stats files from merged vcf files filtered for missing data"""
    input:
        rules.merged_vcf_stats.output.stats,
    output:
        stats="results/all/vcf/" + REF_NAME + "/stats/vcf_merged/multiqc/multiqc_report.html",
    params:
        indir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged/",
        outdir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged/multiqc",
    group:
        "merged_vcf_group"
    log:
        "results/logs/9_merge_vcfs/all/" + REF_NAME + "/merged_vcf_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule filter_vcf_biallelic:
    input:
        bcf="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf",
        index="results/all/vcf/" + REF_NAME + ".all.merged.snps.bcf.csi",
        stats=rules.merged_vcf_stats.output.stats,
        multiqc=rules.merged_vcf_multiqc.output.stats,
    output:
        bcf=temp("results/all/vcf/" + REF_NAME + ".all.merged.biallelic.bcf"),
        index=temp("results/all/vcf/" + REF_NAME + ".all.merged.biallelic.bcf.csi"),
    threads: 2
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_filter_vcf_biallelic.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools view -m2 -M2 -v snps -Ob -o {output.bcf} {input.bcf} 2> {log} &&
        bcftools index -f {output.bcf} 2>> {log}
        """


rule biallelic_filtered_vcf_stats:
    """Obtain summary stats of merged vcf file"""
    input:
        bcf=rules.filter_vcf_biallelic.output.bcf,
        index=rules.filter_vcf_biallelic.output.index,
    output:
        stats="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_biallelic/" + REF_NAME + ".all.merged.biallelic.vcf.stats.txt",
    group:
        "biallelic_filtered_vcf_group"
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_biallelic_filtered_vcf_stats.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools stats {input.bcf} > {output.stats} 2> {log}
        """


rule biallelic_filtered_vcf_multiqc:
    """Collect all stats files from merged vcf files filtered for biallelic sites"""
    input:
        rules.biallelic_filtered_vcf_stats.output.stats,
    output:
        stats="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_biallelic/multiqc/multiqc_report.html",
    params:
        indir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_biallelic/",
        outdir="results/all/vcf/" + REF_NAME + "/stats/vcf_merged_biallelic/multiqc",
    group:
        "biallelic_filtered_vcf_group"
    log:
        "results/logs/9_merge_vcfs/all/" + REF_NAME + "/biallelic_filtered_vcf_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule filter_vcf_missing:
    input:
        bcf=rules.filter_vcf_biallelic.output.bcf,
        index=rules.filter_vcf_biallelic.output.index,
        multiqc=rules.biallelic_filtered_vcf_multiqc.output.stats,
    output:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.genome.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.genome.vcf.gz.csi",
    threads: 2
    params:
        fmiss=config["f_missing"],
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_fmissing{fmiss}_filter_vcf_missing.log",
    singularity:
        bcftools_container
    shell:
        """
        # only include sites with zero missing data
        if [[ `echo 0.0 {params.fmiss} | awk '{{print ($1 == $2)}}'` == 1 ]]
        then
            bcftools view -i 'F_MISSING = {params.fmiss}' -Oz -o {output.vcf} {input.bcf} 2> {log}
        # include all sites
        elif [[ `echo 1.0 {params.fmiss} | awk '{{print ($1 == $2)}}'` == 1 ]]
        then
            bcftools view -i 'F_MISSING <= {params.fmiss}' -Oz -o {output.vcf} {input.bcf} 2> {log}
        # include sites with less than the fraction f_missing of missing data
        elif [[ `echo 0.0 {params.fmiss} 1.0 | awk '{{print ($1 < $2 && $2 < $3)}}'` == 1 ]]
        then 
            bcftools view -i 'F_MISSING < {params.fmiss}' -Oz -o {output.vcf} {input.bcf} 2> {log}
        fi
        
        bcftools index -f {output.vcf} 2>> {log}
        """


rule remove_chromosomes:
    input:
        bcf=rules.filter_vcf_missing.output.vcf,
        index=rules.filter_vcf_missing.output.index,
    output:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.autos.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.autos.vcf.gz.csi",
    threads: 2
    params:
        exclude = ",".join(sexchromosomeList) # parse list with contigs/scaffolds to exclude and convert to format chr1,chr2,chr3 for removal with bcftools view
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_fmissing{fmiss}.autos_remove_chromosomes.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools view {input.bcf} \
        -t ^{params.exclude} \
        -O z -o {output.vcf}

        bcftools index -f {output.vcf} 2>> {log}
        """


rule filtered_vcf2bed:
    """Convert the VCF file after removal of missing data (and optionally sex chromosomes) to BED file containing the remaining sites"""
    input:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    output:
        bed="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.bed",
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".all_fmissing{fmiss}.{chr}_filtered_vcf2bed.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        gzip -cd {input.vcf} | grep -v "^#" | awk -F'\t' '{{print $1, $2-1, $2}}' OFS='\t' > {output.bed} 2> {log}
        """


rule extract_historical_samples:
    input:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
        bed=rules.filtered_vcf2bed.output.bed,
    output:
        vcf="results/historical/vcf/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/historical/vcf/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".historical_fmissing{fmiss}.{chr}_extract_historical_samples.log",
    singularity:
        bcftools_container
    params:
        samples=hist_sm,
        all_samples=ALL_SAMPLES,
    shell:
        """
        samples_edited=`echo {params.samples} | sed 's/ /,/g'`
        samples_len=`echo {params.samples} | wc -w` # count the number of historical samples
        all_samples_len=`echo {params.all_samples} | wc -w` # count the number of all samples

        if [ $samples_len != $all_samples_len ]
        then
            bcftools view -Oz -s $samples_edited -o {output.vcf} {input.vcf} 2> {log} &&
            bcftools index -f {output.vcf} 2>> {log}
        else
            cp {input.vcf} {output.vcf} && touch {output.vcf} 2> {log} &&
            bcftools index -f {output.vcf} 2>> {log}
            echo "Only historical samples present. Copying the input vcf file." >> {log}
        fi
        """


rule extract_modern_samples:
    input:
        vcf="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/all/vcf/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
        bed=rules.filtered_vcf2bed.output.bed,
    output:
        vcf="results/modern/vcf/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/modern/vcf/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".modern_fmissing{fmiss}.{chr}_extract_modern_samples.log",
    singularity:
        bcftools_container
    params:
        samples=mod_sm,
        all_samples=ALL_SAMPLES,
    shell:
        """
        samples_edited=`echo {params.samples} | sed 's/ /,/g'`
        samples_len=`echo {params.samples} | wc -w` # count the number of modern samples
        all_samples_len=`echo {params.all_samples} | wc -w` # count the number of all samples

        if [ $samples_len != $all_samples_len ]
        then
            bcftools view -Oz -s $samples_edited -o {output.vcf} {input.vcf} 2> {log} &&
            bcftools index -f {output.vcf} 2>> {log}
        else
            cp {input.vcf} {output.vcf} && touch {output.vcf} 2> {log} &&
            bcftools index -f {output.vcf} 2>> {log}
            echo "Only modern samples present. Copying the input vcf file." >> {log}
        fi
        """


rule missingness_filtered_vcf_stats:
    """Obtain summary stats of merged vcf file"""
    input:
        merged="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    output:
        stats="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_merged_missing/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
    log:
        "results/logs/9_merge_vcfs/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}_missingness_filtered_vcf_stats.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools stats {input.merged} > {output.stats} 2> {log}
        """


rule missingness_filtered_vcf_multiqc:
    """Collect all stats files from merged vcf files filtered for missing data"""
    input:
        missingness_filtered_vcf_multiqc_inputs,
    output:
        stats=report(
            "results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_merged_missing/multiqc/multiqc_report.html",
            caption="../report/missingness_filtered_vcf_multiqc.rst",
            category="VCF file processing",),
    params:
        indir="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_merged_missing/",
        outdir="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_merged_missing/multiqc",
    log:
        "results/logs/9_merge_vcfs/{dataset}/" + REF_NAME + "/missingness_filtered_vcf_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


# Rules to filter individual BCF files for biallelic sites, missingness and sex-chromosomal scaffolds
# for snpEff and GERP steps. Only triggered when snpEff or GERP is run.

rule repmasked_bcf2vcf:
    """Convert bcf format to vcf.gz for removal of sites"""
    input:
        bcf="results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.bcf",
        index="results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.bcf.csi",
    output:
        vcf=temp("results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.tmp.vcf.gz"),
    log:
        "results/logs/9_merge_vcfs/{dataset}/" + REF_NAME + "/{sample}.{filtered}_repmasked_bcf2vcf.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools convert -O z -o {output.vcf} {input.bcf} 2> {log}
        """


rule filter_biallelic_missing_vcf:
    """Keep only sites with certain upper fraction missingness as specified in config file and sites that are biallelic across all samples (and optionally autosomes) in individual vcf files"""
    input:
        vcf=rules.repmasked_bcf2vcf.output.vcf,
        bed=rules.filtered_vcf2bed.output.bed,
        genomefile=rules.genome_file.output.genomefile,
    output:
        filtered="results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
    threads: 6
    log:
        "results/logs/9_merge_vcfs/{dataset}/" + REF_NAME + "/{sample}.{filtered}_fmissing{fmiss}.{chr}_filter_biallelic_missing_vcf.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools intersect -a {input.vcf} -b {input.bed} -header -sorted -g {input.genomefile} | bgzip -c > {output.filtered} 2> {log}
        """


rule biallelic_missing_filtered_vcf_stats:
    """Obtain summary stats of filtered vcf file"""
    input:
        filtered="results/{dataset}/vcf/" + REF_NAME + "/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
    output:
        stats="results/{dataset}/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/{sample}.{filtered}.snps5.noIndel.QUAL30.dp.AB.repma.biallelic.fmissing{fmiss}.{chr}.vcf.stats.txt",
    log:
        "results/logs/9_merge_vcfs/{dataset}/" + REF_NAME + "/{sample}.{filtered}_fmissing{fmiss}.{chr}_biallelic_missing_filtered_vcf_stats.log",
    singularity:
        bcftools_container
    shell:
        """
        bcftools stats {input.filtered} > {output.stats} 2> {log}
        """


rule historical_biallelic_missing_filtered_vcf_multiqc:
    """Collect all stats files from historical vcf files filtered for biallelic sites and missing data (and optionally sex chromosomes)"""
    input:
        historical_biallelic_missing_filtered_vcf_multiqc_inputs,
    output:
        stats="results/historical/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/multiqc/multiqc_report.html",
    params:
        indir="results/historical/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/",
        outdir="results/historical/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/multiqc",
    log:
        "results/logs/9_merge_vcfs/historical/" + REF_NAME + "/biallelic_missing_{chr}_filtered_vcf_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """


rule modern_biallelic_missing_filtered_vcf_multiqc:
    """Collect all stats files from modern vcf files filtered for biallelic sites and missing data (and optionally sex chromosomes)"""
    input:
        modern_biallelic_missing_filtered_vcf_multiqc_inputs,
    output:
        stats="results/modern/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/multiqc/multiqc_report.html",
    params:
        indir="results/modern/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/",
        outdir="results/modern/vcf/" + REF_NAME + "/stats/vcf_biallelic_missing_{chr}/multiqc",
    log:
        "results/logs/9_merge_vcfs/modern/" + REF_NAME + "/biallelic_missing_{chr}_filtered_vcf_multiqc.log",
    singularity:
        multiqc_container
    shell:
        """
        multiqc -f {params.indir} -o {params.outdir} 2> {log}
        """