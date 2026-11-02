##########################################################################
### 11. Estimate runs of homozygosity (ROHs) in plink based on merged BCF files

# Code collecting output files from this part of the pipeline
ROH_outputs=[]

if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
    ROH_outputs.append(expand("results/all/ROH/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_plot.pdf",
        fmiss=config["f_missing"],
        chr=CHR,
        homsnp=config["homozyg-snp"],
        homkb=config["homozyg-kb"],
        homwinsnp=config["homozyg-window-snp"],
        homwinhet=config["homozyg-window-het"],
        homwinmis=config["homozyg-window-missing"],
        homhet=config["homozyg-het"],))

elif os.path.exists(config["historical_samples"]):
    ROH_outputs.append(expand("results/historical/ROH/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_plot.pdf",
        fmiss=config["f_missing"],
        chr=CHR,
        homsnp=config["homozyg-snp"],
        homkb=config["homozyg-kb"],
        homwinsnp=config["homozyg-window-snp"],
        homwinhet=config["homozyg-window-het"],
        homwinmis=config["homozyg-window-missing"],
        homhet=config["homozyg-het"],))

elif os.path.exists(config["modern_samples"]):
    ROH_outputs.append(expand("results/modern/ROH/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_plot.pdf",
        fmiss=config["f_missing"],
        chr=CHR,
        homsnp=config["homozyg-snp"],
        homkb=config["homozyg-kb"],
        homwinsnp=config["homozyg-window-snp"],
        homwinhet=config["homozyg-window-het"],
        homwinmis=config["homozyg-window-missing"],
        homhet=config["homozyg-het"],))


# snakemake rules
#localrules: FROH_min_2Mb,


rule filter_vcf_hwe:
    """
    Filter the VCF files for ROH analysis
    """
    input:
        vcf="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    output:
        vcf=temp("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.recode.vcf"),
        compressed=temp("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.recode.vcf.gz"),
        index=temp("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.recode.vcf.gz.tbi"),
    threads: 2
    params:
        out="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05",
    log:
        "results/logs/11_ROH/{dataset}/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}_filter_vcf_hwe.log",
    container:
        bcftools_vcftools_container
    shell:
        """
        vcftools --gzvcf {input.vcf} --hwe 0.05 --recode --recode-INFO-all --out {params.out} 2> {log} &&
        bcftools view --threads {threads} -Oz -o {output.compressed} {output.vcf} 2>> {log} &&
        bcftools index -f -t {output.compressed} 2>> {log}
        """


rule ROHs:
    """
    Estimate runs of homozygosity for each individual in the merged BCF file using plink1.9
    """
    input:
        vcf=rules.filter_vcf_hwe.output.compressed,
        index=rules.filter_vcf_hwe.output.index,
    output:
        bed=temp("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.bed"),
        bim=temp("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.bim"),
        fam=temp("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.fam"),
        nosex=temp("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.nosex"),
        roh="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.hom",
        indiv="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.hom.indiv",
        summary="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.hom.summary",
    params:
        homsnp=config["homozyg-snp"],  # Min SNP count per ROH.
        homkb=config["homozyg-kb"],  # Min length of ROH, with min SNP count.
        homwinsnp=config["homozyg-window-snp"],  # Scanning window size.
        homwinhet=config["homozyg-window-het"],  # Max hets in scanning window hit.
        homwinmis=config["homozyg-window-missing"],  # Max missing calls in scanning window hit.
        homhet=config["homozyg-het"],  # By default, a ROH can contain an unlimited number of heterozygous calls; you can impose a limit with --homozyg-het. (This flag was silently ignored by PLINK 1.07.)
        roh="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}",
    log:
        "results/logs/11_ROH/{dataset}/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}_ROHs.log",
    container:
        plink_container
    shell:
        """
        # Convert VCF to plink format
        plink --vcf {input.vcf} --make-bed --allow-extra-chr --out {params.roh} 2> {log} &&

        # Calculate ROH
        plink --bfile {params.roh} --homozyg --homozyg-window-threshold 0.05 --allow-extra-chr \
        --homozyg-snp {params.homsnp} --homozyg-kb {params.homkb} --homozyg-window-snp {params.homwinsnp} \
        --homozyg-window-het {params.homwinhet} --homozyg-window-missing {params.homwinmis} --homozyg-het {params.homhet} --out {params.roh} 2>> {log}
        """


# Code to collect output files for reports
all_ROH_outputs = []

if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
    both_ROH_outputs = expand("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.{filetype}",
        dataset=["modern", "historical"],
        fmiss=config["f_missing"],
        chr=CHR,
        homsnp=config["homozyg-snp"],
        homkb=config["homozyg-kb"],
        homwinsnp=config["homozyg-window-snp"],
        homwinhet=config["homozyg-window-het"],
        homwinmis=config["homozyg-window-missing"],
        homhet=config["homozyg-het"],
        filetype=["hom", "hom.indiv"],)
    all_ROH_outputs.append(both_ROH_outputs)

elif os.path.exists(config["historical_samples"]):
    historical_ROH_outputs = expand("results/historical/ROH/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.{filetype}",
        fmiss=config["f_missing"],
        chr=CHR,
        homsnp=config["homozyg-snp"],
        homkb=config["homozyg-kb"],
        homwinsnp=config["homozyg-window-snp"],
        homwinhet=config["homozyg-window-het"],
        homwinmis=config["homozyg-window-missing"],
        homhet=config["homozyg-het"],
        filetype=["hom", "hom.indiv"],)
    all_ROH_outputs.append(historical_ROH_outputs)

elif os.path.exists(config["modern_samples"]):
    modern_ROH_outputs = expand("results/modern/ROH/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.{filetype}",
        fmiss=config["f_missing"],
        chr=CHR,
        homsnp=config["homozyg-snp"],
        homkb=config["homozyg-kb"],
        homwinsnp=config["homozyg-window-snp"],
        homwinhet=config["homozyg-window-het"],
        homwinmis=config["homozyg-window-missing"],
        homhet=config["homozyg-het"],
        filetype=["hom", "hom.indiv"],)
    all_ROH_outputs.append(modern_ROH_outputs)


rule FROH_min_2Mb:
    """
    Calculate and plot the proportion of the genome in runs of homozygosity, for ROHs >= 2Mb
    """
    input:
        genomefile=rules.index_reference.output.fai,
        ROH=all_ROH_outputs,
    output:
        table=report("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_table.txt",
            caption="../report/FROH_min_2Mb_table.rst",
            category="ROH",),
        plot=report("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_plot.pdf",
            caption="../report/FROH_min_2Mb_plot.rst",
            category="ROH",),
    log:
        "results/logs/11_ROH/{dataset}/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.hwe0.05.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb.log",
    container: shell_container
    shell:
        """
        # Calculate FROH
        python3 workflow/scripts/FROH_min_2Mb_table.py {input.genomefile} {output.table} {input.ROH} 2> {log} &&

        # Plot FROH
        python3 workflow/scripts/FROH_min_2Mb_plot.py {output.table} {output.plot} 2>> {log}
        """
