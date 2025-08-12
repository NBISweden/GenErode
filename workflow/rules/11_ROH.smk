##########################################################################
### 11. Estimate runs of homozygosity (ROHs) in plink based on merged BCF files

# Code collecting output files from this part of the pipeline
ROH_outputs=[]

if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
    ROH_outputs.append(expand("results/all/ROH/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_plot.pdf",
        fmiss=config["f_missing"],
        chr=CHR,
        homsnp=config["homozyg-snp"],
        homkb=config["homozyg-kb"],
        homwinsnp=config["homozyg-window-snp"],
        homwinhet=config["homozyg-window-het"],
        homwinmis=config["homozyg-window-missing"],
        homhet=config["homozyg-het"],))

elif os.path.exists(config["historical_samples"]):
    ROH_outputs.append(expand("results/historical/ROH/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_plot.pdf",
        fmiss=config["f_missing"],
        chr=CHR,
        homsnp=config["homozyg-snp"],
        homkb=config["homozyg-kb"],
        homwinsnp=config["homozyg-window-snp"],
        homwinhet=config["homozyg-window-het"],
        homwinmis=config["homozyg-window-missing"],
        homhet=config["homozyg-het"],))

elif os.path.exists(config["modern_samples"]):
    ROH_outputs.append(expand("results/modern/ROH/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_plot.pdf",
        fmiss=config["f_missing"],
        chr=CHR,
        homsnp=config["homozyg-snp"],
        homkb=config["homozyg-kb"],
        homwinsnp=config["homozyg-window-snp"],
        homwinhet=config["homozyg-window-het"],
        homwinmis=config["homozyg-window-missing"],
        homhet=config["homozyg-het"],))


# snakemake rules
localrules:
    FROH_min_2Mb_plot,

rule vcf2plink_roh:
    """
    Convert VCF files to plink format (version 1.9) for ROH
    """
    input:
        vcf="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    output:
        bed="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.bed",
        bim="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.bim",
        fam="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.fam",
        nosex="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.nosex",
    threads: 2
    params:
        bfile="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}",
    log:
        "results/logs/11_ROH/{dataset}/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}_vcf2plink_roh.log",
    singularity:
        plink_container
    shell:
        """
        plink --vcf {input.vcf} --make-bed --allow-extra-chr --out {params.bfile} 2> {log}
        """


rule ROHs:
    """
    Estimate runs of homozygosity for each individual in the merged BCF file using plink1.9
    """
    input:
        rules.vcf2plink_hwe.output,
    output:
        roh="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.hom",
        indiv="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.hom.indiv",
        summary="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.hom.summary",
    params:
        bfile="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}",
        homsnp=config["homozyg-snp"],  # Min SNP count per ROH.
        homkb=config["homozyg-kb"],  # Min length of ROH, with min SNP count.
        homwinsnp=config["homozyg-window-snp"],  # Scanning window size.
        homwinhet=config["homozyg-window-het"],  # Max hets in scanning window hit.
        homwinmis=config["homozyg-window-missing"],  # Max missing calls in scanning window hit.
        homhet=config["homozyg-het"],  # By default, a ROH can contain an unlimited number of heterozygous calls; you can impose a limit with --homozyg-het. (This flag was silently ignored by PLINK 1.07.)
        roh="results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}",
    log:
        "results/logs/11_ROH/{dataset}/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}_ROHs.log",
    singularity:
        plink_container
    shell:
        """
        plink --bfile {params.bfile} --homozyg --homozyg-window-threshold 0.05 --allow-extra-chr \
        --homozyg-snp {params.homsnp} --homozyg-kb {params.homkb} --homozyg-window-snp {params.homwinsnp} \
        --homozyg-window-het {params.homwinhet} --homozyg-window-missing {params.homwinmis} --homozyg-het {params.homhet} --out {params.roh} 2> {log}
        """


# Code to collect output files for reports
all_ROH_outputs = []

if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
    both_ROH_outputs = expand("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.{filetype}",
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
    historical_ROH_outputs = expand("results/historical/ROH/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.{filetype}",
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
    modern_ROH_outputs = expand("results/modern/ROH/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.{filetype}",
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


rule FROH_min_2Mb_table:
    """
    Calculate the proportion of the genome in runs of homozygosity, for ROHs >= 2Mb
    """
    input:
        genomefile=REF_DIR + "/" + REF_NAME + ".genome",
        ROH=all_ROH_outputs,
    output:
        table=report("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_table.txt",
            caption="../report/FROH_min_2Mb_table.rst",
            category="ROH",),
    log:
        "results/logs/11_ROH/{dataset}/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_table.log",
    script:
        "../scripts/FROH_min_2Mb_table.py"


rule FROH_min_2Mb_plot:
    """Plot the proportion of the genome in runs of homozygosity, for ROHs >= 2Mb"""
    input:
        "results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_table.txt",
    output:
        plot=report("results/{dataset}/ROH/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_plot.pdf",
            caption="../report/FROH_min_2Mb_plot.rst",
            category="ROH",),
    log:
        "results/logs/11_ROH/{dataset}/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.homsnp{homsnp}.homkb{homkb}.homwinsnp{homwinsnp}.homwinhet{homwinhet}.homwinmis{homwinmis}.homhet{homhet}.FROH_min_2Mb_plot.log",
    script:
        "../scripts/FROH_min_2Mb_plot.py"
