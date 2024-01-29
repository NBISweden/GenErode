##########################################################################
### 10. Plot PCAs for historical, modern, and all samples combined

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]):
    all_outputs.append(expand("results/historical/pca/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.pc1_pc2.pdf",
        fmiss=config["f_missing"],
        chr=CHR,))
    all_outputs.append(expand("results/historical/pca/" + REF_NAME + ".historical.merged.biallelic.fmissing{fmiss}.{chr}.pc1_pc3.pdf",
        fmiss=config["f_missing"],
        chr=CHR,))

if os.path.exists(config["modern_samples"]):
    all_outputs.append(expand("results/modern/pca/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.pc1_pc2.pdf",
        fmiss=config["f_missing"],
        chr=CHR,))
    all_outputs.append(expand("results/modern/pca/" + REF_NAME + ".modern.merged.biallelic.fmissing{fmiss}.{chr}.pc1_pc3.pdf",
        fmiss=config["f_missing"],
        chr=CHR,))

if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
    all_outputs.append(expand("results/all/pca/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.pc1_pc2.pdf",
        fmiss=config["f_missing"],
        chr=CHR,))
    all_outputs.append(expand("results/all/pca/" + REF_NAME + ".all.merged.biallelic.fmissing{fmiss}.{chr}.pc1_pc3.pdf",
        fmiss=config["f_missing"],
        chr=CHR,))


# snakemake rules
localrules:
    plot_pc1_pc2,
    plot_pc1_pc3,


rule vcf2plink_pca:
    """Convert VCF files to plink format (version 1.9) for PCA"""
    input:
        vcf="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz",
        index="results/{dataset}/vcf/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.vcf.gz.csi",
    output:
        bed=temp("results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.bed"),
        bim=temp("results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.bim"),
        fam=temp("results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.fam"),
        nosex=temp("results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.nosex"),
    threads: 2
    params:
        bfile="results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}",
    log:
        "results/logs/10_pca/{dataset}/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}_vcf2plink_pca.log",
    singularity:
        "docker://quay.io/biocontainers/plink:1.90b6.12--heea4ae3_0"
    shell:
        """
        plink --vcf {input.vcf} --make-bed --allow-extra-chr --out {params.bfile} 2> {log}
        """


rule plink_eigenvec:
    """Create one PCA per BCF file to find e.g. batch effects"""
    input:
        bed=rules.vcf2plink_pca.output.bed,
        bim=rules.vcf2plink_pca.output.bim,
        fam=rules.vcf2plink_pca.output.fam,
        nosex=rules.vcf2plink_pca.output.nosex,
    output:
        eigenvec="results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.eigenvec",
        eigenval="results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.eigenval",
    params:
        bfile="results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}",
    log:
        "results/logs/10_pca/{dataset}/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}_plink_eigenvec.log",
    singularity:
        "docker://quay.io/biocontainers/plink:1.90b6.12--heea4ae3_0"
    shell:
        """
        samples=`cat {input.fam} | wc -l`
        if [ "$samples" -gt 1 ]
        then
          plink --bfile {params.bfile} --allow-extra-chr --pca --out {params.bfile} 2> {log}
        else
          touch {output.eigenvec} && touch {output.eigenval} 2> {log}
          echo "Not enough samples to calculate a PCA." >> {log}
        fi
        """


rule plot_pc1_pc2:
    """Create a pdf of the PCA axes PC1 and PC2"""
    input:
        eigenvec=rules.plink_eigenvec.output.eigenvec,
        eigenval=rules.plink_eigenvec.output.eigenval,
    output:
        pdf=report("results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.pc1_pc2.pdf",
            caption="../report/pca_plot_pc1_pc2.rst",
            category="PCA",),
    log:
        "results/logs/10_pca/{dataset}/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}_plot_pc1_pc2.log",
    script:
        "../scripts/pc1_vs_pc2_plot.py"


rule plot_pc1_pc3:
    """Create a pdf of the PCA axes PC1 and PC3"""
    input:
        eigenvec=rules.plink_eigenvec.output.eigenvec,
        eigenval=rules.plink_eigenvec.output.eigenval,
    output:
        pdf=report("results/{dataset}/pca/" + REF_NAME + ".{dataset}.merged.biallelic.fmissing{fmiss}.{chr}.pc1_pc3.pdf",
            caption="../report/pca_plot_pc1_pc3.rst",
            category="PCA",),
    log:
        "results/logs/10_pca/{dataset}/" + REF_NAME + ".{dataset}_fmissing{fmiss}.{chr}_plot_pc1_pc3.log",
    script:
        "../scripts/pc1_vs_pc3_plot.py"
