##########################################################################
### 7. Run mlRho on filtered BAM files

# Code collecting output files from this part of the pipeline
if os.path.exists(config["historical_samples"]) and os.path.exists(config["modern_samples"]):
    all_outputs.append("results/all/mlRho/" + REF_NAME + ".all.mlRho_theta_plot.pdf")

elif os.path.exists(config["historical_samples"]):
    all_outputs.append("results/historical/mlRho/" + REF_NAME + ".historical.mlRho_theta_plot.pdf")

elif os.path.exists(config["modern_samples"]):
    all_outputs.append("results/modern/mlRho/" + REF_NAME + ".modern.mlRho_theta_plot.pdf")


# Functions used by rules of this part of the pipeline
def bam_file_mlRho(wildcards):
    """Select correct bam file for each sample"""
    if wildcards.sample in HIST_NOT_RESCALED_NOT_SUBSAMPLED_SAMPLES:
        bam = ("results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam".format(sample=wildcards.sample))
    elif wildcards.sample in HIST_NOT_RESCALED_SUBSAMPLED_SAMPLES:
        bam = "results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.bam".format(sample=wildcards.sample, DP=config["subsampling_depth"])
    elif wildcards.sample in HIST_RESCALED_NOT_SUBSAMPLED_SAMPLES:
        bam = "results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.bam".format(sample=wildcards.sample)
    elif wildcards.sample in HIST_RESCALED_SUBSAMPLED_SAMPLES:
        bam = "results/historical/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.bam".format(sample=wildcards.sample, DP=config["subsampling_depth"])
    elif wildcards.sample in MODERN_NOT_SUBSAMPLED_SAMPLES:
        bam = "results/modern/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.bam".format(sample=wildcards.sample)
    elif wildcards.sample in MODERN_SUBSAMPLED_SAMPLES:
        bam = "results/modern/mapping/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.bam".format(sample=wildcards.sample, DP=config["subsampling_depth"])
    return [bam]

def depth_file_mlRho(wildcards):
    """Select correct depth stats file for each sample"""
    if wildcards.sample in HIST_NOT_SUBSAMPLED_SAMPLES:
        dpstats = "results/historical/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample)
    elif wildcards.sample in HIST_NOT_RESCALED_SUBSAMPLED_SAMPLES:
        dpstats = "results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample, DP=config["subsampling_depth"])
    elif wildcards.sample in HIST_RESCALED_SUBSAMPLED_SAMPLES:
        dpstats = "results/historical/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample, DP=config["subsampling_depth"])
    elif wildcards.sample in MODERN_NOT_SUBSAMPLED_SAMPLES:
        dpstats = "results/modern/mapping/" + REF_NAME + "/stats/bams_indels_realigned/{sample}.merged.rmdup.merged.realn.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample)
    elif wildcards.sample in MODERN_SUBSAMPLED_SAMPLES:
        dpstats = "results/modern/mapping/" + REF_NAME + "/stats/bams_subsampled/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.Q30.bam.dpstats.txt".format(sample=wildcards.sample, DP=config["subsampling_depth"])
    return dpstats

def bed_file_autos_mlRho(wildcards):
    """Select correct bed file for filtering during mlRho analysis"""
    if len(sexchromosomeList) > 0:
        if config["CpG_from_vcf"] == True:
            bed = "results/" + REF_NAME + ".noCpG_vcf.repma.autos.bed"
        elif config["CpG_from_reference"] == True:
            bed = "results/" + REF_NAME + ".noCpG_ref.repma.autos.bed"
        elif config["CpG_from_vcf_and_reference"] == True:
            bed = "results/" + REF_NAME + ".noCpG_vcfref.repma.autos.bed"
        else:
            bed = "results/" + REF_NAME + ".repma.autos.bed"
    return bed

def bed_file_sexchr_mlRho(wildcards):
    """Select correct bed file for filtering during mlRho analysis"""
    if len(sexchromosomeList) > 0:
        if config["CpG_from_vcf"] == True:
            bed = "results/" + REF_NAME + ".noCpG_vcf.repma.sexchr.bed"
        elif config["CpG_from_reference"] == True:
            bed = "results/" + REF_NAME + ".noCpG_ref.repma.sexchr.bed"
        elif config["CpG_from_vcf_and_reference"] == True:
            bed = "results/" + REF_NAME + ".noCpG_vcfref.repma.sexchr.bed"
        else:
            bed = "results/" + REF_NAME + ".repma.sexchr.bed"
    return bed

def bed_file_genome_mlRho(wildcards):
    """Select correct bed file for filtering during mlRho analysis"""
    if len(sexchromosomeList) == 0:
        if config["CpG_from_vcf"] == True:
            bed = "results/" + REF_NAME + ".noCpG_vcf.repma.bed"
        elif config["CpG_from_reference"] == True:
            bed = "results/" + REF_NAME + ".noCpG_ref.repma.bed"
        elif config["CpG_from_vcf_and_reference"] == True:
            bed = "results/" + REF_NAME + ".noCpG_vcfref.repma.bed"
        else:
            bed = "results/" + REF_NAME + ".repma.bed"
    return bed

def all_mlRho_outputs(wildcards):
    """Collect output files of this step of the pipeline"""
    outlist = []
    if os.path.exists(config["historical_samples"]):
        if len(sexchromosomeList) > 0:
            if config["mlRho_autosomes_sexchromosomes"] == True:
                rescaled_not_subsampled_not_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.repma.{chr}.mlRho.txt",
                    sample=HIST_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,
                    chr=["autos", "sexchr"],)
                not_rescaled_not_subsampled_not_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.repma.{chr}.mlRho.txt",
                    sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,
                    chr=["autos", "sexchr"],)
                rescaled_subsampled_not_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.repma.{chr}.mlRho.txt",
                    sample=HIST_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
                    DP=config["subsampling_depth"],
                    chr=["autos", "sexchr"],)
                not_rescaled_subsampled_not_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.{chr}.mlRho.txt",
                    sample=HIST_NOT_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
                    DP=config["subsampling_depth"],
                    chr=["autos", "sexchr"],)
                outlist += (rescaled_not_subsampled_not_CpG_chr_mlRho + not_rescaled_not_subsampled_not_CpG_chr_mlRho + rescaled_subsampled_not_CpG_chr_mlRho + not_rescaled_subsampled_not_CpG_chr_mlRho)
                if config["CpG_from_vcf"] == True:
                    rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.noCpG_vcf.repma.{chr}.mlRho.txt",
                        sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                        chr=["autos", "sexchr"],)
                    not_rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcf.repma.{chr}.mlRho.txt",
                        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                        chr=["autos", "sexchr"],)
                    rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.noCpG_vcf.repma.{chr}.mlRho.txt",
                        sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],
                        chr=["autos", "sexchr"],)
                    not_rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcf.repma.{chr}.mlRho.txt",
                        sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],
                        chr=["autos", "sexchr"],)
                    outlist += (rescaled_not_subsampled_CpG_chr_mlRho + not_rescaled_not_subsampled_CpG_chr_mlRho + rescaled_subsampled_CpG_chr_mlRho + not_rescaled_subsampled_CpG_chr_mlRho)
                elif config["CpG_from_reference"] == True:
                    rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.noCpG_ref.repma.{chr}.mlRho.txt",
                        sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                        chr=["autos", "sexchr"],)
                    not_rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_ref.repma.{chr}.mlRho.txt",
                        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                        chr=["autos", "sexchr"],)
                    rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.noCpG_ref.repma.{chr}.mlRho.txt",
                        sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],
                        chr=["autos", "sexchr"],)
                    not_rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_ref.repma.{chr}.mlRho.txt",
                        sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],
                        chr=["autos", "sexchr"],)
                    outlist += (rescaled_not_subsampled_CpG_chr_mlRho + not_rescaled_not_subsampled_CpG_chr_mlRho + rescaled_subsampled_CpG_chr_mlRho + not_rescaled_subsampled_CpG_chr_mlRho)
                elif config["CpG_from_vcf_and_reference"] == True:
                    rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.noCpG_vcfref.repma.{chr}.mlRho.txt",
                        sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                        chr=["autos", "sexchr"],)
                    not_rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcfref.repma.{chr}.mlRho.txt",
                        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,
                        chr=["autos", "sexchr"],)
                    rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.noCpG_vcfref.repma.{chr}.mlRho.txt",
                        sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],
                        chr=["autos", "sexchr"],)
                    not_rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcfref.repma.{chr}.mlRho.txt",
                        sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],
                        chr=["autos", "sexchr"],)
                    outlist += (rescaled_not_subsampled_CpG_chr_mlRho + not_rescaled_not_subsampled_CpG_chr_mlRho + rescaled_subsampled_CpG_chr_mlRho + not_rescaled_subsampled_CpG_chr_mlRho)
            elif config["mlRho_autosomes_sexchromosomes"] == False:
                rescaled_not_subsampled_not_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.repma.autos.mlRho.txt",
                    sample=HIST_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
                not_rescaled_not_subsampled_not_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.repma.autos.mlRho.txt",
                    sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
                rescaled_subsampled_not_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.repma.autos.mlRho.txt",
                    sample=HIST_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                not_rescaled_subsampled_not_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.autos.mlRho.txt",
                    sample=HIST_NOT_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                outlist += (rescaled_not_subsampled_not_CpG_chr_mlRho + not_rescaled_not_subsampled_not_CpG_chr_mlRho + rescaled_subsampled_not_CpG_chr_mlRho + not_rescaled_subsampled_not_CpG_chr_mlRho)
                if config["CpG_from_vcf"] == True:
                    rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.noCpG_vcf.repma.autos.mlRho.txt",
                        sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                    not_rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcf.repma.autos.mlRho.txt",
                        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                    rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.noCpG_vcf.repma.autos.mlRho.txt",
                        sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],)
                    not_rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcf.repma.autos.mlRho.txt",
                        sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],)
                    outlist += (rescaled_not_subsampled_CpG_chr_mlRho + not_rescaled_not_subsampled_CpG_chr_mlRho + rescaled_subsampled_CpG_chr_mlRho + not_rescaled_subsampled_CpG_chr_mlRho)
                elif config["CpG_from_reference"] == True:
                    rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.noCpG_ref.repma.autos.mlRho.txt",
                        sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                    not_rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_ref.repma.autos.mlRho.txt",
                        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                    rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.noCpG_ref.repma.autos.mlRho.txt",
                        sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],)
                    not_rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_ref.repma.autos.mlRho.txt",
                        sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],)
                    outlist += (rescaled_not_subsampled_CpG_chr_mlRho + not_rescaled_not_subsampled_CpG_chr_mlRho + rescaled_subsampled_CpG_chr_mlRho + not_rescaled_subsampled_CpG_chr_mlRho)
                elif config["CpG_from_vcf_and_reference"] == True:
                    rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.noCpG_vcfref.repma.autos.mlRho.txt",
                        sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                    not_rescaled_not_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcfref.repma.autos.mlRho.txt",
                        sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                    rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.noCpG_vcfref.repma.autos.mlRho.txt",
                        sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],)
                    not_rescaled_subsampled_CpG_chr_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcfref.repma.autos.mlRho.txt",
                        sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],)
                    outlist += (rescaled_not_subsampled_CpG_chr_mlRho + not_rescaled_not_subsampled_CpG_chr_mlRho + rescaled_subsampled_CpG_chr_mlRho + not_rescaled_subsampled_CpG_chr_mlRho)
        elif len(sexchromosomeList) == 0:
            rescaled_not_subsampled_not_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.repma.genome.mlRho.txt",
                sample=HIST_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
            not_rescaled_not_subsampled_not_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.repma.genome.mlRho.txt",
                sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
            rescaled_subsampled_not_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.repma.genome.mlRho.txt",
                sample=HIST_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
                DP=config["subsampling_depth"],)
            not_rescaled_subsampled_not_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.genome.mlRho.txt",
                sample=HIST_NOT_RESCALED_SUBSAMPLED_NOT_CpG_SAMPLES,
                DP=config["subsampling_depth"],)
            outlist += (rescaled_not_subsampled_not_CpG_mlRho + not_rescaled_not_subsampled_not_CpG_mlRho + rescaled_subsampled_not_CpG_mlRho + not_rescaled_subsampled_not_CpG_mlRho)
            if config["CpG_from_vcf"] == True:
                rescaled_not_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.noCpG_vcf.repma.genome.mlRho.txt",
                    sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                not_rescaled_not_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcf.repma.genome.mlRho.txt",
                    sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                rescaled_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.noCpG_vcf.repma.genome.mlRho.txt",
                    sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                not_rescaled_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcf.repma.genome.mlRho.txt",
                    sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                outlist += (rescaled_not_subsampled_CpG_mlRho + not_rescaled_not_subsampled_CpG_mlRho + rescaled_subsampled_CpG_mlRho + not_rescaled_subsampled_CpG_mlRho)
            elif config["CpG_from_reference"] == True:
                rescaled_not_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.noCpG_ref.repma.genome.mlRho.txt",
                    sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                not_rescaled_not_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_ref.repma.genome.mlRho.txt",
                    sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                rescaled_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.noCpG_ref.repma.genome.mlRho.txt",
                    sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                not_rescaled_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_ref.repma.genome.mlRho.txt",
                    sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                outlist += (rescaled_not_subsampled_CpG_mlRho + not_rescaled_not_subsampled_CpG_mlRho + rescaled_subsampled_CpG_mlRho + not_rescaled_subsampled_CpG_mlRho)
            elif config["CpG_from_vcf_and_reference"] == True:
                rescaled_not_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.noCpG_vcfref.repma.genome.mlRho.txt",
                    sample=HIST_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                not_rescaled_not_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcfref.repma.genome.mlRho.txt",
                    sample=HIST_NOT_RESCALED_NOT_SUBSAMPLED_CpG_SAMPLES,)
                rescaled_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.rescaled.mapped_q30.subs_dp{DP}.noCpG_vcfref.repma.genome.mlRho.txt",
                    sample=HIST_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                not_rescaled_subsampled_CpG_mlRho = expand("results/historical/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcfref.repma.genome.mlRho.txt",
                    sample=HIST_NOT_RESCALED_SUBSAMPLED_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                outlist += (rescaled_not_subsampled_CpG_mlRho + not_rescaled_not_subsampled_CpG_mlRho + rescaled_subsampled_CpG_mlRho + not_rescaled_subsampled_CpG_mlRho)
    if os.path.exists(config["modern_samples"]):
        if len(sexchromosomeList) > 0:
            if config["mlRho_autosomes_sexchromosomes"] == True:
                not_subsampled_not_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.repma.{chr}.mlRho.txt",
                    sample=MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,
                    chr=["autos", "sexchr"],)
                subsampled_not_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.{chr}.mlRho.txt",
                    sample=MODERN_SUBSAMPLED_NOT_CpG_SAMPLES,
                    DP=config["subsampling_depth"],
                    chr=["autos", "sexchr"],)
                outlist += (not_subsampled_not_CpG_chr_mlRho + subsampled_not_CpG_chr_mlRho)
                if config["CpG_from_vcf"] == True:
                    not_subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcf.repma.{chr}.mlRho.txt",
                        sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,
                        chr=["autos", "sexchr"],)
                    subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcf.repma.{chr}.mlRho.txt",
                        sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],
                        chr=["autos", "sexchr"],)
                    outlist += (not_subsampled_CpG_chr_mlRho + subsampled_CpG_chr_mlRho)
                elif config["CpG_from_reference"] == True:
                    not_subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_ref.repma.{chr}.mlRho.txt",
                        sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,
                        chr=["autos", "sexchr"],)
                    subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_ref.repma.{chr}.mlRho.txt",
                        sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],
                        chr=["autos", "sexchr"],)
                    outlist += (not_subsampled_CpG_chr_mlRho + subsampled_CpG_chr_mlRho)
                elif config["CpG_from_vcf_and_reference"] == True:
                    not_subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcfref.repma.{chr}.mlRho.txt",
                        sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,
                        chr=["autos", "sexchr"],)
                    subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcfref.repma.{chr}.mlRho.txt",
                        sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],
                        chr=["autos", "sexchr"],)
                    outlist += (not_subsampled_CpG_chr_mlRho + subsampled_CpG_chr_mlRho)
            elif config["mlRho_autosomes_sexchromosomes"] == False:
                not_subsampled_not_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.repma.autos.mlRho.txt",
                    sample=MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
                subsampled_not_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.autos.mlRho.txt",
                    sample=MODERN_SUBSAMPLED_NOT_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                outlist += (not_subsampled_not_CpG_chr_mlRho + subsampled_not_CpG_chr_mlRho)
                if config["CpG_from_vcf"] == True:
                    not_subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcf.repma.autos.mlRho.txt",
                        sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
                    subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcf.repma.autos.mlRho.txt",
                        sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],)
                    outlist += (not_subsampled_CpG_chr_mlRho + subsampled_CpG_chr_mlRho)
                elif config["CpG_from_reference"] == True:
                    not_subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_ref.repma.autos.mlRho.txt",
                        sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
                    subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_ref.repma.autos.mlRho.txt",
                        sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],)
                    outlist += (not_subsampled_CpG_chr_mlRho + subsampled_CpG_chr_mlRho)
                elif config["CpG_from_vcf_and_reference"] == True:
                    not_subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcfref.repma.autos.mlRho.txt",
                        sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
                    subsampled_CpG_chr_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcfref.repma.autos.mlRho.txt",
                        sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                        DP=config["subsampling_depth"],)
                    outlist += (not_subsampled_CpG_chr_mlRho + subsampled_CpG_chr_mlRho)
        elif len(sexchromosomeList) == 0:
            not_subsampled_not_CpG_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.repma.genome.mlRho.txt",
                sample=MODERN_NOT_SUBSAMPLED_NOT_CpG_SAMPLES,)
            subsampled_not_CpG_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.repma.genome.mlRho.txt",
                sample=MODERN_SUBSAMPLED_NOT_CpG_SAMPLES,
                DP=config["subsampling_depth"],)
            outlist += (not_subsampled_not_CpG_mlRho + subsampled_not_CpG_mlRho)
            if config["CpG_from_vcf"] == True:
                not_subsampled_CpG_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcf.repma.genome.mlRho.txt",
                    sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
                subsampled_CpG_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcf.repma.genome.mlRho.txt",
                    sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                outlist += (not_subsampled_CpG_mlRho + subsampled_CpG_mlRho)
            elif config["CpG_from_reference"] == True:
                not_subsampled_CpG_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_ref.repma.genome.mlRho.txt",
                    sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
                subsampled_CpG_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_ref.repma.genome.mlRho.txt",
                    sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                outlist += (not_subsampled_CpG_mlRho + subsampled_CpG_mlRho)
            elif config["CpG_from_vcf_and_reference"] == True:
                not_subsampled_CpG_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.noCpG_vcfref.repma.genome.mlRho.txt",
                    sample=MODERN_NOT_SUBSAMPLED_CpG_SAMPLES,)
                subsampled_CpG_mlRho = expand("results/modern/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.realn.mapped_q30.subs_dp{DP}.noCpG_vcfref.repma.genome.mlRho.txt",
                    sample=MODERN_SUBSAMPLED_CpG_SAMPLES,
                    DP=config["subsampling_depth"],)
                outlist += (not_subsampled_CpG_mlRho + subsampled_CpG_mlRho)
    return outlist


# snakemake rules
localrules:
    mlRho_table,
    mlRho_theta_plot,


rule bam2pro_autos:
    """Generate pro files from bam files"""
    """Note that the depth filter is recalculated for subsampled bam files, according to the target depth for subsampling"""
    input:
        bam=bam_file_mlRho,
        dp=depth_file_mlRho,
        bed=bed_file_autos_mlRho,
    output:
        pro=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.autos.pro"),
    log:
        "results/logs/7_mlRho/{dataset}/" + REF_NAME + "/{sample}.{processed}_bam2pro_autos.log",
    singularity:
        "oras://community.wave.seqera.io/library/mlrho_samtools:f1b883982ae4586d"
    shell:
        """
        minDP=`head -n 1 {input.dp} | cut -d' ' -f 2`
        maxDP=`head -n 1 {input.dp} | cut -d' ' -f 3`

        # check minimum depth threshold
        if awk "BEGIN{{exit ! ($minDP < 3)}}"
        then
          minDP=3
        fi

        samtools mpileup -q 30 -Q 30 -B -l {input.bed} {input.bam[0]} | awk -v minDP="$minDP" -v maxDP="$maxDP" '$4 >=minDP && $4 <=maxDP' | \
        sam2pro -c 5 > {output.pro} 2> {log}
        """


rule mlRho_autos:
    """Format the pro file and run mlRho"""
    """Note that the depth filter is recalculated for subsampled bam files, according to the target depth for subsampling"""
    input:
        pro=rules.bam2pro_autos.output,
        dp=depth_file_mlRho,
    output:
        mlRho="results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.autos.mlRho.txt",
        con=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.autos_profileDb.con"),
        lik=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.autos_profileDb.lik"),
        pos=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.autos_profileDb.pos"),
        sum=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.autos_profileDb.sum"),
    params:
        db="results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.autos_profileDb",
    log:
        "results/logs/7_mlRho/{dataset}/" + REF_NAME + "/{sample}.{processed}_mlRho_autos.log",
    singularity:
        "oras://community.wave.seqera.io/library/mlrho_samtools:f1b883982ae4586d"
    shell:
        """
        minDP=`head -n 1 {input.dp} | cut -d' ' -f 2`

        # check minimum depth threshold
        if awk "BEGIN{{exit ! ($minDP < 3)}}"
        then
          minDP=3
        fi

        # Further format the pro file
        formatPro -c $minDP -n {params.db} {input.pro} 2> {log} &&

        # run mlRho
        mlRho -M 0 -I -n {params.db} > {output.mlRho} 2>> {log}
        """


rule bam2pro_sexchr:
    """Generate pro files from bam files"""
    """Note that the depth filter is recalculated for subsampled bam files, according to the target depth for subsampling"""
    input:
        bam=bam_file_mlRho,
        dp=depth_file_mlRho,
        bed=bed_file_sexchr_mlRho,
    output:
        pro=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.sexchr.pro"),
    log:
        "results/logs/7_mlRho/{dataset}/" + REF_NAME + "/{sample}.{processed}_bam2pro_sexchr.log",
    singularity:
        "oras://community.wave.seqera.io/library/mlrho_samtools:f1b883982ae4586d"
    shell:
        """
        minDP=`head -n 1 {input.dp} | cut -d' ' -f 2`
        maxDP=`head -n 1 {input.dp} | cut -d' ' -f 3`

        # check minimum depth threshold
        if awk "BEGIN{{exit ! ($minDP < 3)}}"
        then
          minDP=3
        fi

        samtools mpileup -q 30 -Q 30 -B -l {input.bed} {input.bam[0]} | awk -v minDP="$minDP" -v maxDP="$maxDP" '$4 >=minDP && $4 <=maxDP' | \
        sam2pro -c 5 > {output.pro} 2> {log}
        """


rule mlRho_sexchr:
    """Format the pro file and run mlRho"""
    """Note that the depth filter is recalculated for subsampled bam files, according to the target depth for subsampling"""
    input:
        pro=rules.bam2pro_sexchr.output,
        dp=depth_file_mlRho,
    output:
        mlRho="results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.sexchr.mlRho.txt",
        con=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.sexchr_profileDb.con"),
        lik=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.sexchr_profileDb.lik"),
        pos=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.sexchr_profileDb.pos"),
        sum=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.sexchr_profileDb.sum"),
    params:
        db="results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.sexchr_profileDb",
    log:
        "results/logs/7_mlRho/{dataset}/" + REF_NAME + "/{sample}.{processed}_mlRho_sexchr.log",
    singularity:
        "oras://community.wave.seqera.io/library/mlrho_samtools:f1b883982ae4586d"
    shell:
        """
        minDP=`head -n 1 {input.dp} | cut -d' ' -f 2`

        # check minimum depth threshold
        if awk "BEGIN{{exit ! ($minDP < 3)}}"
        then
          minDP=3
        fi

        # Further format the pro file
        formatPro -c $minDP -n {params.db} {input.pro} 2> {log} &&

        # run mlRho
        mlRho -M 0 -I -n {params.db} > {output.mlRho} 2>> {log}
        """


rule bam2pro_genome:
    """Generate pro files from bam files"""
    """Note that the depth filter is recalculated for subsampled bam files, according to the target depth for subsampling"""
    input:
        bam=bam_file_mlRho,
        dp=depth_file_mlRho,
        bed=bed_file_genome_mlRho,
    output:
        pro=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.genome.pro"),
    log:
        "results/logs/7_mlRho/{dataset}/" + REF_NAME + "/{sample}.{processed}_bam2pro_genome.log",
    singularity:
        "oras://community.wave.seqera.io/library/mlrho_samtools:f1b883982ae4586d"
    shell:
        """
        minDP=`head -n 1 {input.dp} | cut -d' ' -f 2`
        maxDP=`head -n 1 {input.dp} | cut -d' ' -f 3`

        # check minimum depth threshold
        if awk "BEGIN{{exit ! ($minDP < 3)}}"
        then
          minDP=3
        fi

        samtools mpileup -q 30 -Q 30 -B -l {input.bed} {input.bam[0]} | awk -v minDP="$minDP" -v maxDP="$maxDP" '$4 >=minDP && $4 <=maxDP' | \
        sam2pro -c 5 > {output.pro} 2> {log}
        """


rule mlRho_genome:
    """Format the pro file and run mlRho"""
    """Note that the depth filter is recalculated for subsampled bam files, according to the target depth for subsampling"""
    input:
        pro=rules.bam2pro_genome.output,
        dp=depth_file_mlRho,
    output:
        mlRho="results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.genome.mlRho.txt",
        con=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.genome_profileDb.con"),
        lik=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.genome_profileDb.lik"),
        pos=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.genome_profileDb.pos"),
        sum=temp("results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.genome_profileDb.sum"),
    params:
        db="results/{dataset}/mlRho/" + REF_NAME + "/{sample}.merged.rmdup.merged.{processed}.genome_profileDb",
    log:
        "results/logs/7_mlRho/{dataset}/" + REF_NAME + "/{sample}.{processed}_mlRho_genome.log",
    singularity:
        "oras://community.wave.seqera.io/library/mlrho_samtools:f1b883982ae4586d"
    shell:
        """
        minDP=`head -n 1 {input.dp} | cut -d' ' -f 2`

        # check minimum depth threshold
        if awk "BEGIN{{exit ! ($minDP < 3)}}"
        then
          minDP=3
        fi

        # Further format the pro file
        formatPro -c $minDP -n {params.db} {input.pro} 2> {log} &&

        # run mlRho
        mlRho -M 0 -I -n {params.db} > {output.mlRho} 2>> {log}
        """


rule mlRho_table:
    """
    Concatenate mlRho output files for all samples
    """
    input:
        all_mlRho_outputs,
    output:
        table=report("results/{dataset}/mlRho/" + REF_NAME + ".{dataset}.mlRho_table.txt",
            caption="../report/mlRho_table.rst",
            category="mlRho",),
    log:
        "results/logs/7_mlRho/{dataset}/" + REF_NAME + ".{dataset}.mlRho_table.log",
    message:
        "The following files are concatenated: {input}"
    script:
        "../scripts/mlRho_table.py"


rule mlRho_theta_plot:
    """
    Plot theta for all samples plus confidence intervals
    """
    input:
        table="results/{dataset}/mlRho/" + REF_NAME + ".{dataset}.mlRho_table.txt",
    output:
        plot=report("results/{dataset}/mlRho/" + REF_NAME + ".{dataset}.mlRho_theta_plot.pdf",
            caption="../report/mlRho_theta_plot.rst",
            category="mlRho",),
    log:
        "results/logs/7_mlRho/{dataset}/" + REF_NAME + ".{dataset}.mlRho_theta_plot.log",
    script:
        "../scripts/mlRho_theta_plot.py"
