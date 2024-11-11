##########################################################################
### 0.2 Repeat prediction and repeat masking of the reference genome

# Code collecting output files from this part of the pipeline
all_outputs.append(REF_DIR + "/" + REF_NAME + ".repeats.sorted.bed")
all_outputs.append(REF_DIR + "/" + REF_NAME + ".repma.bed")


# snakemake rules
rule ref_upper:
    """Reference assembly file preparation for RepeatModeler and RepeatMasker: change all lowercase bases to uppercase"""
    input:
        ref=config["ref_path"],
    output:
        ref_upper=REF_DIR + "/" + REF_NAME + ".upper.fasta",
    log:
        "results/logs/0.2_repeat_identification/" + REF_NAME + "_ref_upper.log",
    shell:
        """
        awk '{{ if ($0 !~ />/) {{print toupper($0)}} else {{print $0}} }}' {input.ref} > {output.ref_upper} 2> {log}
        """


rule repeatmodeler:
    """RepeatModeler for de novo repeat prediction from a reference assembly"""
    input:
        ref_upper=rules.ref_upper.output,
    output:
        repmo=REF_DIR + "/repeatmodeler/" + REF_NAME + "/RM_raw.out/consensi.fa.classified",
        stk=REF_DIR + "/repeatmodeler/" + REF_NAME + "/RM_raw.out/families-classified.stk",
    params:
        dir=REF_DIR + "/repeatmodeler/" + REF_NAME + "/",
        name=REF_NAME,
        ref_upper="../../" + REF_NAME + ".upper.fasta",
        abs_tmp=os.path.abspath("tmpConsensi.fa"),
    log:
        os.path.abspath("results/logs/0.2_repeat_identification/" + REF_NAME + "_repeatmodeler.log"),
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/repeatmodeler:2.0.4--pl5321hdfd78af_0"
    shell:
        """
        cd {params.dir}

        # Build repeat database
        BuildDatabase -engine ncbi -name {params.name} {params.ref_upper} 2> {log} &&

        # Run RepeatModeler
        RepeatModeler -engine ncbi -threads {threads} -database {params.name} 2>> {log} &&

        # copy the output files to a new directory
        cp RM_*.*/consensi.fa.classified RM_raw.out/ 2>> {log} &&
        cp RM_*.*/families-classified.stk RM_raw.out/ 2>> {log}

        # remove temporary file
        if [ -f {params.abs_tmp} ]
        then
          rm {params.abs_tmp} 2>> {log}
        fi
        """


rule repeatmasker:
    """Repeat mask the full genome assembly using raw de novo predicted repeats"""
    input:
        ref_upper=rules.ref_upper.output,
        repmo=rules.repeatmodeler.output.repmo,
    output:
        rep_masked=REF_DIR + "/repeatmasker/" + REF_NAME + "/" + REF_NAME + ".upper.fasta.masked",
        rep_tbl=REF_DIR + "/repeatmasker/" + REF_NAME + "/" + REF_NAME + ".upper.fasta.tbl",
        rep_out=REF_DIR + "/repeatmasker/" + REF_NAME + "/" + REF_NAME + ".upper.fasta.out",
        rep_cat=REF_DIR + "/repeatmasker/" + REF_NAME + "/" + REF_NAME + ".upper.fasta.cat.gz",
    params:
        dir=REF_DIR + "/repeatmasker/" + REF_NAME + "/",
        repmo="../../repeatmodeler/" + REF_NAME + "/RM_raw.out/consensi.fa.classified",
        ref_upper="../../" + REF_NAME + ".upper.fasta",
        rep_cat_unzip=REF_NAME + ".upper.fasta.cat",
    log:
        os.path.abspath("results/logs/0.2_repeat_identification/" + REF_NAME + "_repeatmasker.log"),
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/repeatmodeler:2.0.4--pl5321hdfd78af_0"
    shell:
        """
        cd {params.dir} &&
        RepeatMasker -pa {threads} -xsmall -gccalc -dir ./ -lib {params.repmo} {params.ref_upper} 2> {log} &&

        # Check if *.cat file is compressed or uncompressed
        if [ ! -f {output.rep_cat} ]
        then
          gzip {params.rep_cat_unzip}
        fi
        """


rule make_repeats_bed:
    """Generate a BED file from repeatmasker out file"""
    input:
        rep_out=rules.repeatmasker.output.rep_out,
    output:
        rep_bed=REF_DIR + "/" + REF_NAME + ".repeats.bed",
    group:
        "reference_group"
    log:
        os.path.abspath("results/logs/0.2_repeat_identification/" + REF_NAME + "_make_repeats_bed.log"),
    run:
        with open(input.rep_out, "r") as f, open(output.rep_bed, "w") as b:
            for line in f:
                edited = line.strip().split()
                if len(edited) > 6 and edited[0].isdigit():
                    start = str(int(edited[5]) - 1)
                    b.write(edited[4] + "\t" + start + "\t" + edited[6] + "\n")


rule sort_repeats_bed:
    input:
        rep_bed=rules.make_repeats_bed.output.rep_bed,
        genomefile=rules.genome_file.output.genomefile,
    output:
        sorted_rep_bed=REF_DIR + "/" + REF_NAME + ".repeats.sorted.bed",
    group:
        "reference_group"
    log:
        "results/logs/0.2_repeat_identification/" + REF_NAME + "_sort_repeats_bed.log",
    singularity:
        "community.wave.seqera.io/library/bedtools_htslib:62540682b559998a"
    shell:
        """
        bedtools sort -g {input.genomefile} -i {input.rep_bed} > {output.sorted_rep_bed} 2> {log}
        """


rule make_no_repeats_bed:
    input:
        ref_bed=rules.make_reference_bed.output,
        sorted_rep_bed=rules.sort_repeats_bed.output.sorted_rep_bed,
    output:
        no_rep_bed=REF_DIR + "/" + REF_NAME + ".repma.bed",
    group:
        "reference_group"
    log:
        "results/logs/0.2_repeat_identification/" + REF_NAME + "_make_no_repeats_bed.log",
    singularity:
        "community.wave.seqera.io/library/bedtools_htslib:62540682b559998a"
    shell:
        """
        bedtools subtract -a {input.ref_bed} -b {input.sorted_rep_bed} > {output.no_rep_bed} 2> {log}
        """
