##########################################################################
### 0.2 Repeat prediction and repeat masking of the reference genome

# Code collecting output files from this part of the pipeline
repeat_id_outputs=[]
repeat_id_outputs.append("results/references/" + REF_NAME + "/" + REF_NAME + ".repeats.sorted.bed")
repeat_id_outputs.append("results/references/" + REF_NAME + "/" + REF_NAME + ".repma.bed")


# snakemake rules
rule ref_upper:
    """Reference assembly file preparation: change all lowercase bases to uppercase"""
    input:
        ref=config["ref_path"],
    output:
        ref_upper="results/references/" + REF_NAME + "/" + REF_NAME + ".upper.fasta",
    log:
        "results/logs/0.2_repeat_identification/" + REF_NAME + "_ref_upper.log",
    shell:
        """
        awk '{{ if ($0 !~ />/) {{print toupper($0)}} else {{print $0}} }}' {input.ref} > {output.ref_upper} 2> {log}
        """


if os.path.exists(config["repeat_bed_file"]):
    rule sort_userprovided_repeat_bed:
        """Sort userprovided repeat BED file"""
        input:
            rep_bed=config["repeat_bed_file"],
            genomefile=rules.genome_file.output.genomefile,
        output:
            sorted_rep_bed="results/references/" + REF_NAME + "/" + REF_NAME + ".repeats.sorted.bed",
        log:
            "results/logs/0.2_repeat_identification/" + REF_NAME + "_sort_userprovided_repeat_bed.log",
        singularity:
            bedtools_htslib_container
        shell:
            """
            bedtools sort -g {input.genomefile} -i {input.rep_bed} > {output.sorted_rep_bed} 2> {log}
            """

else:
    rule repeatmodeler_builddatabase:
    """Build database for RepeatModeler"""
        input:
            ref_upper=rules.ref_upper.output,
        output:
            nhr=temp("results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + ".nhr"),
            nin=temp("results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + ".nin"),
            nnd=temp("results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + ".nnd"),
            nni=temp("results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + ".nni"),
            nog=temp("results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + ".nog"),
            nsq=temp("results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + ".nsq"),
            translation=temp("results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + ".translation"),
        container:
            repeatmodeler_container
        params:
            db="results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME,
        log:
            os.path.abspath("results/logs/0.2_repeat_identification/" + REF_NAME + "_repeatmodeler_builddatabase.log"),
        threads: 16
        shadow:
            "minimal"
        singularity:
            repeatmodeler_container
        shell:
            """
            BuildDatabase -engine ncbi -name {params.db} {input.ref_upper} &> {log}
            """


    rule repeatmodeler:
        """RepeatModeler for de novo repeat prediction from a reference assembly"""
        input:
            database=rules.repeatmodeler_builddatabase.output,
        output:
            repmo="results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + "-families.fa",
            stk="results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + "-families.stk",
            log="results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME + "-rmod.log",
        log:
            os.path.abspath("results/logs/0.2_repeat_identification/" + REF_NAME + "_repeatmodeler.log"),
        container:
            repeatmodmask_container
        params:
            db="results/references/" + REF_NAME + "/repeatmodeler/" + REF_NAME,
        threads: 16
        shadow:
            "minimal"
        singularity:
            repeatmodeler_container
        shell:
            """
            RepeatModeler -engine ncbi -threads {threads} -database {params.db} -quick &> {log}
            """


    rule repeatmasker:
        """Repeat mask the full genome assembly using raw de novo predicted repeats"""
        input:
            ref_upper=rules.ref_upper.output,
            repmo=rules.repeatmodeler.output.repmo,
        output:
            rep_tbl="results/references/" + REF_NAME + "/repeatmasker/" + REF_NAME + ".upper.fasta.tbl",
            rep_out="results/references/" + REF_NAME + "/repeatmasker/" + REF_NAME + ".upper.fasta.out",
        log:
            os.path.abspath("results/logs/0.2_repeat_identification/" + REF_NAME + "_repeatmasker.log"),
        params:
            dir="results/references/" + REF_NAME + "/repeatmasker/",
        threads: 16
        shadow:
            "minimal"
        singularity:
            repeatmodeler_container
        shell:
            """
            RepeatMasker -pa {threads} -xsmall -gccalc -dir {params.dir} -lib {input.repmo} {input.ref_upper} 2> {log}
            """


    rule make_repeats_bed:
        """Generate a BED file from repeatmasker out file"""
        input:
            rep_out=rules.repeatmasker.output.rep_out,
        output:
            rep_bed="results/references/" + REF_NAME + "/" + REF_NAME + ".repeats.bed",
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
            sorted_rep_bed="results/references/" + REF_NAME + "/" + REF_NAME + ".repeats.sorted.bed",
        log:
            "results/logs/0.2_repeat_identification/" + REF_NAME + "_sort_repeats_bed.log",
        singularity:
            bedtools_htslib_container
        shell:
            """
            bedtools sort -g {input.genomefile} -i {input.rep_bed} > {output.sorted_rep_bed} 2> {log}
            """


rule make_no_repeats_bed:
    input:
        ref_bed=rules.make_reference_bed.output,
        sorted_rep_bed="results/references/" + REF_NAME + "/" + REF_NAME + ".repeats.sorted.bed",
        genomefile=rules.genome_file.output.genomefile,
    output:
        no_rep_bed="results/references/" + REF_NAME + "/" + REF_NAME + ".repma.bed",
    log:
        "results/logs/0.2_repeat_identification/" + REF_NAME + "_make_no_repeats_bed.log",
    singularity:
        bedtools_htslib_container
    shell:
        """
        bedtools subtract -a {input.ref_bed} -b {input.sorted_rep_bed} \
        -sorted -g {input.genomefile} > {output.no_rep_bed} 2> {log}
        """
