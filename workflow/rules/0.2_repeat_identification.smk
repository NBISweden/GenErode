##########################################################################
### 0.2 Repeat prediction and repeat masking of the reference genome

# Code collecting output files from this part of the pipeline
all_outputs.append(REF_DIR + "/" + REF_NAME + ".repeats.bed")
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


rule cp_repeatmasker_libs:
    """Copy RepeatMasker libraries from container"""
    output:
        art=temp("workflow/resources/RepeatMasker/Libraries/Artefacts.embl"),
        embl=temp("workflow/resources/RepeatMasker/Libraries/Dfam.embl"),
        hmm=temp("workflow/resources/RepeatMasker/Libraries/Dfam.hmm"),
        repann=temp("workflow/resources/RepeatMasker/Libraries/RepeatAnnotationData.pm"),
        phr=temp("workflow/resources/RepeatMasker/Libraries/RepeatPeps.lib.phr"),
        psq=temp("workflow/resources/RepeatMasker/Libraries/RepeatPeps.lib.psq"),
        lib=temp("workflow/resources/RepeatMasker/Libraries/RepeatPeps.lib"),
        pin=temp("workflow/resources/RepeatMasker/Libraries/RepeatPeps.lib.pin"),
        peprm=temp("workflow/resources/RepeatMasker/Libraries/RepeatPeps.readme"),
        meta=temp("workflow/resources/RepeatMasker/Libraries/RMRBMeta.embl"),
        rm=temp("workflow/resources/RepeatMasker/Libraries/README.meta"),
        tax=temp("workflow/resources/RepeatMasker/Libraries/taxonomy.dat"),
    log:
        "results/logs/0.2_repeat_identification/" + REF_NAME + "_cp_repeatmasker_libs.log",
    singularity:
        "docker://quay.io/biocontainers/repeatmodeler:2.0.1--pl526_0"
    shell:
        """
        cp /usr/local/share/RepeatMasker/Libraries/* workflow/resources/RepeatMasker/Libraries/ 2> {log}
        """


rule embl2fasta:
    """Convert Dfam embl to fasta format"""
    input:
        dfam_embl=rules.cp_repeatmasker_libs.output.embl,
    output:
        rm_lib=temp("workflow/resources/RepeatMasker/Libraries/RepeatMasker.lib"),
    log:
        "results/logs/0.2_repeat_identification/" + REF_NAME + "_embl2fasta.log",
    run:
        from Bio import SeqIO
        with open(input.dfam_embl, "rU") as input_handle, open(output.rm_lib, "w") as output_handle:
            sequences = SeqIO.parse(input_handle, "embl")
            count = SeqIO.write(sequences, output_handle, "fasta")
        print("Converted %i records" % count)


rule make_repma_blast_db:
    input:
        rm_lib=rules.embl2fasta.output.rm_lib,
    output:
        nhr=temp("workflow/resources/RepeatMasker/Libraries/RepeatMasker.lib.nhr"),
        nin=temp("workflow/resources/RepeatMasker/Libraries/RepeatMasker.lib.nin"),
        nsq=temp("workflow/resources/RepeatMasker/Libraries/RepeatMasker.lib.nsq"),
    params:
        dir="workflow/resources/RepeatMasker/Libraries/",
        rm_lib="RepeatMasker.lib",
    log:
        os.path.abspath("results/logs/0.2_repeat_identification/" + REF_NAME + "_make_repma_blast_db.log"),
    singularity:
        "docker://quay.io/biocontainers/repeatmodeler:2.0.1--pl526_0"
    shell:
        """
        cd {params.dir}
        makeblastdb -dbtype nucl -in {params.rm_lib} 2> {log}
        """


rule repeatmodeler:
    """RepeatModeler for de novo repeat prediction from a reference assembly"""
    input:
        ref_upper=rules.ref_upper.output,
    output:
        repmo=REF_DIR + "/repeatmodeler/" + REF_NAME + "/RM_raw.out/consensi.fa",
        stk=REF_DIR + "/repeatmodeler/" + REF_NAME + "/RM_raw.out/families.stk",
    params:
        dir=REF_DIR + "/repeatmodeler/" + REF_NAME + "/",
        name=REF_NAME,
        ref_upper="../../" + REF_NAME + ".upper.fasta",
        abs_tmp=os.path.abspath("tmpConsensi.fa"),
    log:
        os.path.abspath("results/logs/0.2_repeat_identification/" + REF_NAME + "_repeatmodeler.log"),
    threads: 16
    singularity:
        "docker://quay.io/biocontainers/repeatmodeler:2.0.1--pl526_0"
    shell:
        """
        cd {params.dir}

        # Build repeat database
        BuildDatabase -engine ncbi -name {params.name} {params.ref_upper} 2> {log} &&

        # Run RepeatModeler
        RepeatModeler -engine ncbi -pa {threads} -database {params.name} 2>> {log} &&

        # copy the output files to a new directory
        cp RM_*.*/consensi.fa RM_raw.out/ 2>> {log} &&
        cp RM_*.*/families.stk RM_raw.out/ 2>> {log}

        # remove temporary file
        if [ -f {params.abs_tmp} ]
        then
          rm {params.abs_tmp} 2>> {log}
        fi
        """


rule repeatclassifier:
    """Create final RepeatModeler output files"""
    input:
        repmo=rules.repeatmodeler.output.repmo,
        stk=rules.repeatmodeler.output.stk,
        rm_lib=rules.embl2fasta.output.rm_lib,
        rm_db=rules.make_repma_blast_db.output,
        rm_libs=rules.cp_repeatmasker_libs.output,
    output:
        repmo=REF_DIR + "/repeatmodeler/" + REF_NAME + "/RM_raw.out/consensi.fa.classified",
        stk=REF_DIR + "/repeatmodeler/" + REF_NAME + "/RM_raw.out/families-classified.stk",
    params:
        repma_dir="workflow/resources/RepeatMasker",
    log:
        "results/logs/0.2_repeat_identification/" + REF_NAME + "_repeatclassifier.log",
    threads: 2
    singularity:
        "docker://quay.io/biocontainers/repeatmodeler:2.0.1--pl526_0"
    shell:
        """
        RepeatClassifier -repeatmasker_dir {params.repma_dir} -consensi {input.repmo} -stockholm {input.stk} 2> {log}
        """


rule repeatmasker:
    """Repeat mask the full genome assembly using raw de novo predicted repeats"""
    input:
        ref_upper=rules.ref_upper.output,
        repmo=rules.repeatclassifier.output.repmo,
        rm_lib=rules.embl2fasta.output.rm_lib,
        rm_db=rules.make_repma_blast_db.output,
        rm_libs=rules.cp_repeatmasker_libs.output,
    output:
        rep_masked=REF_DIR + "/repeatmasker/" + REF_NAME + "/" + REF_NAME + ".upper.fasta.masked",
        rep_align=REF_DIR + "/repeatmasker/" + REF_NAME + "/" + REF_NAME + ".upper.fasta.align",
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
        "docker://quay.io/biocontainers/repeatmasker:4.0.9_p2--pl526_2"
    shell:
        """
        export REPEATMASKER_LIB_DIR=$PWD/workflow/resources/RepeatMasker/Libraries &&
        cd {params.dir} &&
        RepeatMasker -pa {threads} -a -xsmall -gccalc -dir ./ -lib {params.repmo} {params.ref_upper} 2> {log} &&

        # Check if *.cat file is compressed or uncompressed
        if [ ! -f {output.rep_cat} ] # check if there are at least 2 files for merging. If there is only one file, copy the sorted bam file.
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
        "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
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
        no_rep_bed_dir="results/" + REF_NAME + ".repma.bed",
    group:
        "reference_group"
    log:
        "results/logs/0.2_repeat_identification/" + REF_NAME + "_make_no_repeats_bed.log",
    singularity:
        "docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
    shell:
        """
        bedtools subtract -a {input.ref_bed} -b {input.sorted_rep_bed} > {output.no_rep_bed} 2> {log} &&
        cp {output.no_rep_bed} {output.no_rep_bed_dir} 2>> {log}
        """
