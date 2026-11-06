#!/usr/bin/python3

"""
Script to generate a BED file from a RepeatMasker .out file.

Input and output files refer to Snakemake directives.
"""

infile = snakemake.input.rep_out  # input file from Snakemake rule
outbed = snakemake.output.rep_bed  # output file from Snakemake rule

with open(infile, "r") as f, open(outbed, "w") as b:
    for line in f:
        edited = line.strip().split()
        if len(edited) > 6 and edited[0].isdigit():
            start = str(int(edited[5]) - 1)
            b.write(edited[4] + "\t" + start + "\t" + edited[6] + "\n")