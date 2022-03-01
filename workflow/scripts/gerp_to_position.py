#!/usr/bin/python3

# Author: Tom van der Valk

from sys import argv
from itertools import islice
import operator
import gzip

def parse_gerp(fasta_file, gerp_file, species_name):

    seq = ""
    gerp_score = []
    outputfile = open(gerp_file + ".parsed", "w")

    adder = False
    with open(fasta_file) as f1:
        for line in f1:
            if line.startswith(">"):
                if line.startswith(">" + species_name):
                    adder = True
                else:
                    adder = False
            if adder and not line.startswith(">"):
                seq += line.strip()

    N_count = seq.count("N")

    with open(gerp_file) as f2:
        for line in f2:
            splitted = line.strip().split("\t")
            gerp_score += [splitted[1]]


    position = 0
    gerp_position = 0
    for i in (seq):
        position += 1
        if i.upper() == "N":
            outputfile.write("0" + "\n")
            gerp_position += 1
        else:
            outputfile.write(gerp_score[position - gerp_position -1] + "\n")

    outputfile.close()

if __name__ == "__main__":
    filename = argv[1]
    species_file = argv[2]
    species_name = argv[3]
    parse_gerp(filename, species_file, species_name)

