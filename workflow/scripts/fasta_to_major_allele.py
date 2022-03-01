#!/usr/bin/python3

# Author: Tom van der Valk

from sys import argv
import operator

def fasta_to_major(filename, sample_to_exclude):

    outputfile = open(filename + ".parsed", "w")
    fasta_dict = {}
    counter = 0
    with open(filename) as f1:
        for line in f1:
            if line.startswith(">"):
                sample = line.strip()
                if sample == ">" + sample_to_exclude:
                    quiter = True
                else:
                    quiter = False
                    fasta_dict[sample] = ""
            else:
                if not quiter:
                    fasta_dict[sample] += line.strip()

    chrom_length = len(fasta_dict[list(fasta_dict.keys())[0]])
    for i in range(chrom_length):
        nucleotide = ""
        for value in list(fasta_dict.values()):
            if value[i] != "N":
                nucleotide += value[i]
        if len(nucleotide) < 3:
            outputfile.write(str(i + 1) + "\t" + "N" + "\n")
        else:
            A,C,T,G = nucleotide.count("A"), nucleotide.count("C"), nucleotide.count("T"), nucleotide.count("G")
            allel_dic = {}
            allel_dic["A"] = A
            allel_dic["C"] = C
            allel_dic["G"] = G
            allel_dic["T"] = T
            major_nucleotide = max(iter(allel_dic.items()), key=operator.itemgetter(1))[0]
            outputfile.write(str(i + 1) + "\t" + major_nucleotide + "\n")

    outputfile.close()

if __name__ == "__main__":
    fasta_to_major(argv[1], argv[2])

