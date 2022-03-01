#!/usr/bin/python
from __future__ import division
from Bio import SeqIO
import operator
from sys import argv
import math

"""
Author: Tom van der Valk
tom.vandervalk@ebc.uu.se

Script modifications: Verena Kutschera

Script identifies CpG sites in a fasta reference
output is in BED-format
"""

def find_cpg(reference,outputname):

    outputfile = open(outputname + ".bed", "w")
    fasta_sequences = SeqIO.parse(open(reference),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq).upper()
        for i in range(len(sequence)-1):
            position = i
            if (sequence[i] + sequence[i+1]) == "CG":
                outputfile.write(name + "\t" + str(position) + "\t" + str(position + 2) + "\n")

    outputfile.close()

if __name__ == "__main__":
    find_cpg(argv[1],argv[2])
