#!/usr/bin/python3

# Author: Tom van der Valk

from sys import argv
from itertools import islice
import operator
import gzip
from Bio import SeqIO

def make_fastq(fasta_file,outputfile):
    "open outputfile"
    outputfile = gzip.open(outputfile,"wt")
    "load fasta file into python dictionary"
    fasta_file = gzip.open(fasta_file,"rt")
    fasta_dict = {rec.id : rec.seq for rec in SeqIO.parse(fasta_file, "fasta")}
    "set fastq header name to 0"
    header = 0
    "loop through fasta dictionary in windows of 35bp"
    for key,value in fasta_dict.items():
        for i in range(0,len(value),35):
            fastq_read = str(value[i:i+35])
            "make sure no missing bases are in the sequence and length is correct"
            if fastq_read.count("N") == 0 and len(fastq_read) == 35:
                "assign header"
                header += 1
                if header % 100000 == 0:
                    print(header)
                outputfile.write("@_" + str(header) + "\n" + fastq_read + "\n" + "+" + "\n" + "J" * len(fastq_read) + "\n")
    outputfile.close()

if __name__ == "__main__":
    fasta_file = argv[1]
    outputfile = argv[2]
    make_fastq(fasta_file,outputfile)
