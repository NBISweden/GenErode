#!/usr/bin/python3

# Author: Tom van der Valk

from sys import argv
from itertools import islice
import operator
import gzip

def parse_mpile(mpile_file, chromname, refname):

    outputfile = open(mpile_file[:-6] + ".fasta","w")
    seq = ""
    with open(mpile_file, "r") as f1:
        for line in f1:
            splitted = line.strip()
            seq += splitted

    #to_substract = len(str(chromname)) + 7
    #outputfile.write(">" + mpile_file[:-to_substract] + "\n" + seq + "\n")
    outputfile.write(">" + refname + "\n" + seq + "\n")
    outputfile.close()


if __name__ == "__main__":
    filename = argv[1]
    chromname = argv[2]
    refname = argv[3]
    parse_mpile(filename, chromname, refname)

