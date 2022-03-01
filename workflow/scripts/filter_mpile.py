#!/usr/bin/python3

# Author: Tom van der Valk

import sys
import operator

for line in sys.stdin:
    splitted = line.strip().split("\t")
    nucleotide = splitted[4].upper()
    A,C,T,G = nucleotide.count("A"), nucleotide.count("C"), nucleotide.count("T"), nucleotide.count("G")
    if A+C+T+G == 1 and len(nucleotide) == 1:
        print (nucleotide)
    elif A+C+T+G == 0:
        print ("N")
    else:
        print ("N")

