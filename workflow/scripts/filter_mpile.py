#!/usr/bin/python3

# Author: Tom van der Valk (modified by Diana Robledo-Ruiz and Verena Kutschera)

# Filters a mpileup file for further processing in the GERP pipeline, produced with 
# `samtools mpileup -aa -r contigname --no-output-ends sample.bam`

import sys

for line in sys.stdin:
    splitted = line.strip().split("\t")

    # Check if at least one sample is present in the mpileup
    if len(splitted) < 5:
        print("Error: Invalid input line:", line)
        sys.stderr.write("Error: Invalid input line:" + line + '\n')
        break  # Stop the loop if invalid line

    # Process column 5 containing the mapped read bases
    nucleotide = splitted[4].upper()
    A,C,T,G = nucleotide.count("A"), nucleotide.count("C"), nucleotide.count("T"), nucleotide.count("G")
    
    # Print the nucleotide if it is A, C, T, or G, and if only one read mapped at the position
    if A+C+T+G == 1 and len(nucleotide) == 1:
        print(nucleotide)

    # Remove any positions other than one of the four bases, e.g. deletions indicated by "*"
    elif A+C+T+G == 0:
        print("N")

    # Remove any other read bases, including those followed by indels (indicated by "+" and an integer and sequence, or "-" and an integer)
    else:
        print("N")

