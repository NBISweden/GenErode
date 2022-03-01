#!/usr/bin/env python3

'''
This script takes a gzipped VCF file containing variant and monomorphic/invariant sites from one single individual, 
and finds the positions of CpG sites within the resequencing data. Not only REF and ALT, but also genotype information 
(individual GT) are considered. Insertion/deletions (indels) and multiallelic sites are considered, too.
Site pairs with positions not directly adjacent to each other will be ignored. The output is printed to a bed file.

Author: Verena Kutschera

Usage:
python find_CpG_genotypes.py resequencedIndividual.vcf.gz resequencedIndividual_CpG.bed
'''

import gzip
import sys
import re
from collections import deque

vcfRead = gzip.open(sys.argv[1], 'rt') #open a gzipped VCF file
bedWrite = open(sys.argv[2], 'w') #write to outputfile in bed format (see https://genome.ucsc.edu/FAQ/FAQformat.html#format1)

### a set of functions to be called in the for loop through the VCF file
# splits each line and returns a list of arguments taken by the other functions
def getSiteInfo(site):
	col = site.strip().split()
	chrom = col[0] # the scaffold name
	pos = col[1] # the position
	ref = col[3] # the reference base
	alt = col[4] # the alt base
	info = col[7].strip().split(";") # all info for the site in one list
	geno = [] # all genotypes in a list to be checked by the functions
	for gt in col[9:]:
		geno.append(gt.split(':')[0])
	return chrom, pos, ref, alt, info, geno

# finds sites with "C" in the genotypes, and return scaffold name, position and fixed allele
def Cvariant(chrom, pos, ref, alt, info, geno):
	if ref=="C":
		if any("0/0" in g for g in geno) or any("0/1" in g for g in geno) or any("0/2" in g for g in geno): # any genotype with a "C"
			return chrom, pos
	elif alt=="C": # capturing sites that have "C" as the alt allele but have a different allele as ref
		if any("0/1" in g for g in geno) or any("1/1" in g for g in geno):
			return chrom, pos
	elif re.search(r"^C,[AGCTN]", alt): # capturing multiallelic sites including a "C"
		if any("0/1" in g for g in geno) or any("1/1" in g for g in geno) or any("1/2" in g for g in geno):
			return chrom, pos
	elif re.search(r"[AGCTN],C$", alt): # capturing multiallelic sites including a "C"
		if any("0/2" in g for g in geno) or any("1/2" in g for g in geno) or any("2/2" in g for g in geno):
			return chrom, pos
	elif re.search(r"[AGCTN]+C$", ref): # reference insertion ending with "C"
		if any("0/0" in g for g in geno) or any("0/1" in g for g in geno) or any("0/2" in g for g in geno): # any genotype with the reference insertion ending with "C"
			return chrom, pos
	elif re.search(r"[AGCTN]+C$", alt): # alt insertion ending with "C"
		if any("0/1" in g for g in geno) or any("1/1" in g for g in geno): # any genotype with the alt insertion ending with "C"
			return chrom, pos
	elif re.search(r"[AGCTN]*C,[AGCTN]+", alt): # capturing multiallelic sites including an insertion ending with "C"
		if any("0/1" in g for g in geno) or any("1/1" in g for g in geno) or any("1/2" in g for g in geno):
			return chrom, pos
	elif re.search(r"[AGCTN]+,[AGCTN]*C$", alt): # capturing multiallelic sites including an insertion ending with "C"
		if any("0/2" in g for g in geno) or any("1/2" in g for g in geno) or any("2/2" in g for g in geno):
			return chrom, pos


# finds sites with "G" in the genotypes, and return scaffold name, position and fixed allele
def Gvariant(chrom, pos, ref, alt, info, geno):
	if ref=="G":
		if any("0/0" in g for g in geno) or any("0/1" in g for g in geno) or any("0/2" in g for g in geno): # any genotype with a "G"
			return chrom, pos
	elif alt=="G": # capturing sites that have "G" as the alt allele but have a different allele as ref
		if any("0/1" in g for g in geno) or any("1/1" in g for g in geno):
			return chrom, pos
	elif re.search(r"^G,[AGCTN]", alt): # capturing multiallelic sites including a "G"
		if any("0/1" in g for g in geno) or  any("1/1" in g for g in geno) or any("1/2" in g for g in geno):
			return chrom, pos
	elif re.search(r"[AGCTN],G$", alt): # capturing multiallelic sites including a "G"
		if any("0/2" in g for g in geno) or any("1/2" in g for g in geno) or any("2/2" in g for g in geno):
			return chrom, pos
	elif re.search(r"^G[AGCTN]+", ref): # reference insertion starting with "G"
		if any("0/0" in g for g in geno) or any("0/1" in g for g in geno) or any("0/2" in g for g in geno): # any genotype with the reference insertion starting with "G"
			return chrom, pos
	elif re.search(r"^G[AGCTN]+", alt): # alt insertion starting with "G"
		if any("0/1" in g for g in geno) or any("1/1" in g for g in geno): # any genotype with the alt insertion starting with "G"
			return chrom, pos
	elif re.search(r"^G[AGCTN]*,[AGCTN]+", alt): # capturing multiallelic sites including an insertion starting with "G"
		if any("0/1" in g for g in geno) or any("1/1" in g for g in geno) or any("1/2" in g for g in geno):
			return chrom, pos
	elif re.search(r"[AGCTN]+,G[AGCTN]*", alt): # capturing multiallelic sites including an insertion starting with "G"
		if any("0/2" in g for g in geno) or any("1/2" in g for g in geno) or any("2/2" in g for g in geno):
			return chrom, pos


# converts line from VCF into bed format
def bedFormat(chrom, pos):
	startpos = int(pos)-1
	return str(chrom) + "\t" + str(startpos) + "\t" + str(pos)
### end of functions


# Loop through VCF file, store subsequent lines in a deque to be checked by the functions
lineDeque = deque(maxlen=2) # store the deque of a maximum length of two (sliding window of two sites)
for line in vcfRead:
	if line.startswith('#'): # skip the header section
		continue
	else:
		lineDeque.append(line) # move the sliding window one line forward
		if len(lineDeque)==2: # only continue if there are at least two lines in the deque
			siteInfo = getSiteInfo(lineDeque[0]) # split the first line in the deque into its parts to get the relevant info needed for the other functions
			nextSiteInfo = getSiteInfo(lineDeque[1]) # split the next site so that it's parameters can be taken by the functions
			if nextSiteInfo[0]==siteInfo[0] and (int(nextSiteInfo[1])-1==int(siteInfo[1])): # check that the sites are on the same scaffold and consecutive
				Csite = Cvariant(*siteInfo) # run the function to test if the resequenced individual has a "C" at that site. The arguments are taken from the list siteInfo.
				if Csite!=None: # if the site contains a "C"
					GnextSite = Gvariant(*nextSiteInfo) # run the function to check if the next site contains a "G" based on the list nextSiteInfo.
					if GnextSite!=None: # check if the next site contains a "G"
						print(bedFormat(*Csite) + "\n" + bedFormat(*GnextSite), file=bedWrite) # print the sites in bed format

vcfRead.close()
bedWrite.close()