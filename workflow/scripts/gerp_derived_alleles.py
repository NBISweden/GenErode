#!/usr/bin/python3

"""
Script to add the number of derived alleles per sample to the gerp output file, run separately for modern and historical VCF files.
Homozygous sites are coded as 2, heterozygous sites as 1.
Written to be run per window, by providing contig (or scaffold/chromosome) name, start and end position (1-based) on the command line.

Author: Verena Kutschera

Usage:
python3 gerp_derived_alleles.py gerp_file vcf_file_gzipped contig start_position end_position out_file_name
"""

from sys import argv
import gzip
import pandas as pd
import numpy as np
import datetime
from collections import deque
import warnings
warnings.filterwarnings("ignore")

gerp=argv[1] # gerp output file with ancestral state and gerp score per site
vcf=argv[2] # VCF file to be merged
contig=argv[3] # contig to be processed
start=argv[4] # start position for window (1-based)
end=argv[5] # end position for window
outfile=argv[6] # output file path

# Function to read the GERP output file per window for further processing
def read_gerp_windows(gerpFile, chrom, start, end):
    skip_rows = 0 # counter to track number of lines to skip when reading from the GERP file
    n_rows = 0 # counter to track number of lines to read from GERP file
    lineDeque = deque(maxlen=2) # store the deque of a maximum length of two
    with open(gerpFile, 'r') as f:
        for line in f:
            lineDeque.append(line) # move the sliding window one line forward
            if len(lineDeque) == 2: # only continue if there are exactly two lines in the deque
                currentGerpChrom = lineDeque[0].strip().split('\t')[0] # get the chromosome name of the first line in the deque
                currentGerpPos = int(lineDeque[0].strip().split('\t')[1]) # get the position of the first line
                nextGerpChrom = lineDeque[1].strip().split('\t')[0] # get the chromosome name of the second line
                if currentGerpChrom != chrom:
                    skip_rows += 1
                elif currentGerpChrom == chrom and currentGerpPos < int(start):
                    skip_rows += 1
                elif currentGerpChrom == chrom and nextGerpChrom == chrom and currentGerpPos >= int(start) and currentGerpPos <= int(end):
                    n_rows += 1
                elif currentGerpChrom == chrom and nextGerpChrom != chrom and currentGerpPos >= int(start) and currentGerpPos <= int(end):
                    n_rows += 1
                    break
                elif currentGerpChrom == chrom and nextGerpChrom == chrom and currentGerpPos >= int(start) and currentGerpPos > int(end):
                    break
                elif currentGerpChrom == chrom and nextGerpChrom != chrom and currentGerpPos >= int(start) and currentGerpPos > int(end):
                    break
            if len(lineDeque) == 1: # handle the last line in the file
                lastGerpChrom = lineDeque[0].strip().split('\t')[0] # get the chromosome name of the line in the deque
                lastGerpPos = int(lineDeque[0].strip().split('\t')[1]) # get the position of the line
                if lastGerpChrom == chrom and lastGerpPos >= int(start) and lastGerpPos == int(end):
                    n_rows += 1
                    break
    if n_rows > 0:
        gerpDF = pd.read_csv(gerpFile, sep='\t', skiprows=skip_rows, nrows=n_rows, 
                                names=['#CHROM', 'POS', 'ancestral_state', 'gerp_score'], 
                                dtype={'#CHROM': 'string', 'POS': int, 'ancestral_state': 'string', 'gerp_score': float}).set_index(['#CHROM', 'POS'])
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Chromosome", chrom, "from position", start, "to position", end, "from GERP file read into memory")
    else: # Add a row with "NaN" in case the window is missing from the file
        gerpDF = pd.DataFrame(columns=['#CHROM', 'POS', 'ancestral_state', 'gerp_score']).set_index(['#CHROM', 'POS'])
        gerpDF.loc[(chrom,start),:] = (np.nan)
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Chromosome", chrom, "from position", start, "to position", end, "not present in GERP file, will be set to 'NaN'")
    #print(gerpDF.info(memory_usage="deep"))
    return gerpDF

# Function to read the VCF file per window for further processing
def read_vcf_windows(vcfFile, chrom, start, end):
    skip_rows = 0 # counter to track number of lines to skip when reading from the VCF
    n_rows = 0 # counter to track number of lines to read from VCF
    lineDeque = deque(maxlen=2) # store the deque of a maximum length of two
    with gzip.open(vcfFile, 'r') as f:
        for line in f:
            if line.decode('utf8').startswith('##'):
                skip_rows += 1
            elif line.decode('utf8').startswith('#CHROM'):
                header = line.decode('utf8').strip().split('\t')
                skip_rows += 1
            else:
                lineDeque.append(line) # move the sliding window one line forward
                if len(lineDeque) == 2: # only continue if there are exactly two lines in the deque
                    currentVcfChrom = lineDeque[0].decode('utf8').strip().split('\t')[0] # get the chromosome name of the first line in the deque
                    currentVcfPos = int(lineDeque[0].decode('utf8').strip().split('\t')[1]) # get the position of the first line 
                    nextVcfChrom = lineDeque[1].decode('utf8').strip().split('\t')[0] # get the chromosome name of the second line 
                    if currentVcfChrom != chrom:
                        skip_rows += 1
                    elif currentVcfChrom == chrom and currentVcfPos < int(start):
                        skip_rows += 1
                    elif currentVcfChrom == chrom and nextVcfChrom == chrom and currentVcfPos >= int(start) and currentVcfPos <= int(end):
                        n_rows += 1
                    elif currentVcfChrom == chrom and nextVcfChrom != chrom and currentVcfPos >= int(start) and currentVcfPos <= int(end):
                        n_rows += 1
                        break
                    elif currentVcfChrom == chrom and nextVcfChrom == chrom and currentVcfPos >= int(start) and currentVcfPos > int(end):
                        break
                    elif currentVcfChrom == chrom and nextVcfChrom != chrom and currentVcfPos >= int(start) and currentVcfPos > int(end):
                        break
                if len(lineDeque) == 1: # handle the last line in the file
                    lastVcfChrom = lineDeque[0].strip().split('\t')[0] # get the chromosome name of the line in the deque
                    lastVcfPos = int(lineDeque[0].strip().split('\t')[1]) # get the position of the line
                    if lastVcfChrom == chrom and lastVcfPos >= int(start) and lastVcfPos == int(end):
                        n_rows += 1
                        break
    usecols_list = [0,1,3,4] + [*range(9,len(header))]
    usecols_header = [header[i] for i in usecols_list]
    if n_rows > 0:
        vcfDF = pd.read_csv(vcfFile, sep='\t', skiprows=skip_rows, nrows=n_rows, usecols=usecols_list, 
                              names=usecols_header, dtype = 'string', converters = {'POS': int}).set_index(['#CHROM', 'POS'])
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Chromosome", chrom, "from position", start, "to position", end, "from VCF file read into memory")
    else: # Add a row with "NaN" in case the window is missing from the file
        vcfDF = pd.DataFrame(columns=usecols_header).set_index(['#CHROM', 'POS'])
        vcfDF.loc[(chrom,start),:] = (np.nan)
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Chromosome", chrom, "from position", start, "to position", end, "not present in VCF file, will be set to 'NaN'")
    #print(vcfDF.info(memory_usage="deep"))
    return vcfDF

# Join GERP and VCFs for each window
def join_DFs_per_window(chrom, start, end, gerpFile, vcfFile):
    df_list = [] # List to collect DataFrames that can be removed from memory
    gerpDF = read_gerp_windows(gerpFile, chrom, start, end)
    df_list.append(gerpDF)
    vcfDF = read_vcf_windows(vcfFile, chrom, start, end)
    # Fix the genotype column to only contain data until the first ":"
    for i in range(2, len(vcfDF.columns)):
        vcfDF.iloc[:, i] = vcfDF.iloc[:, i].str.split(":").str[0]
    df_list.append(vcfDF)
    # Merge the DataFrames
    joinedDF = gerpDF.join(vcfDF, how='left')
    #print(joinedDF.info(memory_usage="deep"))
    # Remove the two DataFrames in the list from memory
    del df_list
    return joinedDF

# Function to calculate the number of derived alleles per site
def num_derived_alleles(row, samplename=""):
    no_derived_alleles = 0
    if pd.isna(row[samplename]): # if genotype is NaN
        no_derived_alleles = np.nan
    elif row['ancestral_state'] == "N" or pd.isna(row['ancestral_state']): # if the ancestral state is unknown
        no_derived_alleles = np.nan
    elif row['ancestral_state'] != "N": # if ancestral state is not "N"
        # split the genotype and create a list with possible values: "." (missing), "0" (REF), "1" (ALT), "2" (ALT), etc.
        focal_genotype = row[samplename].split("/")
        # create a list of REF and ALT alleles with REF at index 0, e.g. "C C,G" becomes ["C","C","G"] or "T G" becomes ["T","G"]
        ref_alt_list = list(row['REF']) + row['ALT'].split(",")
        # check for missing genotypes
        if "." in focal_genotype:
            no_derived_alleles = np.nan
        else:
            # extract the bases from the list of REF + ALT based on the index
            focal_alleleA = ref_alt_list[int(focal_genotype[0])]
            focal_alleleB = ref_alt_list[int(focal_genotype[1])]
            # check if the genotype is derived
            if focal_alleleA != row['ancestral_state']:
                no_derived_alleles += 1
            if focal_alleleB != row['ancestral_state']:
                no_derived_alleles += 1
            no_derived_alleles = str(no_derived_alleles)
    return no_derived_alleles

# Create new columns with numbers of derived alleles per sample
def add_num_derived_alleles(joinedDF):
    samplelist = list(joinedDF.columns[4:])
    for sample in samplelist:
        joinedDF[sample + '_no_derived_alleles'] = joinedDF.apply(num_derived_alleles, axis=1, samplename=sample)
        joinedDF.drop([sample], axis=1, inplace=True)
    # Drop columns that are no longer needed
    joinedDF.drop(['REF', 'ALT'], axis=1, inplace=True)
    #print(joinedDF.info(memory_usage="deep"))
    return joinedDF


# Apply the functions to the data
joinedWindowDF=join_DFs_per_window(contig, start, end, gerp, vcf)
print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "GERP output file and VCF file successfully joined for chromosome", contig, "from position", start, "to position", end)

finalWindowDF=add_num_derived_alleles(joinedWindowDF)
print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "GERP output file and VCF file successfully processed for chromosome", contig, "from position", start, "to position", end)

# Write the window to file for further processing
finalWindowDF.to_csv(outfile, sep='\t', na_rep='NaN', encoding='utf-8')
print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Result for chromosome", contig, "from position", start, "to position", end, "successfully written to", outfile, "\n")