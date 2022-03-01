#!/usr/bin/python3

"""
Script to calculate relative mutational load per sample.
Implemented equation (see von Seth et al. 2021, https://doi.org/10.1038/s41467-021-22386-8): 

for all sites with GERP > minimum GERP score and < maximum GERP score (specified on command line, see below),
(sum of GERP scores, multiplied with 2 at homozygous sites) / (sum of derived alleles, homozygous sites counted as 2)

Author: Verena Kutschera

Usage:
python3 gerp_rel_mut_load_sample.py gerp_file min_gerp_score max_gerp_score out_file
"""

from sys import argv
import pandas as pd
import datetime

infile=argv[1] # input file with ancestral state, GERP score, and numbers of derived alleles for each individual per site
min_gerp_score=argv[2] # minimum GERP score for a site to be considered for calculations
max_gerp_score=argv[3] # maximum GERP score for a site to be considered for calculations
outtable=argv[4] # output file

# Function to fill a dictionary with the required values, summing them up per chunk
def chunk_mut_load_dict(chunk_df, min_gerp, max_gerp):
    mut_load_dict = {}
    for chunk in chunk_df:
        # prepare data to calculate the relative mutational load per sample: 
        chunk['gerp_score'] = pd.to_numeric(chunk['gerp_score'],errors='coerce')
        chunk.dropna(subset=['gerp_score'], inplace=True)
        gerp_filter = chunk.loc[chunk['ancestral_state'] != 'N'].loc[chunk['gerp_score'] > float(min_gerp)].loc[chunk['gerp_score'] < float(max_gerp)] # sites with GERP larger than minimum GERP and smaller than maximum GERP
        #derived_filter = chunk.loc[chunk['ancestral_state'] != 'N'] # uncomment to get total number of derived alleles per sample instead of only number of derived alleles per sample at GERP sites
        samplelist = gerp_filter.columns[4:]
        for sample in samplelist:
            gerp_filter[sample] = pd.to_numeric(gerp_filter[sample],errors='coerce')
            #derived_filter[sample] = pd.to_numeric(derived_filter[sample],errors='coerce') # uncomment to get total number of derived alleles per sample instead of only number of derived alleles per sample at GERP sites
            sampleID = sample.rstrip('_no_derived_alleles')
            gerp_filter[sampleID + '_gerp_score_corrected'] = gerp_filter[sample] * gerp_filter['gerp_score'] # correct the gerp score for each sample, for homozygous sites *2 and for heterozygous sites *1
            if sampleID not in mut_load_dict:
                mut_load_dict[sampleID] = [0.0, 0.0]
            sum_filtered_gerp_score = gerp_filter.loc[gerp_filter[sample] > 0][sampleID + '_gerp_score_corrected'].sum() # sum of corrected gerp scores per chromosome for sites with derived alleles
            sum_filtered_derived_alleles = gerp_filter[sample].sum() # sum of number of derived alleles per chromosome with ancestral state not N and GERP > min_gerp (homozygous sites counted as 2, heterozygous sites counted as 1)
            #sum_filtered_derived_alleles = derived_filter[sample].sum() # uncomment to get the sum of number of derived alleles per chromosome with ancestral state not N
            mut_load_dict[sampleID][0] += sum_filtered_gerp_score
            mut_load_dict[sampleID][1] += sum_filtered_derived_alleles
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "GERP scores and derived alleles successfully summed up:")
    print(mut_load_dict)
    return mut_load_dict

# Function to calculate the relative mutational load per sample
def mut_load_sample(mut_load_dict):
    mut_load_list = []
    for sampleID in mut_load_dict.keys():
        mut_load = (mut_load_dict[sampleID][0] / mut_load_dict[sampleID][1]) # relative mutational load: sum of corrected GERP scores > min GERP and < max GERP/ sum of derived alleles at these sites
        sample_tup = (sampleID, mut_load)
        mut_load_list.append(sample_tup)
    mut_load_df = pd.DataFrame(mut_load_list, columns=['sample', 'relative_mutational_load'])
    print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Relative mutational load successfully calculated:")
    print(mut_load_df)
    return mut_load_df

# read the input file in chunks of 1 million sites
gerp_chunk = pd.read_csv(infile, sep='\t', chunksize=100000, dtype = {'#CHROM': str, 'POS': int, 'ancestral_state': str, 'gerp_score': float})

# apply the functions to the data
mut_load_dict = chunk_mut_load_dict(gerp_chunk, min_gerp_score, max_gerp_score)
mut_load_df = mut_load_sample(mut_load_dict)

# Write the output to file
mut_load_df.to_csv(outtable, sep='\t', na_rep='NaN', encoding='utf-8', index=False)
