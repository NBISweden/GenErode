#!/usr/bin/python3

"""
Script to calculate FROH, the proportion of the genome in runs of homozygosity (ROH) for ROHs >= 2MB, for historical and modern samples.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import pandas as pd
import matplotlib.pyplot as plt
import datetime

ingenfil=snakemake.input.genomefile # input genomefile from Snakemake rule
inroh=snakemake.input.ROH # input ROH files from Snakemake rule: per dataset, *.hom and *.hom.indiv
outtable=snakemake.output.table # output file from Snakemake rule

# read in the genome file and calculate the genome length
genfil=pd.read_table(ingenfil, header=None)
genlen=genfil.iloc[:,1].sum()

# create one dataframe from all ROH input files for plotting
if len(inroh)==4:
    dataset1=str(inroh[0].split("/")[1])
    dataset1_df=pd.read_table(inroh[0], delim_whitespace=True)
    if len(dataset1_df) == 0: # for samples without any ROH, add a row with data from the *.hom.indiv file
        dataset1_ind_df=pd.read_table(inroh[1], delim_whitespace=True)
        dataset1_df = pd.merge(dataset1_df, dataset1_ind_df, how="right", on=['FID', 'IID', 'PHE', 'KB'])
    dataset1_df['DATASET']=dataset1
    dataset2=str(inroh[2].split("/")[1])
    dataset2_df=pd.read_table(inroh[2], delim_whitespace=True)
    if len(dataset2_df) == 0: # for samples without any ROH, add a row with data from the *.hom.indiv file
        dataset2_ind_df=pd.read_table(inroh[3], delim_whitespace=True)
        dataset2_df = pd.merge(dataset2_df, dataset2_ind_df, how="right", on=['FID', 'IID', 'PHE', 'KB'])
    dataset2_df['DATASET']=dataset2   
    data_df = pd.concat([dataset1_df, dataset2_df], ignore_index=True)
elif len(inroh)==2:
    dataset=str(inroh[0].split("/")[1])
    data_df=pd.read_table(inroh[0], delim_whitespace=True)
    if len(data_df) == 0: # for samples without any ROH, add a row with data from the *.hom.indiv file
        data_ind_df=pd.read_table(inroh[1], delim_whitespace=True)
        data_df = pd.merge(data_df, data_ind_df, how="right", on=['FID', 'IID', 'PHE', 'KB'])
    data_df['DATASET']=dataset

data_df.KB = pd.to_numeric(data_df.KB, errors='coerce') # make sure the KB column is numeric
data_df['BP'] = data_df['KB'].map(lambda x: x * 1000) # add a column with length of ROHs in base pairs (BP)
samplelen1=data_df['BP'][data_df.BP >= 2000000.0].groupby([data_df['DATASET'], data_df['FID']]).sum() # sum up all ROH length (in bp) for each sample for ROH >= 2 Mb, incl. DATASET info
samplelen2=data_df['BP'][data_df.BP < 2000000.0].groupby([data_df['DATASET'], data_df['FID']]).sum() # sum up all ROH length (in bp) for each sample for ROH < 2 Mb , incl. DATASET info
samplelen2[samplelen2 >= 0] = 0.0 # set ROH length for samples with ROH < 2 Mb to zero
samplelen = samplelen1.add(samplelen2, fill_value=0) # sum up the two series so that all samples are represented, filling missing data with zero

two_MB = samplelen.div(genlen) # calculate proportion of ROHs, FROH, for ROH >= 2Mb for each sample (samples with ROH < 2 Mb are shown as zero)
two_MB.name = "FROH"
two_MB_df = two_MB.to_frame().reset_index() # convert to dataframe
two_MB_df = two_MB_df.sort_values(by = ['DATASET', 'FID'], ascending = [True, True]) # sort the dataframe
two_MB_df.rename(columns = {'DATASET':'dataset', 'FID':'sample'}, inplace = True)

print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Proportion of ROHs, FROH, for ROH >= 2Mb for each sample (samples with ROH < 2 Mb are shown as zero):")
print(two_MB_df)

# Write dataframe to output table
two_MB_df.to_csv(outtable, sep='\t', na_rep='0.0', encoding='utf-8', index=False)