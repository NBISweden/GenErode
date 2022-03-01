#!/usr/bin/python3

"""
Script to extract numbers of variants of different impact categories from snpEff analysis, for historical and modern samples.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import pandas as pd
import datetime

infiles=snakemake.input # input files from Snakemake rule. Path has to be structured like so for the script to work: "results/[dataset]/snpEff/[reference]/[sample]*_stats.csv"
outtable=snakemake.output[0] # output file from Snakemake rule

# define a function to get the input data
def impact_dataframe(inputfiles):
    # initialize a dictionary with lists containing placeholder values for variant impact for each sample
    counts_dict = {"dataset": ['NA' for i in range(len(inputfiles))], "sample": ['NA' for i in range(len(inputfiles))],
                "HIGH": [0 for i in range(len(inputfiles))], "LOW": [0 for i in range(len(inputfiles))], 
                "MODERATE": [0 for i in range(len(inputfiles))], "MODIFIER": [0 for i in range(len(inputfiles))]}
    # fill the dictionary with values for each sample
    count = 0
    for i in inputfiles:
        path_list = i.strip().split("/")
        dataset = path_list[1]
        sample = path_list[4].split(".")[0]
        counts_dict["dataset"][count] = dataset # add the dataset to the dict
        counts_dict["sample"][count] = sample # add the sample name to the dict
        with open(i, 'r') as f:
            for line in f:
                line_list = line.strip().split(" , ")
                if line.startswith("HIGH"): # get the counts of high impact variants per sample
                    counts_dict[line_list[0]][count] = int(line_list[1])
                elif line.startswith("LOW"): # get the counts of low impact variants per sample
                    counts_dict[line_list[0]][count] = int(line_list[1])
                elif line.startswith("MODERATE"): # get the counts of moderate impact variants per sample
                    counts_dict[line_list[0]][count] = int(line_list[1])
                elif line.startswith("MODIFIER"): # get the counts of modifier impact variants per sample
                    counts_dict[line_list[0]][count] = int(line_list[1])
                else:
                    continue
        count += 1
    counts_df = pd.DataFrame(counts_dict) # convert the dictionary to a pandas DataFrame
    counts_df = counts_df.sort_values(by = ['dataset', 'sample'], ascending = [True, True]) # sort by datatype and sample ID
    counts_df.rename(columns = {'HIGH':'high', 'LOW':'low', 'MODERATE':'moderate', 'MODIFIER':'modifier'}, inplace = True)
    return counts_df

impact_df = impact_dataframe(infiles) # read in the data

print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Counts of SNPs from different impact categories for each sample:")
print(impact_df)

impact_df.to_csv(outtable, sep='\t', na_rep='0.0', encoding='utf-8', index=False) # write to file
