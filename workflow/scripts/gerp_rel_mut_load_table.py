#!/usr/bin/python3

"""
Script to combine relative mutational load results for all samples into one table.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import pandas as pd
import datetime

infiles=snakemake.input # input files from Snakemake rule
outtable=snakemake.output[0] # output file from Snakemake rule

# define a function to get the input data
def rel_load_dataframe(inputfiles):
    dataframes = []
    for i in inputfiles: # create one dataframe per input file
        path_list = i.split("/")
        dataset = path_list[2]
        sample_df=pd.read_table(i, delim_whitespace=True)
        sample_df.insert(0, 'dataset', [dataset])
        dataframes.append(sample_df)
    conc_df = pd.concat(dataframes) # concatenate all dataframes
    conc_df.reset_index(drop=True, inplace=True)
    conc_df = conc_df.sort_values(by = ['dataset', 'sample'], ascending = [True, True])
    return conc_df

# apply the function
rel_load_df = rel_load_dataframe(infiles)

print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "Relative mutational load results for each sample:")
print(rel_load_df)

rel_load_df.to_csv(outtable, sep='\t', na_rep='0.0', encoding='utf-8', index=False) # write to file