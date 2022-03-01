#!/usr/bin/python3

"""
Script to combine mlRho results for all samples into one table.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import pandas as pd
import datetime

infiles=snakemake.input # input files from Snakemake rule
outtable=snakemake.output.table # output file from Snakemake rule

# define a function to get the input data
def mlRho_dataframe(inputfiles):
    dataframes = []
    for i in inputfiles: # create one dataframe per input file
        path_list = i.split("/")
        dataset = path_list[1]
        sample = path_list[4].split(".")[0]
        genomeregion = path_list[4].split(".")[-3]
        if genomeregion == "all":
            genomeregion = "genomewide"
        elif genomeregion == "autos":
            genomeregion = "autosomes"
        elif genomeregion == "sexchr":
            genomeregion = "sexchromosomes"
        sample_df=pd.read_table(i, delim_whitespace=True, comment='#')
        sample_df.insert(0, 'dataset', [dataset]) 
        sample_df.insert(1, 'sample', [sample])
        sample_df.insert(2, 'genomeregion', [genomeregion])
        dataframes.append(sample_df)
    conc_df = pd.concat(dataframes) # concatenate all dataframes
    conc_df.reset_index(drop=True, inplace=True)
    conc_df = conc_df.sort_values(by = ['dataset', 'sample'], ascending = [True, True])
    return conc_df

# apply the function
mlRho_df = mlRho_dataframe(infiles)

print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "mlRho results for each sample:")
print(mlRho_df)

mlRho_df.to_csv(outtable, sep='\t', na_rep='0.0', encoding='utf-8', index=False) # write to file