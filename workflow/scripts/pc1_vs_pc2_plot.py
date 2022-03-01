#!/usr/bin/python3

"""
Script to plot the PCA axes PC1 and PC2.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

eigenvec=snakemake.input.eigenvec # input file from Snakemake rule
eigenval=snakemake.input.eigenval # input file from Snakemake rule
outplot=snakemake.output[0] # output file from Snakemake rule

if os.stat(eigenvec).st_size > 0:
    eig_vec = pd.read_table(eigenvec, delim_whitespace=True, header=None)
    eig_val = pd.read_table(eigenval, delim_whitespace=True, header=None) # read in to calculate percentage variance explained
    if len(eig_vec.columns) > 3:
        pva = eig_val[0] / eig_val[0].sum() * 100 # calculate percentage variance explained
        eig_vec.rename(columns={0: 'FID', 1: 'IID', 2: 'PC1', 3: 'PC2'}, inplace=True)
        grouped = eig_vec.groupby('FID')
        cividis = cm.get_cmap('cividis', len(grouped)) # create the cividis colormap
        fig,ax = plt.subplots()
        n=0 # counter to take different colors
        for k,d in grouped:
            ax.scatter(d['PC1'], d['PC2'], label=k, c=cividis.colors[n].reshape(1,-1))
            n+=1
        plt.legend(bbox_to_anchor=(1,1), loc="upper left")
        plt.xlabel('PC1' + ' (' + str(round(pva[0],1)) + ' %)')
        plt.ylabel('PC2' + ' (' + str(round(pva[1],1)) + ' %)')
else:
    fig,ax = plt.subplots()
    plt.xlabel('PC1')
    plt.ylabel('PC2')
plt.savefig(outplot, bbox_inches='tight', format='pdf')
