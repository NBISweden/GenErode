#!/usr/bin/python3

"""
Script to plot numbers of variants of different effects from snpEff analysis, for historical and modern samples.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

infile=snakemake.input[0] # input file from Snakemake rule
outplot=snakemake.output[0] # output file from Snakemake rule

# define a function for the subplots
def gerp_plot(dataframe, dataset, ax, color):
    ax.bar(dataframe['sample'][dataframe.dataset == dataset], dataframe['relative_mutational_load'][dataframe.dataset == dataset], label=dataset, color=color)
    ax.set_xticks(range(len(dataframe['sample'][dataframe.dataset == dataset])))
    ax.set_xticklabels(list(dataframe['sample'][dataframe.dataset == dataset]), rotation = 45) # rotate x-axis tick labels

# define a function to find maximum value
def max_load(dataframe):
    if dataframe['relative_mutational_load'].isnull().values.all():
        return 0
    else:
        return dataframe['relative_mutational_load'].max()

# create one dataframe from the GERP input file for plotting
data_df=pd.read_table(infile, sep='\t')

# plot the data
cividis = cm.get_cmap('cividis', 256) # create the cividis colormap

if len(data_df['dataset'].unique())==2: # create barplots for each modern and historical sample
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True) # create the plot
    ax1.set_xlabel('historical', fontsize=12)
    ax2.set_xlabel('modern', fontsize=12)
    gerp_plot(data_df, 'historical', ax1, cividis.colors[230])
    gerp_plot(data_df, 'modern', ax2, cividis.colors[40])
    fig.supxlabel('samples') # x axis label
elif len(data_df['dataset'].unique())==1: # if only historical or only modern samples are available,
    fig, ax = plt.subplots(nrows=1, ncols=1) # create the plot
    if data_df['dataset'].str.contains('historical').any(): # create barplots for each historical sample
        ax.set_xlabel('historical samples', fontsize=12)
        gerp_plot(data_df, 'historical', ax, cividis.colors[230])
    elif data_df['dataset'].str.contains('modern').any(): # create barplots for each modern sample
        ax.set_xlabel('modern samples', fontsize=12)
        gerp_plot(data_df, 'modern', ax, cividis.colors[40])

# fix the plot layout
max_ylim = max_load(data_df) + (0.05 * max_load(data_df))
if max_ylim != 0:
    plt.ylim(0, max_ylim)
else:
    plt.ylim(0, None)

fig.set_figheight(6) # fix figure height

if len(data_df['sample'].unique()) > 4: # fix figure width
    widthscale = len(data_df['sample'].unique()) * 0.5
else:
    widthscale = len(data_df['sample'].unique()) * 1.5
fig.set_figwidth(widthscale)

fig.supylabel('relative mutational load', fontsize=12) # y axis label
fig.align_labels() # align axis labels

plt.tight_layout(rect=[0.01, 0, 1, 1]) # fix figure layout, allowing for more space between supylabel and ylabel

fig.savefig(outplot, bbox_inches='tight', format='pdf') # save the figure