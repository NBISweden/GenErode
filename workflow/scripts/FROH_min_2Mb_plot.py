#!/usr/bin/python3

"""
Script to plot FROH, the proportion of the genome in runs of homozygosity (ROH) for ROHs >= 2MB, for historical and modern samples.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

intable=snakemake.input[0] # input table from Snakemake rule
outplot=snakemake.output.plot # output file from Snakemake rule

# define a function for the subplots
def froh_plot(dataframe, dataset, ax, color):
    ax.bar(dataframe['sample'][dataframe.dataset == dataset], dataframe['FROH'][dataframe.dataset == dataset], label=dataset, color=color)
    ax.set_xticks(range(len(dataframe['sample'][dataframe.dataset == dataset])))
    ax.set_xticklabels(list(dataframe['sample'][dataframe.dataset == dataset]), rotation = 45) # rotate x-axis tick labels

# define a function to find maximum value
def max_froh(dataframe):
    return dataframe['FROH'].max()

# plot the data
two_MB_df = pd.read_table(intable, sep="\t") # read in the FROH table

cividis = cm.get_cmap('cividis', 256) # create the cividis colormap

if len(two_MB_df['dataset'].unique())==2: # if there are both modern and historical samples, plot both
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True) # create the plot
    ax1.set_xlabel('historical', fontsize=12)
    ax2.set_xlabel('modern', fontsize=12)
    froh_plot(two_MB_df, 'historical', ax1, cividis.colors[230])
    froh_plot(two_MB_df, 'modern', ax2, cividis.colors[40])
    fig.supxlabel('samples', fontsize=12) # x axis label for both pplots combined
elif len(two_MB_df['dataset'].unique())==1: # if there is only one type of data, plot only this
    fig, ax = plt.subplots(nrows=1, ncols=1) # create the plot
    if two_MB_df['dataset'].str.contains("historical").any():
        ax.set_xlabel('historical samples', fontsize=12)
        froh_plot(two_MB_df, 'historical', ax, cividis.colors[230])
    elif two_MB_df['dataset'].str.contains("modern").any():
        ax.set_xlabel('modern samples', fontsize=12)
        froh_plot(two_MB_df, 'modern', ax, cividis.colors[40])

# fix the plot layout
max_ylim = max_froh(two_MB_df) + (0.05 * max_froh(two_MB_df))
if max_ylim != 0:
    plt.ylim(0, max_ylim)
else:
    plt.ylim(0, None)

fig.set_figheight(6) # fix figure height

if len(two_MB_df['sample'].unique()) > 4: # fix figure width
    widthscale = len(two_MB_df['sample'].unique()) * 0.5
else:
    widthscale = len(two_MB_df['sample'].unique()) * 1.5
fig.set_figwidth(widthscale)

fig.supylabel('$F_{ROH}$', fontsize=12) # y axis label for botpysubplots combined
fig.align_labels() # align axis labels

plt.tight_layout() # fix figure layout

# save the figure
fig.savefig(outplot, bbox_inches='tight', format='pdf') # save the figure