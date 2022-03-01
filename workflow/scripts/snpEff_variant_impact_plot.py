#!/usr/bin/python3

"""
Script to plot numbers of variants from different impact categories from snpEff analysis, for historical and modern samples.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

infile=snakemake.input[0] # input files from Snakemake rule. Path has to be structured like so for the script to work: "results/[dataset]/snpEff/[reference]/[sample]*_stats.csv"
outplot=snakemake.output[0] # output file from Snakemake rule


# define a function for the subplots
def impact_plot(dataframe, dataset, impact, ax, color):
    ax.bar(dataframe['sample'][dataframe.dataset == dataset], dataframe[impact][dataframe.dataset == dataset], label=dataset, color=color)
    ax.set_xticks(range(len(dataframe['sample'][dataframe.dataset == dataset])))
    ax.set_xticklabels(list(dataframe['sample'][dataframe.dataset == dataset]), rotation = 45)

# define a function to find maximum value
def max_impact(dataframe, impact):
    return dataframe[impact].max()

# read in the data
impact_df = pd.read_table(infile, sep="\t")

# create a plot with four subplots, one for each variant impact
cividis = cm.get_cmap('cividis', 256) # create the cividis colormap
impact = ['high', 'low', 'moderate', 'modifier']

if len(impact_df['dataset'].unique())==2: # create barplots for each modern and historical sample
    fig, axes = plt.subplots(nrows=4, ncols=2, sharex='col', sharey='row')
    for index, row in enumerate(axes):
        ylabel = "{} impact variants".format(impact[index])
        row[0].set_ylabel(ylabel, fontsize = 12)
        impact_plot(impact_df, 'historical', impact[index], row[0], cividis.colors[230])
        impact_plot(impact_df, 'modern', impact[index], row[1], cividis.colors[40])
        max_ylim = max_impact(impact_df, impact[index]) + (0.05 * max_impact(impact_df, impact[index]))
        if max_ylim != 0:
            row[0].set_ylim(0,max_ylim)
            row[1].set_ylim(0,max_ylim)
        else:
            row[0].set_ylim(0,None)
            row[1].set_ylim(0,None)
    axes[3][0].set_xlabel('historical', fontsize=12)
    axes[3][1].set_xlabel('modern', fontsize=12)
    fig.supxlabel('samples') # x axis label

elif len(impact_df['dataset'].unique())==1: # if only historical or only modern samples are available,
    fig, axes = plt.subplots(nrows=4, ncols=1, sharex='col')
    for index, row in enumerate(axes):
        ylabel = "{} impact variants".format(impact[index])
        row.set_ylabel(ylabel.lower(), fontsize = 12)
        if impact_df['dataset'].str.contains("historical").any(): # create barplots for each historical sample
            impact_plot(impact_df, 'historical', impact[index], row, cividis.colors[230])
            axes[3].set_xlabel('historical samples', fontsize=12)
        elif impact_df['dataset'].str.contains("modern").any(): # create barplots for each modern sample
            impact_plot(impact_df, 'modern', impact[index], row, cividis.colors[40])
            axes[3].set_xlabel('modern samples', fontsize=12)
        max_ylim = max_impact(impact_df, impact[index]) + (0.05 * max_impact(impact_df, impact[index]))
        if max_ylim != 0:
            row.set_ylim(0,max_ylim)
        else:
            row.set_ylim(0,None)

# fix the plot layout
fig.set_figheight(12) # fix figure height

if len(impact_df['sample'].unique()) > 4: # fix figure width
    widthscale = len(impact_df['sample'].unique()) * 0.8
else:
    widthscale = len(impact_df['sample'].unique()) * 1.5
fig.set_figwidth(widthscale)

fig.supylabel('number of') # common y axis label
fig.align_labels() # align axis labels

plt.tight_layout() # fix figure layout

fig.savefig(outplot, bbox_inches='tight', format='pdf') # save the figure