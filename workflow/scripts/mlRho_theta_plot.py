#!/usr/bin/python3

"""
Script to plot theta values estimated in mlRho, for historical and modern samples.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

intable=snakemake.input[0] # input file from Snakemake rule
outplot=snakemake.output[0] # output file from Snakemake rule

# define a function for the subplots
def theta_plot(dataframe, grouped, dataset, genomeregion, ax, color):
    data_df = grouped.get_group((dataset, genomeregion))
    for lower,theta,upper,x in zip(data_df['lower_theta'],data_df['theta_est'],data_df['upper_theta'],range(len(data_df))):
        ax.plot((x,x),(lower,upper),'-',color=color)
        ax.plot((x),(theta),'o',color=color)
    ax.set_xticks(range(len(data_df)))
    ax.set_xticklabels(list(data_df['sample']), rotation = 45)

# read in and process the data
table_df = pd.read_table(intable, sep="\t")
table_df['lower_theta'] = table_df['theta'].str.split('<').str[0].astype(float)
table_df['theta_est'] = table_df['theta'].str.split('<').str[1].astype(float)
table_df['upper_theta'] = table_df['theta'].str.split('<').str[2].astype(float)
grouped_df = table_df.groupby(['dataset', 'genomeregion'])

# create the plot, with different layout depending on which samples, datasets and genome regions are plotted
cividis = cm.get_cmap('cividis', 256) # create the cividis colormap

if len(table_df['dataset'].unique())==2:
    if len(table_df['genomeregion'].unique())==1:
        fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, sharey=True)
        ax1.set_xlabel("historical", fontsize=12)
        ax2.set_xlabel("modern", fontsize=12)
        if table_df['genomeregion'].unique() == "genomewide":
            theta_plot(table_df, grouped_df, "historical", "genomewide", ax1, cividis.colors[230])
            theta_plot(table_df, grouped_df, "modern", "genomewide", ax2, cividis.colors[40])
        elif table_df['genomeregion'].unique() == "autosomes":
            theta_plot(table_df, grouped_df, "historical", "autosomes", ax1, cividis.colors[230])
            theta_plot(table_df, grouped_df, "modern", "autosomes", ax2, cividis.colors[40])
    else:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2, sharey='row', sharex='col')
        ax3.set_xlabel("historical", fontsize=12)
        ax4.set_xlabel("modern", fontsize=12)
        ax1.set_ylabel("autosomes", fontsize = 12)
        ax3.set_ylabel("sexchromosomes", fontsize = 12)
        theta_plot(table_df, grouped_df, "historical", "autosomes", ax1, cividis.colors[230])
        theta_plot(table_df, grouped_df, "modern", "autosomes", ax2, cividis.colors[40])
        theta_plot(table_df, grouped_df, "historical", "sexchromosomes", ax3, cividis.colors[230])
        theta_plot(table_df, grouped_df, "modern", "sexchromosomes", ax4, cividis.colors[40])
    fig.supxlabel('samples') # x axis label

elif len(table_df['dataset'].unique())==1:
    if len(table_df['genomeregion'].unique())==1:
        fig, ax = plt.subplots(ncols=1, nrows=1)
        if table_df['dataset'].str.contains("historical").any():
            ax.set_xlabel('historical samples', fontsize = 12) # x axis label
            if table_df['genomeregion'].unique() == "genomewide":
                theta_plot(table_df, grouped_df, "historical", "genomewide", ax, cividis.colors[230])
            elif table_df['genomeregion'].unique() == "autosomes":
                theta_plot(table_df, grouped_df, "historical", "autosomes", ax, cividis.colors[230])
        elif table_df['dataset'].str.contains("modern").any():
            ax.set_xlabel('modern samples', fontsize = 12) # x axis label
            if table_df['genomeregion'].unique() == "genomewide":
                theta_plot(table_df, grouped_df, "modern", "genomewide", ax, cividis.colors[40])
            elif table_df['genomeregion'].unique() == "autosomes":
                theta_plot(table_df, grouped_df, "modern", "autosomes", ax, cividis.colors[40])
    else:
        fig, (ax1, ax2) = plt.subplots(ncols=1, nrows=2, sharex=True)
        ax1.set_ylabel('autosomes', fontsize = 12)
        ax2.set_ylabel('sexchromosomes', fontsize = 12)
        if table_df['dataset'].str.contains("historical").any():
            ax2.set_xlabel('historical samples', fontsize = 12) # x axis label
            theta_plot(table_df, grouped_df, "historical", "autosomes", ax1, cividis.colors[230])
            theta_plot(table_df, grouped_df, "historical", "sexchromosomes", ax2, cividis.colors[230])
        elif table_df['dataset'].str.contains("modern").any():
            ax2.set_xlabel('modern samples', fontsize = 12) # x axis label
            theta_plot(table_df, grouped_df, "modern", "autosomes", ax1, cividis.colors[40])
            theta_plot(table_df, grouped_df, "modern", "sexchromosomes", ax2, cividis.colors[40])

# fix figure width and height based on data shape
if len(table_df['sample'].unique()) > 4:
    widthscale = len(table_df['sample'].unique()) / 1.5
else:
    if len(table_df['dataset'].unique()) == 1:
        widthscale = len(table_df['sample'].unique()) * 1.5
    else:
        widthscale = len(table_df['sample'].unique()) * 3
fig.set_figwidth(widthscale)

heightscale = len(table_df['genomeregion'].unique()) * 4
fig.set_figheight(heightscale)

fig.supylabel('theta', fontsize=12) # common y axis label
fig.align_labels() # align axis labels

plt.tight_layout()

# save the figure as PDF
fig.savefig(outplot, bbox_inches='tight', format='pdf')