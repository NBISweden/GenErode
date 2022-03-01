#!/usr/bin/python3

"""
Script to create a histogram of depths from the output of samtools depth.

Input and output files refer to Snakemake directives.

Author: Verena Kutschera
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

dpstatsfile = snakemake.input.dp # file containing genome-wide average depth, minimum and maximum depth thresholds
dpsitesfile = snakemake.input.tmp # file containing the depth per site
histpdf = snakemake.output.pdf # output file for plotting

# read the file containing genome-wide average depth, minimum and maximum depth thresholds
with open(dpstatsfile, "r") as stats:
    line = stats.readline().split()
    minDP = float(line[1])
    maxDP = float(line[2])
    avg = float(line[0])
    if minDP < 3: # set minimum depth to 3 if lower
        minDP = 3

# read the file containing the depth per site in chunks of 10 million sites
dp_chunk = pd.read_csv(dpsitesfile, sep='\t', chunksize=10000000)

totalHistList = [0 for i in range(50)] # list containing empty histogram, will be filled when chunks are processed
bin_edges = list(np.arange(0,102,2)) # bin edges for histogram calculation and plotting

# each chunk is processed in df format
for chunk in dp_chunk:
    dpsites = chunk.iloc[:, 2] # extract column with depth
    hist, newbins = np.histogram(dpsites, bins=bin_edges) # calculate histogram in numpy
    newHistList = hist.tolist() # convert numpy array to list
    sumHistList = [sum(pair) for pair in zip(newHistList, totalHistList)] # sum up histograms
    totalHistList = sumHistList # replace empty histogram with the updated histogram numbers

# plot the histogram
cividis = cm.get_cmap('cividis', 3) # create the cividis colormap
plt.bar(bin_edges[:-1], totalHistList, width = 4, color = cividis.colors[1]) # bar plot
plt.xlim([0,maxDP+(maxDP/2)]) # x axis defined by maximum depth threshold
plt.axvline(minDP, color=cividis.colors[0], linewidth=2) # horizontal line with minimum depth threshold
plt.axvline(maxDP, color=cividis.colors[0], linewidth=2) # horizontal line with maximum depth threshold
plt.axvline(avg, color=cividis.colors[2], linewidth=2) # horizontal line with average genome-wide depth
plt.xlabel('Depth of coverage') # x axis
plt.ylabel('Number of sites') # y axis
plt.savefig(histpdf, bbox_inches='tight', format='pdf')
