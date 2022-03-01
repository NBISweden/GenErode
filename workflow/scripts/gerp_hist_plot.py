#!/usr/bin/python3

"""
Script to plot a histogram of GERP scores.

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

gerpfile = snakemake.input.gerp_out # file containing gerp scores along with ancestral state
histpdf = snakemake.output.pdf # output file for plotting

# read the file in chunks of 1 million sites containing gerp scores to find minimum and maximum value
low = np.inf
high = -np.inf
sites = 0

# find the overall min/max and number of sites
for chunk in pd.read_csv(gerpfile, sep='\t', chunksize=1000000, header=None, usecols=[3]):
    chunk = chunk.astype(float)
    low = np.minimum(chunk.min(), low)
    high = np.maximum(chunk.max(), high)
    sites += chunk.shape[0]

# define the shape of the histogram
nbins = 50
bin_edges = np.linspace(low, high, nbins + 1).flatten()
total = np.zeros(nbins, np.uint)

# iterate again over the file in chunks of 1 million sites
for chunk in pd.read_csv(gerpfile, sep='\t', chunksize=1000000, header=None, usecols=[3]):
    chunk = chunk.astype(float)
    # compute bin counts
    subtotal, e = np.histogram(chunk, bins=bin_edges)
    # accumulate bin counts over chunks
    total += subtotal.astype(np.uint)

# plot the histogram
cividis = cm.get_cmap('cividis', 3) # create the cividis colormap
plt.hist(bin_edges[:-1], bins=bin_edges, weights=total, color=cividis.colors[1]) # plot the histogram
plt.xlabel('GERP score') # x axis
plt.ylabel('Number of sites') # y axis
plt.savefig(histpdf, bbox_inches='tight', format='pdf')
