#!/usr/bin/python3

"""
Script to approximate any percentile of GERP scores from the entire genome
using histogram bins, to be used as minimum or maximum GERP value to calculate 
relative mutational load.

Usage:
Provide path to gzipped GERP ancestral.rates file and the desired percentile on the command line. 
Example to approximate the 99th percentile (for the top 1% of GERP scores):

python get_gerp_score_percentile.py ../results/gerp/GCF_000283155.1_CerSimSim1.0_genomic.Sc9M7eS_2_HRSCAF_41.ancestral.rates.gz 99
"""

import sys
import pandas as pd
import numpy as np

gerpfile=sys.argv[1]
percentile=sys.argv[2]

def min_max_gerp(file):
    # read the file in chunks of 1 million sites to find minimum and maximum gerp score
    low = np.inf
    high = -np.inf
    for chunk in pd.read_csv(file, sep='\t', chunksize=1000000, header=None, usecols=[3]):
        chunk = chunk.astype(float)
        low = np.minimum(chunk.min(), low)
        high = np.maximum(chunk.max(), high)
    return low, high

def gerp_hist(file, low, high):
    # create a histogram
    # define the shape of the histogram (i.e. number of bins to be filled with counts)
    nbins = 1000
    # minimum and maximum gerp score per bin. Add 1 because there will be 2 values per bin.
    bin_edges = np.linspace(low, high, nbins + 1).flatten()
    # counts of scores per bin, starting with zero
    total = np.zeros(nbins, np.uint)
    # iterate over the file in chunks of 1 million sites to fill the histogram with values
    for chunk in pd.read_csv(gerpfile, sep='\t', chunksize=1000000, header=None, usecols=[3]):
        chunk = chunk.astype(float)
        # compute bin counts
        subtotal, e = np.histogram(chunk, bins=bin_edges)
        # accumulate bin counts over chunks
        total += subtotal.astype(np.uint)
    return bin_edges, total

def approx_percentile(p, bin_edges, total):
    print(p, "th percentile:")
    # calculate the cumulative sums of histogram counts in percentages
    cs = np.cumsum(total)/np.sum(total) * 100
    # how many values of the cumulative sum in percentage are lower or equal to p? This is used as index in the next step.
    i = len(cs[cs <= int(p)])
    # return the maximum gerp score of the bin that is equal or lower to p% with i as index. All gerp scores higher than that are the top 1-p%.
    print(bin_edges[i])

min_gerp, max_gerp = min_max_gerp(gerpfile)
bins, counts = gerp_hist(gerpfile, min_gerp, max_gerp)
approx_percentile(percentile, bins, counts)