#!/usr/bin/python3

"""
Author: Julius Falck

Script to concatenate all depth statistics files from the bam_processing step into one file,
including mean genome-wide depth, minimum depth threshold and maximum depth threshold as 
calculated in GenErode for depth filtering.

Usage:
Move into the directory "GenErode/" and start the script from the command line with 
the following command, replacing "[dataset]" with modern or historical and "[reference_name]" 
with the name of the reference genome used for mapping:

python summarize_depth.py results/[dataset]/mapping/[reference_name]/stats/bams_indels_realigned/

e.g. 
python summarize_depth.py results/modern/mapping/sumatran_rhinoceros/stats/bams_indels_realigned/
"""

import os
import sys

directory = sys.argv[1]
dir_list = os.listdir(directory)

depth_table = "depth_table.txt"

# Add header line to the depth_table file
with open(depth_table, 'w') as file:
    file.write("sample mean_depth minimum_depth_threshold maximum_depth_threshold\n")

for f in dir_list:
    if f.endswith('dpstats.txt'):
        with open(directory + '/' + f, 'r') as file:
            line = file.readlines()[0]
            with open(depth_table, 'a') as file:
                file.write(f.split('.')[0] + ' ' + line + '\n')

