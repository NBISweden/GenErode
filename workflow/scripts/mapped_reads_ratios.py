#!/usr/bin/python3

"""
Script to calculate the ratios of reads mapped to other species' mitogenomes
vs. reads mapped to the target species mitogenome, based on samtools flagstat
files. Outliers might indicate contamination.

Input and output files refer to Snakemake directives.
"""

import os

species_files = snakemake.input.species  # flagstat files vs. target species mitogenome
other_files = snakemake.input.other  # flagstat files vs. other species' mitogenomes
outtable = snakemake.output.ratiotables  # output file from Snakemake rule

def get_mapped_reads(inlist):
    """
    Parse the samtools flagstat files from a python list to return a 
    dictionary containing the filename and the number of mapped reads
    """
    readsDict = {}
    for i in inlist:
        bn = os.path.basename(i)
        name = "_".join(bn.split("_", 3)[:3])
        with open(i, "r") as f:
            for line in f:
                if " mapped (" in line:
                    reads = int(line.strip().split()[0])
                    if reads == 0:
                        readsDict[name] = reads + 1
                    else:
                        readsDict[name] = reads
    return readsDict

def ratio_of_mapped_reads(readsDictOtherSpecies, readsDictTargetSpecies):
    """
    Take two dictionaries containing the filename and the number of mapped reads for the same set of reads 
    mapped to a different species and to the target species and return the ratio of mapped reads
    """
    ratio = { k: round(readsDictOtherSpecies[k] / readsDictTargetSpecies[k], 2) for k in readsDictTargetSpecies if k in readsDictOtherSpecies }
    return ratio


with open(outtable, "w") as out:
    # Parse the samtools flagstat files from the different input lists to return a dictionary containing the filename and the number of mapped reads for each of them
    speciesDict = get_mapped_reads(species_files)
    otherDict = get_mapped_reads(other_files)

    # Calculate the ratios of mapped reads for each species to the target species reference
    other_species = ratio_of_mapped_reads(otherDict, speciesDict)

    # Write to output file
    for k, v in other_species.items():
        out.write("{} {}\n".format(k, v))