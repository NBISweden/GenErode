# Slurm configuration file templates

This folder contains a few templates to run GenErode with 
a slurm profile that can be used as starting points to 
define compute resources. Please go through the *yaml 
file you are planning to use and adjust the number of 
threads per rul under `set-threads` and memory and run time 
under `set-resources` if needed. For example, to map a 
fastq file with 100 million reads to a 6 Gb reference 
genome on Dardel, 60 threads and 60000 Mb memory has 
worked well in the past.

`config_plugin_dardel.yaml` is a template for running 
GenErode on Dardel (PDC/KTH) and Snakemake version 8. 

`config_plugin_rackham.yaml` is a template for the cluster 
Rackham (UPPMAX) and Snakemake version 8. 

`config_plugin_dardel_testdata.yaml` is a template for 
runs with small test datasets on Dardel (PDC/KTH) and 
Snakemake version 8. 
