# Slurm configuration file templates

This folder contains a few templates to run GenErode with 
a slurm profile that can be used as starting points to 
define compute resources. 

`config_plugin_dardel.yaml` is a template for running 
GenErode on Dardel (PDC/KTH) and Snakemake version 8. 

`config_plugin_rackham.yaml` is a template for the cluster 
Rackham (UPPMAX) and Snakemake version 8. 

`config_plugin_dardel_testdata.yaml` is a template for 
runs with small test datasets on Dardel (PDC/KTH) and 
Snakemake version 8. 

`config_v7.yaml` is a template for running older GenErode 
versions with Snakemake version 7 and the slurm profile from 
https://github.com/Snakemake-Profiles/slurm, based on the 
cluster Rackham. 
