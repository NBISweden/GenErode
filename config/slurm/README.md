# GenErode execution on SLURM clusters

With the switch to Snakemake version 8, GenErode can be run 
the following on SLURM clusters:

1) Create the GenErode conda environment or update an earlier 
version. The latest conda environment contains the Snakemake 
executor plugin for slurm:

```
conda create -f environment.yaml -n generode
```

2) Copy one of the example configuration files `config/slurm/profile/config_plugin_rackham.yaml` 
or `config/slurm/profile/config_plugin_dardel.yaml` to 
`slurm/config.yaml`. This file specifies compute resources 
for each rule or group jobs. Any rule or group job that is 
not listed under `set-threads` or `set-resources` uses 
default resources specified under `default-resources`. If 
any rule or group jobs fail due to too little memory or run 
time, their compute resources can be updated in this file. 

> Note that the current configuration files were adjusted to the 
HPC clusters Rackham from UPPMAX and Dardel from PDC/KTH. Details 
on how to configure and run GenErode on Dardel are provided below. 
The configuration file for Snakemake version 7 was kept for comparison 
which was also written for Rackham/UPPMAX. 

3) Start GenErode the following:

- Open a tmux or screen session
- Activate the GenErode conda environment
- Start the dry run:

```
snakemake --profile slurm -np &> YYMMDD_dry.out
```

- Start the main run:

```
snakemake --profile slurm &> YYMMDD_main.out
```

> Useful flags for running the pipeline: `--ri` to re-run 
incomplete jobs and `-k` to keep going in case a job fails. 

## Specific instructions for Dardel

1) Load the following modules on Dardel:

```
module load PDC UPPMAX bioinfo-tools conda singularity tmux
```

2) After cloning the repository, change permissions for the 
Snakefile:

```
chmod 755 Snakefile
```

3) Create the GenErode conda environment or update an earlier 
version. The latest conda environment contains the Snakemake 
executor plugin for slurm:

```
conda create -f environment.yaml -n generode
```

4) Copy the configuration file `config/slurm/profile/config_plugin_dardel.yaml` 
to `slurm/config.yaml`. This file specifies compute resources 
for each rule or group jobs to be run on Dardel. Any rule or 
group job that is not listed under `set-threads` or `set-resources` 
uses default resources specified under `default-resources`. If 
any rule or group jobs fail due to too little memory or run 
time, their compute resources can be updated in this file. 

> Note that the current version of `config/slurm/profile/config_plugin_dardel.yaml` 
is still being tested. Threads are currently specified under 
`set-threads` and under `set-resources` as `cpus_per_task`.

5) Start GenErode the following:

- Open a tmux session (alternatively, you can use screen)

- Activate the GenErode conda environment (create or update 
from `environment.yaml`), replacing the path to the location 
of the conda environment:

```
export CONDA_ENVS_PATH=/cfs/klemming/home/.../
conda activate generode
```

- Start the dry run:

```
snakemake --profile slurm -np &> YYMMDD_dry.out
```

- Start the main run:

```
snakemake --profile slurm &> YYMMDD_main.out
```

> Useful flags for running the pipeline: `--ri` to re-run 
incomplete jobs and `-k` to keep going in case a job fails. 
