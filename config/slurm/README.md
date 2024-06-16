# GenErode execution on SLURM clusters

With the switch to Snakemake version 8, GenErode can be run 
the following on SLURM clusters:

1) Install the Snakemake SLURM plugin the first time GenErode 
is run

```
pip install snakemake-executor-plugin-slurm
```

2) Copy the configuration file `config/slurm/profile/config_plugin.yaml` 
to `slurm/config.yaml`. This file specifies compute resources for each 
rule or group jobs. Any rule or group job that is not listed under 
`set-threads` or `set-resources` uses default resources specified under 
`default-resources`. If any rule or group jobs fail due to too little 
memory or run time, their compute resources can be updated in this file. 

> Note that the current configuration file was adjusted to the 
HPC cluster Dardel from PDC/KTH. The configuration file for 
Snakemake version 7 was kept for comparison, which was set up 
for Rackham/UPPMAX with a different cluster configuration.

3) Start GenErode the following:

- Open a tmux or screen session
- Activate the GenErode conda environment (create or update 
from `environment.yaml`)
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

1) Install the Snakemake SLURM plugin the first time GenErode 
is run

```
pip install snakemake-executor-plugin-slurm
```

2) Copy the configuration file `config/slurm/profile/config_plugin.yaml` 
to `slurm/config.yaml`. This file specifies compute resources for each 
rule or group jobs. Any rule or group job that is not listed under 
`set-threads` or `set-resources` uses default resources specified under 
`default-resources`. If any rule or group jobs fail due to too little 
memory or run time, their compute resources can be updated in this file. 

3) Start GenErode the following:

- Load the following modules:

```
module load PDC UPPMAX bioinfo-tools conda singularity tmux
```

- Open a tmux or screen session

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
