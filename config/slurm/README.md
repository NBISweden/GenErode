# GenErode execution on SLURM clusters

With the switch to Snakemake version 8, GenErode can be run 
the following on SLURM clusters:

1) Install the Snakemake SLURM plugin the first time GenErode 
is run

```
pip install snakemake-executor-plugin-slurm
```

2) Copy the configuration file `config/slurm/profile/config_plugin.yaml` 
to `slurm/config.yaml` and adjust the compute resources according to 
the local cluster

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

## Specific instructions for Dardel

1) Install the Snakemake SLURM plugin the first time GenErode 
is run

```
pip install snakemake-executor-plugin-slurm
```

2) Copy the configuration file `config/slurm/profile/config_plugin.yaml` 
to `slurm/config.yaml` and adjust the compute resources according to 
the local cluster

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