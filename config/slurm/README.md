# GenErode execution on Dardel (PDC/KTH)

1) Load the following modules:

```
module load PDC UPPMAX bioinfo-tools conda singularity tmux
```

> Note that tmux is only available as a module on Dardel 
but the equivalent tool screen is pre-installed and does 
not need to be loaded. 

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
any rule or group job fails due to too little memory or run 
time, their compute resources can be updated in this file. 

> Note that memory requirements are specified three times in 
the configuration file: 1) under `set-threads` (used by Snakemake 
to specify threads in rules), 2) under `set-resources` and therein 
under `mem_mb`, specifying the memory in Megabytes (multiplying 
the number of threads with the available memory per thread), 
and 3) under `set-resources` and therein under `cpus-per-task` 
(the same number as specified under `set-threads`, required for 
correct memory assignment on Dardel). 

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
snakemake --profile slurm -n &> YYMMDD_dry.out
```

- Start the main run:

```
snakemake --profile slurm &> YYMMDD_main.out
```

> Useful flags for running the pipeline: `--ri` to re-run 
incomplete jobs and `-k` to keep going in case a job fails. 
