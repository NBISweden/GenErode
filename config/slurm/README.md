# GenErode execution on Dardel (PDC/KTH)

1) Load the following modules:

```
module load PDC bioinfo-tools apptainer tmux
```

> Note that `tmux` is only available as a module on Dardel 
but the equivalent tool `screen` is pre-installed and does 
not need to be loaded. 

> Apptainer (former singularity) will store its cache per 
default in your home directory which will quickly run out of 
storage space. You can tell it to use your `scratch` instead, a 
temporary directory with unlimited space by adding this row 
to your `~/.bashrc`: `export APPTAINER_CACHEDIR=$PDC_TMP`. 
The files in this temporary directory are deleted if they have 
not been used for 30 days. Alternatively, you can set the cache 
directory to a directory in your storage project, adding this 
line to your `~/.bashrc` (replacing the path to an existing 
directory in your storage project):
`export APPTAINER_CACHEDIR=/cfs/klemming/projects/supr/sllstore.../.../apptainer-cache` 

2) After cloning the repository, change permissions for the 
Snakefile:

```
chmod 755 Snakefile
```

3) Create the GenErode conda environment or update an earlier 
version. The latest conda environment contains the Snakemake 
executor plugin for slurm. Since home directories on Dardel 
are limited in storage space, you need to create a directory in 
your storage project for the conda environment to be installed 
in, and run the following command instead of the command above
(replacing the path with the path to the directory in your storage
project for the conda environment): 

```
export CONDA_ENVS_PATH=/cfs/klemming/projects/supr/sllstore.../.../conda-envs
conda env create -f environment.yml -p /cfs/klemming/projects/supr/sllstore.../.../conda-envs/generode
```

> Note that you can save storage space in your storage project 
on Dardel by creating a common GenErode conda environment for 
several people. 

4) Copy the configuration file `config/slurm/profile/config_plugin_dardel.yaml` 
to `slurm/config.yaml`. This file specifies Snakemake and apptainer 
parameters that are used to run the pipeline, as well as compute 
resources for each rule or group jobs to be run on Dardel. When 
you start the pipeline on the command line with `--profile slurm` 
(as described below), it expects a folder `slurm` with the file
`config.yaml` therein. Please make sure to add your slurm compute 
project ID in line 13 of `slurm/config.yaml` (`slurm_account`). 

> If a rule or group job fails due to too little memory or run time,
their compute resources can be updated in `slurm/config.yaml`. 
Rule or group jobs are using `default-resources` unless more threads
(corresponding to cpus-per-task) or longer run times are required,
which are specified per rule or group job under `set-threads` or
`set-resources`, respectively. Memory in MB is automatically calculated
from the number of threads specified under `default-resources` or
`set-threads`, respectively.  

5) Start GenErode the following:

- Open a tmux session (alternatively, you can use screen)

- Activate the GenErode conda environment with the 
following commands, replacing the path to the location of 
the conda environment accordingly:

```
export CONDA_ENVS_PATH=/cfs/klemming/projects/supr/sllstore.../.../conda-envs/
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
