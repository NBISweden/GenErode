# GenErode execution on Dardel (PDC/KTH)

1) Load the following modules:

```
module load PDC bioinfo-tools apptainer tmux
```

> Note that `tmux` is only available as a module on Dardel 
but the equivalent tool `screen` is pre-installed and does 
not need to be loaded. 

> The Snakemake cache is stored per default in your home directory
that can run out of storage space. There are two ways to set a 
directory in your storage project as cache directory: 
1) by setting the cache with the following commands (replacing 
the path with the path to your GenErode directory on Dardel): 
`mkdir -p /cfs/klemming/projects/supr/sllstore.../.../GenErode/.cache` 
`export XDG_CACHE_HOME=/cfs/klemming/projects/supr/sllstore.../.../GenErode/.cache` 
You will have to run the `export` command every time you start a new 
tmux session. Alternatively, you can add the export command to 
your `~/.bashrc` file. 
2) Another option is to start the pipeline run with the following 
command (replacing the path with the path to your GenErode 
directory on Dardel): 
`snakemake --profile slurm --directory /cfs/klemming/projects/supr/sllstore.../.../GenErode/`. 
This will set the GenErode directory as the working directory 
and will also store the cache there. 

> Apptainer (former singularity) will also store its cache per 
default in your home directory. You can tell it to use your 
`scratch` instead, a temporary directory with unlimited spac, 
by adding this row to your `~/.bashrc`: `export APPTAINER_CACHEDIR=$PDC_TMP`. 
Note that the files in this temporary directory are deleted if 
they have not been used for 30 days. Alternatively, you can set 
the cache directory to a directory in your storage project, 
adding this line to your `~/.bashrc` (replacing the path to an existing 
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
(corresponding to cpus-per-task), more memory, or longer run times 
are required,which are specified per rule or group job under `set-threads`
or `set-resources`, respectively. 

> Note that software that can be run with multi-threading (e.g. `fastp`, 
`bwa`, `samtools`) require 2-4X more memory (specified with `mem_mb` 
under `set-resources`) than the number of threads. Java tools like 
`Picard` can be optimized by providing a very small number of 
threads and a large amount of memory. 

5) Once the setup (1-4) is complete, start GenErode the following:

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
