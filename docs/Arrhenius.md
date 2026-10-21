# GenErode execution on Arrhenius (NAISS)

*Note that this document is still in preparation*

On Arrhenius, Apptainer, tmux and screen are pre-installed 
and do not have to be loaded with module load.

1) After cloning the repository, change permissions for the 
Snakefile:

```
chmod 755 Snakefile
```

2) Install Miniconda in your home directory if you haven't 
done so yet. 

Note that conda environments and packages take a lot of 
storage space. You can configure conda and set storage 
directories for environments and packages by adding the 
following lines to the `.condarc` file in your home directory: 

```
envs_dirs: 
  - /nfs/volume/user/conda_envs  
pkgs_dirs: 
  - /nfs/volume/user/conda_pkgs
```

Create the GenErode conda environment, or update an earlier 
version. 

```
conda env create -f environment.yml -n generode
```

> If you do not want to set a common location for all packages 
> and environments via `.condarc`, you can specify a storage 
> location for the GenErode conda environment when creating 
> it using the command
> `conda env create -f environment.yml --prefix /nobackup/proj/disk/.../conda-envs/generode` 
> instead of the command above. In that case, you need to 
> provide the full path to the conda environment when activating 
> it (see below). 

3) Adjust the configuration file `profiles/arrhenius/config.yaml` 
according to your compute resource needs. This file specifies 
Snakemake and apptainer parameters that are used to run the pipeline, 
as well as compute resources for each rule or group jobs to be run 
on Arrhenius. When you start the pipeline on the command line with 
`--profile profiles/arrhenius` (as described below), it expects the file 
`config.yaml` in the folder `profiles/arrhenius`. 

- Please make sure to add your slurm compute project ID to `slurm_account` 
in `profiles/arrhenius/config.yaml`. 
- Specify a directory in which to store container `*.simg` 
image files using the parameter `apptainer-prefix` in 
`profiles/arrhenius/config.yaml`. 

> Note that in case you are starting the pipeline from within 
> an interactive job, Snakemake downloads the container images 
> per default to `/scratch/local`, which is deleted after the 
> interactive job has ended. 

- If a rule or group job fails due to too little memory or 
run time, their compute resources can be updated in 
`profiles/arrhenius/config.yaml` as follows. Rule or group 
jobs are using `default-resources` unless more threads 
(corresponding to cpus-per-task), more memory, or longer 
run times are required, which are specified per rule or 
group job under `set-threads` or `set-resources`, respectively. 

> Note that Java tools like `Picard` can be optimized by providing a 
> small number of threads and a large amount of memory, but this will 
> bill the number of CPU hours corresponding to the requested amount 
> of memory.

4) Once the setup is complete, start GenErode the following:

- Open a tmux session (alternatively, you can use screen)

- Activate the GenErode conda environment:

```
conda activate generode
```

- Start the dry run:

```
snakemake --profile profiles/arrhenius -n &> YYMMDD_dry.out
```

> If you notice that jobs would be run that were already 
> complete and should not be re-run, try to add the flag 
> `--rerun-triggers mtime`.

- Start the main run:

```
snakemake --profile profiles/arrhenius &> YYMMDD_main.out
```

> Useful flags for running the pipeline: `--ri` to re-run 
> incomplete jobs and `-k` to keep going in case a job fails. 

# Snakemake cache location

The Snakemake cache is stored per default in your home directory
that can run out of storage space. There are two ways to set a 
directory in your storage project as cache directory: 

1) by setting the cache with the following commands (replacing 
the path with the path to your GenErode directory on Arrhenius): 
`mkdir -p /nobackup/proj/disk/.../.../GenErode/.cache` 
`export XDG_CACHE_HOME=/nobackup/proj/disk/.../.../GenErode/.cache`. 
You will have to run the `export` command every time you start a new 
tmux session. Alternatively, you can add the export command to 
your `~/.bashrc` file. 

2) Another option is to start the pipeline run with the following 
command (replacing the path with the path to your GenErode 
directory on Arrhenius): 
`snakemake --profile slurm --directory /nobackup/proj/disk/.../.../GenErode/`. 
This will set the GenErode directory as the working directory 
and will also store the cache there. 
