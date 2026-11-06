# GenErode execution on Arrhenius (NAISS)

On Arrhenius, Apptainer, tmux and screen are pre-installed 
and do not have to be loaded with module load.

1) After cloning the repository, change permissions for the 
Snakefile:

```
chmod 755 Snakefile
```

2) Install Miniconda in your home directory if you haven't 
done so yet. 

Create the GenErode conda environment, or update an earlier 
version. Note that since GenErode version 0.7.1, this 
environment only contains Snakemake and the slurm plugin, 
so it is much smaller than the conda environment from 
earlier pipeline versions.

```
conda env create -f environment.yml -n generode
```

> Conda environments and packages take a lot of storage space 
> in your home directory. You can specify a storage location 
> for the GenErode conda environment when creating it using 
> the command
> `conda env create -f environment.yml --prefix /nobackup/proj/disk/.../conda-envs/generode` 
> instead of the command above. In that case, you need to 
> provide the full path to the conda environment when activating 
> it (see below). 

> A `pixi.toml` file is provided in the repository so that pixi 
> can be used as an alternative to conda. 

3) Adjust the configuration file `profiles/arrhenius/config.yaml` 
according to your compute resource needs. This file specifies 
Snakemake and apptainer parameters that are used to run the pipeline, 
as well as compute resources for each jobs to be run 
on Arrhenius. When you start the pipeline on the command line with 
`--profile profiles/arrhenius` (as described below), it expects the file 
`config.yaml` in the folder `profiles/arrhenius`. 

The following parameters need to be updated in 
`profiles/arrhenius/config.yaml` before the first run: 

- Add your slurm compute project ID to `slurm_account`. 
On Arrhenius, the account name must end with `-cpu` 
(e.g. `naissXXXX-X-XX-cpu`). 
- Specify a directory in your storage project in which to 
store container `*.simg` image files using the parameter 
`apptainer-prefix`. 
- Mount your project directory and the local scratch to 
the containers using the parameter `apptainer-args` (e.g. 
`"-B /scratch/local -B /nobackup/proj/disk/..."`), so that 
the tools run inside the containers can access input and 
output files. 

Compute resources are defined as follows in 
`profiles/arrhenius/config.yaml`: 

- All jobs use `default-resources`, unless more threads 
(corresponding to cpus-per-task) are specified per job 
under `set-threads`, or longer run times or more memory 
are specified under `set-resources`. 
- Memory is calculated automatically from the number of 
threads, with 2700 MB per thread (the amount of memory per 
core on Arrhenius). If a job fails due to too little memory, 
increase its number of threads under `set-threads`. 
- If a job fails due to too little run time, increase its 
`runtime` (in minutes) under `set-resources`. The maximum 
run time on the Arrhenius `cpu` partition is 4320 minutes 
(3 days). 
- Do not set `nodes` (`-N`) unless you need a full node. 
- Temporary files are written to `/scratch/local` (`tmpdir`), 
as `/tmp`, `/var/tmp` and `/dev/shm` are limited in size on 
Arrhenius.

> Note that Java tools like `Picard` can be optimized by providing a 
> small number of threads and a large amount of memory, but this will 
> bill the number of CPU hours corresponding to the requested amount 
> of memory.

> To run GenErode on the small test dataset (in `.test/data`), 
> use the profile `profiles/arrhenius_test` instead, which uses 
> shorter run times. The same parameters need to be updated in 
> `profiles/arrhenius_test/config.yaml` as described above. 

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
