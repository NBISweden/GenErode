#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

# Usage: move to the directory with QualiMap results and start the script from there

module load bioinfo-tools MultiQC/1.9
multiqc .
