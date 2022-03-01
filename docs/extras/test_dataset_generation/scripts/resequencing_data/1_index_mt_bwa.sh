#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 1
#SBATCH -t 01:00:00
#SBATCH -M snowy

########## Modify the following parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
fasta="rhino_NC_012684.fasta" # enter correct reference genome file name (FASTA format)
##########

cd ${refdir}
module load bioinfo-tools bwa/0.7.17

bwa index ${fasta}

