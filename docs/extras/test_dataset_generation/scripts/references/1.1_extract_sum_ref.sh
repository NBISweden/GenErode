#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 1
#SBATCH -t 01:00:00

########## Modify the following parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
fasta="sumatran_rhino_22Jul2017_9M7eS_haploidified_headersFixed.fasta" # enter correct reference genome file name (fasta format)
out="sumatran_rhino_22Jul2017_9M7eS_haploidified_headersFixed"
##########

cd ${refdir}
module load bioinfo-tools samtools/1.9

samtools faidx ${fasta} Sc9M7eS_2_HRSCAF_41 > ${out}_Sc9M7eS_2_HRSCAF_41.fasta
