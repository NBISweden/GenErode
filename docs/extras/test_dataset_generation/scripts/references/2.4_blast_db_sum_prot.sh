#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 1
#SBATCH -t 02:00:00


########## Input parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
ref="SR_mespa_WR90_jan2019_prots.fa" # enter correct reference genome file name (fasta format)
##########

cd ${refdir}
module load bioinfo-tools blast/2.11.0+
makeblastdb -in ${ref} -dbtype 'prot' -input_type fasta
