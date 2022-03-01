#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 1
#SBATCH -t 01:00:00

########## Modify the following parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
fasta="SR_mespa_WR90_jan2019_prots.fa" # enter correct reference genome file name (fasta format)
out="SR_mespa_WR90_jan2019_prots"
##########

cd ${refdir}
module load bioinfo-tools seqtk/1.2-r101

seqtk subseq ${fasta} ${out}.Sc9M7eS_2_HRSCAF_41.list.txt > ${out}.Sc9M7eS_2_HRSCAF_41.fa
