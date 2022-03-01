#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00


########## Input parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
whi="GCF_000283155.1_CerSimSim1.0_genomic" # enter reference fasta file prefix
##########

cd ${refdir}
module load bioinfo-tools seqtk/1.2-r101
seqtk subseq ${whi}.fna ${whi}.Sc9M7eS_2_HRSCAF_41.list.txt > ${whi}.Sc9M7eS_2_HRSCAF_41.fasta