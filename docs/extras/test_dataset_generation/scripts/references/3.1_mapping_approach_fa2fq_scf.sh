#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 2
#SBATCH -t 05:00:00


########## Input parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
scriptdir="../../GenErode/workflow/scripts" # enter correct path to scripts directory
sum="sumatran_rhino_22Jul2017_9M7eS_haploidified_headersFixed" # reference fasta file prefix
##########

cd ${refdir}
conda activate generode

#cp ${sum}_Sc9M7eS_2_HRSCAF_41.fasta ${sum}_Sc9M7eS_2_HRSCAF_41.fa && 
#gzip ${sum}_Sc9M7eS_2_HRSCAF_41.fa &&

python3 ${scriptdir}/fa2fq.py ${sum}_Sc9M7eS_2_HRSCAF_41.fa.gz ${sum}_Sc9M7eS_2_HRSCAF_41.fq.gz
