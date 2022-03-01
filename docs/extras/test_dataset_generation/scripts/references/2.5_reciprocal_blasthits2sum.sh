#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 19
#SBATCH -t 4-00:00:00


########## Input parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
db="GCF_000283155.1_CerSimSim1.0_genomic.fna" # enter correct reference genome database name
query="SR_mespa_WR90_jan2019_prots" # query file prefix
##########

cd ${refdir}
module load bioinfo-tools blast/2.11.0+
blastx -db ${query}.Sc9M7eS_2_HRSCAF_41.fa -query tblastn_${db}_${query}.Sc9M7eS_2_HRSCAF_41.tophits.fasta -outfmt 6 -evalue 10e-5 -num_threads 19 -out blastx_${db}_${query}.Sc9M7eS_2_HRSCAF_41.tophits.recipr.out
