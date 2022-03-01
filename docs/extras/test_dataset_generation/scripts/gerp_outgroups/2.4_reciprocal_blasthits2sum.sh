#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 19
#SBATCH -t 4-00:00:00

### Usage:
# while IFS="," read -r name accession; do sbatch 2.4_reciprocal_blasthits2sum.sh $name $accession; done < gerp_30_outgroups_accessions.list

########## Input parameters:
db=${1} # database name from for loop on command line
query_path="../../testdata/reference/SR_mespa_WR90_jan2019_prots.Sc9M7eS_2_HRSCAF_41.fa" # path to query file
query="SR_mespa_WR90_jan2019_prots.Sc9M7eS_2_HRSCAF_41"
outdir="../../testdata/gerp/reciprocal_blast/"
##########

cd ${outdir}
module load bioinfo-tools blast/2.11.0+
blastx -db ${query_path} -query tblastn_${db}_${query}.tophits.fasta -outfmt 6 -evalue 10e-5 -num_threads 19 -out blastx_${db}_${query}.tophits.recipr.out
