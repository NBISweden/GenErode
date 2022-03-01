#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 1
#SBATCH -t 00:30:00


########## Input parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
db="GCF_000283155.1_CerSimSim1.0_genomic.fna" # enter correct reference genome database name
query="SR_mespa_WR90_jan2019_prots" # query file prefix
##########

cd ${refdir}

# Extract top blast hit per protein and convert to BED format
awk '!x[$1]++' tblastn_${db}_${query}.Sc9M7eS_2_HRSCAF_41.out | awk '{print $2, $9, $10, $1}' OFS='\t' > tblastn_${db}_${query}.Sc9M7eS_2_HRSCAF_41.tophits.out &&
awk '{ if ( $2<=$3 ) {print $1, $2-1, $3, $4} else if ( $2>$3 ) {print $1, $3-1, $2, $4} }' OFS='\t' tblastn_${db}_${query}.Sc9M7eS_2_HRSCAF_41.tophits.out > tblastn_${db}_${query}.Sc9M7eS_2_HRSCAF_41.tophits.bed

# Extract sequence in fasta format
module load bioinfo-tools BEDTools/2.29.2
bedtools getfasta -fi ${db} -bed tblastn_${db}_${query}.Sc9M7eS_2_HRSCAF_41.tophits.bed -name > tblastn_${db}_${query}.Sc9M7eS_2_HRSCAF_41.tophits.fasta
