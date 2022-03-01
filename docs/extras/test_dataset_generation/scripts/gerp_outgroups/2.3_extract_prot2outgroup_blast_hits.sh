#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 1
#SBATCH -t 00:30:00

### Usage:
# while IFS="," read -r name accession; do sbatch 2.3_extract_prot2outgroup_blast_hits.sh $name $accession; done < gerp_30_outgroups_accessions.list

########## Input parameters:
db=${1} # database name from for loop on command line
query="SR_mespa_WR90_jan2019_prots.Sc9M7eS_2_HRSCAF_41"
outdir="../../testdata/gerp/reciprocal_blast"
##########

cd ${outdir}

# Extract top blast hit per protein and convert to BED format
awk '!x[$1]++' tblastn_${db}_${query}.out | awk '{print $2, $9, $10, $1}' OFS='\t' > tblastn_${db}_${query}.tophits.out &&
awk '{ if ( $2<=$3 ) {print $1, $2-1, $3, $4} else if ( $2>$3 ) {print $1, $3-1, $2, $4} }' OFS='\t' tblastn_${db}_${query}.tophits.out > tblastn_${db}_${query}.tophits.bed &&

# Extract sequence in fasta format
module load bioinfo-tools BEDTools/2.29.2
gzip -cd ../ncbi_datasets_downloads/${db}.fa.gz > ${db}.fa &&
bedtools getfasta -fi ${db}.fa -bed tblastn_${db}_${query}.tophits.bed -name > tblastn_${db}_${query}.tophits.fasta &&
rm ${db}.fa
