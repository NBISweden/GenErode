#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 1
#SBATCH -t 01:00:00

### Usage:
# while IFS="," read -r name accession; do sbatch 3_extract_scaffolds.sh $name $accession; done < gerp_30_outgroups_accessions.list

########## Input parameters:
db=${1} # database name from for loop on command line
query="SR_mespa_WR90_jan2019_prots.Sc9M7eS_2_HRSCAF_41"
outdir="../../testdata/gerp/outgroup_Sc9M7eS_2_HRSCAF_41/all_scaffolds"
##########

cd ${outdir}

awk '{print $1}' ../../reciprocal_blast/blastx_${db}_${query}.tophits.recipr.tophits.matches.out.hist > ../../reciprocal_blast/blastx_${db}_${query}.tophits.recipr.tophits.matches.out.list.txt &&

module load bioinfo-tools seqtk/1.2-r101
seqtk subseq ../../ncbi_datasets_downloads/${db}.fa.gz ../../reciprocal_blast/blastx_${db}_${query}.tophits.recipr.tophits.matches.out.list.txt > ${db}.fa
