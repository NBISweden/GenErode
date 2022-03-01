#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00

### Usage:
# while IFS="," read -r name accession; do sbatch 1_download_gerp_30_outgroups.sh $name $accession; done < gerp_30_outgroups_accessions.list

ncbi_dir="../../testdata/gerp/ncbi_datasets_downloads/"
name=$1
accession=$2

conda activate datasets

cd ${ncbi_dir}
mkdir -p ${name} && cd ${name}
datasets download genome accession ${accession} --exclude-gff3 --exclude-protein --exclude-rna --filename ${accession}.zip