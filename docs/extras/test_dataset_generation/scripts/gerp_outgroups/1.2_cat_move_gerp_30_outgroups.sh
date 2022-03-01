#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00

### Usage: move to testdata/gerp/ncbi_datasets_downloads and run (scripts located in test_dataset_generation/scripts/gerp_outgroups)
# for i in */ncbi_dataset/data/*/; do sbatch 1.2_cat_move_gerp_30_outgroups.sh $i; done

name=`echo ${1} | cut -d"/" -f1`
cat ${1}/*.fna | gzip - > ${name}.fa.gz
