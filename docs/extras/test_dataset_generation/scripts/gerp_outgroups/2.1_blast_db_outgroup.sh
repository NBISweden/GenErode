#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 1
#SBATCH -t 1-00:00:00

### Usage: move to testdata/gerp/ncbi_datasets_downloads and start the script
### by looping through each genome *fa.gz file
### (scripts located in test_dataset_generation/scripts/gerp_outgroups)
# for i in *.fa.gz; do sbatch 2.1_blast_db_outgroup.sh $i; done

db=`echo ${1} | sed 's/.fa.gz/_blastDB/g'`

module load bioinfo-tools blast/2.11.0+
gunzip -c ${1} | makeblastdb -in - -dbtype 'nucl' -input_type fasta -title ${db} -out ${db}
