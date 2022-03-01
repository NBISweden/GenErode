#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -M snowy

### Usage: Move to the directory containing the fastq files, and loop through the fastq files to start the script
### (scripts located in test_dataset_generation/scripts/resequencing_data)
# for i in *.fastq.gz; do sbatch 1_fastqc_raw.sh $i; done

module load bioinfo-tools FastQC/0.11.5

if [ ! -d fastqc_raw ]; then mkdir fastqc_raw; fi
fastqc -o fastqc_raw --extract ${1}
