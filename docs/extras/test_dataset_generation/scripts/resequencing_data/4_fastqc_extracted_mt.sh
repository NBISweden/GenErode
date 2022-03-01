#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

### Usage: Move to the directory containing the fastq files, 
### and loop through the fastq files to start the script (adjust "*.fastq" accordingly if only run on a specific file version). 
# for i in *mt.*extr*.fastq.gz; do sbatch 4_fastqc_extracted_mt.sh $i; done

module load bioinfo-tools FastQC/0.11.5

if [ ! -d fastqc_extracted_mt ]; then mkdir fastqc_extracted_mt; fi
fastqc -o fastqc_extracted_mt --extract ${1}
