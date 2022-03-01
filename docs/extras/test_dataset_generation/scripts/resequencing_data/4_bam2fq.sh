#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 05:00:00

### Usage: Move to the directory containing the BAM files, and loop through all BAM files.
### (scripts located in test_dataset_generation/scripts/resequencing_data)
# for i in $(ls *.bam| sed 's/.bam//g'); do sbatch 4_bam2fq.sh $i; done

module load bioinfo-tools samtools/1.9 FastQC/0.11.5

# Convert from BAM to FASTQ format, removing secondary alignments
samtools collate -@ 2 -Ou ${1}.bam | samtools fastq -@ 2 -1 ${1}.extr.R1.fastq -2 ${1}.extr.R2.fastq -0 /dev/null -s /dev/null -n - &&

# Zip the fastq files
gzip ${1}.extr*R*.fastq
