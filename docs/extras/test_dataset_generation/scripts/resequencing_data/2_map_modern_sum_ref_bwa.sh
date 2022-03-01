#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 4-00:00:00
#SBATCH -M snowy

### Usage: loop through sample IDs to start the script (from within the scripts directory).
# for i in $(cat 2_rhino_modern_samples.txt); do sbatch 2_map_modern_sum_ref_bwa.sh $i; done

########## Modify the following parameters:
fastq="../../testdata/modern" # path to directory with fastq reads. 
bam="../../testdata/modern" # path to output directory for bam files. 
ref="../../testdata/reference/sumatran_rhino_22Jul2017_9M7eS_haploidified_headersFixed.fasta" # enter correct path to reference genome (directory and file name)
refname="sumatran_rhino_full" # insert a short and descriptive name for the reference genome used for mapping. Will be inserted into the file names of all downstream files.
##########

module load bioinfo-tools bwa/0.7.17 samtools/1.10

bwa mem -t 4 ${ref} ${fastq}/${1}_R1.fastq.gz ${fastq}/${1}_R2.fastq.gz | samtools sort -o ${bam}/${1}.${refname}.bam
