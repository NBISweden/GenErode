#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2-00:00:00

### Usage: loop through sample names to extract random reads
# for i in $(cat 2_rhino_modern_samples.txt); do sbatch 5_sample_random_paired_reads_mt_modern.sh $i; done

module load bioinfo-tools seqtk/1.2-r101

cd ../../testdata/modern

seqtk sample -s100 ${1}.mt.filtered.extr.R1.fastq.gz 5000 > ${1}.5k.mt.filtered.extr.R1.fastq.gz &&
seqtk sample -s100 ${1}.mt.filtered.extr.R2.fastq.gz 5000 > ${1}.5k.mt.filtered.extr.R2.fastq.gz
