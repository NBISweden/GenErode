#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2-00:00:00

### Usage: loop through sample names to extract random reads
# for i in $(cat 2_rhino_modern_samples.txt); do sbatch 5_sample_random_paired_reads_scf_modern.sh $i; done

module load bioinfo-tools seqtk/1.2-r101

cd ../../testdata/modern

# Sample random reads from original raw data to add some random noise
seqtk sample -s999 ${1}_R1.fastq.gz 42000 > ${1}.42k.extr.R1.fastq.gz &&
seqtk sample -s999 ${1}_R2.fastq.gz 42000 > ${1}.42k.extr.R2.fastq.gz