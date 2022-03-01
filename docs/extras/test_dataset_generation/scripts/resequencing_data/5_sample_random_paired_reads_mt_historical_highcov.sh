#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2-00:00:00

### Usage: loop through sample names to extract random reads
# for i in JvS008_08_L2 JvS008_08_L6 JvS008_10_L2 JvS008_10_L6 JvS008_11_L2 JvS008_11_L6 JvS009_09_L3 JvS009_09_L7 JvS009_15_L3 JvS009_15_L7 JvS009_19_L3 JvS009_19_L7; do sbatch 5_sample_random_paired_reads_mt_historical_highcov.sh $i; done

module load bioinfo-tools seqtk/1.2-r101

cd ../../testdata/historical

seqtk sample -s100 ${1}.mt.filtered.extr.R1.fastq.gz 1000 > ${1}.1k.mt.filtered.extr.R1.fastq.gz &&
seqtk sample -s100 ${1}.mt.filtered.extr.R2.fastq.gz 1000 > ${1}.1k.mt.filtered.extr.R2.fastq.gz