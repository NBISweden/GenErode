#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2-00:00:00

### Usage: loop through sample names to extract random reads
# for i in JvS022_01_L1 JvS022_02_L1 JvS022_03_L1 JvS022_04_L1 JvS022_05_L1 JvS022_06_L1 JvS022_07_L1 JvS022_08_L1 JvS022_09_L1 JvS022_10_L1 JvS022_11_L1 JvS022_12_L1 JvS022_74_L8 JvS022_75_L8 JvS022_76_L8 JvS022_77_L8 JvS022_78_L8 JvS022_79_L8 JvS022_80_L8 JvS022_81_L8 JvS022_82_L8 JvS022_83_L8 JvS022_84_L8 JvS022_85_L8; do sbatch 5_sample_random_paired_reads_scf_historical_lowcov.sh $i; done

module load bioinfo-tools seqtk/1.2-r101

cd ../../testdata/historical

# Sample random reads from original raw data to add some random noise
seqtk sample -s999 ${1}_R1.fastq.gz 2000 > ${1}.2k.extr.R1.fastq.gz &&
seqtk sample -s999 ${1}_R2.fastq.gz 2000 > ${1}.2k.extr.R2.fastq.gz