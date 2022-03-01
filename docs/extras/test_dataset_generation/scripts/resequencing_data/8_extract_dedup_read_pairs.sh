#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

### Usage: 
# for i in $(cat 2_rhino_modern_samples.txt); do sbatch 8_extract_dedup_read_pairs.sh $i modern; done
# for i in $(cat 2_rhino_historical_samples.txt); do sbatch 8_extract_dedup_read_pairs.sh $i historical; done

##########
sample=${1}
dir=${2}
##########

cd ../../testdata/${dir}

module load bioinfo-tools seqtk/1.2-r101

seqtk subseq duplicated_reads/${1}.extr.R1.fastq.gz duplicated_reads/${1}.extr.R1.unique_readnames.txt > ${1}.extr.R1.fastq &&
seqtk subseq duplicated_reads/${1}.extr.R2.fastq.gz duplicated_reads/${1}.extr.R1.unique_readnames.txt > ${1}.extr.R2.fastq