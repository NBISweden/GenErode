#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

### Usage (move all FASTQ files that should be concatenated to the directory `concatenation_input/` before running this script): 
# for i in $(cat 2_rhino_modern_samples.txt); do sbatch 6_concatenate_per_sample.sh $i modern; done
# for i in $(cat 2_rhino_historical_samples.txt); do sbatch 6_concatenate_per_sample.sh $i historical; done


##########
sample=${1}
dir=${2}
##########

cd ../../testdata/${dir}

cat concatenation_input/${1}.*.extr.R1.fastq.gz > ${1}.extr.R1.fastq.gz &&
cat concatenation_input/${1}.*.extr.R2.fastq.gz > ${1}.extr.R2.fastq.gz