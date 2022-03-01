#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:15:00
#SBATCH -M snowy

# Check if any read is duplicated from merging with randomly sampled reads
### Usage: 
# for i in $(cat 2_rhino_modern_samples.txt); do sbatch 7_check_dedup_read_pairs.sh $i modern; done
# for i in $(cat 2_rhino_historical_samples.txt); do sbatch 7_check_dedup_read_pairs.sh $i historical; done

##########
sample=${1}
dir=${2}
##########

cd ../../testdata/${dir}

# Get a list of all read names from the forward reads FASTQ file
gzip -cd ${sample}.extr.R1.fastq.gz | awk 'NR%4==1 {print substr($1,2)}' | sort > ${sample}.extr.R1.readnames.txt &&

# Extract duplicate read names
uniq -d ${sample}.extr.R1.readnames.txt | sort > ${sample}.extr.R1.duplicated_readnames.txt &&

# Create a directory for samples with duplicated read names
if [ ! -d duplicated_reads ]; then mkdir duplicated_reads; fi

# If there are any duplicates, create a list of read names to keep (i.e. subtract list of duplicated read names from all read names)
if [ -s ${sample}.extr.R1.duplicated_readnames.txt ]; then
    # The file is not-empty.
    grep -Fvx -f ${sample}.extr.R1.duplicated_readnames.txt ${sample}.extr.R1.readnames.txt > ${sample}.extr.R1.unique_readnames.txt
    # Move the read list and FASTQ file with duplicated read names into the directory "duplicated_reads/"
    mv ${sample}.extr.R1.fastq.gz ${sample}.extr.R1.unique_readnames.txt duplicated_reads/ &&
    # Remove intermediate files
    rm ${sample}.extr.R1.readnames.txt
else
    # The file is empty.
    echo "There are no duplicated read names" &&
    # Remove intermediate files
    rm ${sample}.extr.R1.*readnames.txt
fi

