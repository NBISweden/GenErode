#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00

### Usage: Move to the directory containing the BAM files, and loop through all BAM files.
### (scripts located in test_dataset_generation/scripts/resequencing_data)
# for i in $(ls *.sumatran_rhino_full.bam | sed 's/.bam//g'); do sbatch 3_extract_sum_ref_scf_from_bams.sh $i; done

module load bioinfo-tools samtools/1.9

samtools index ${1}.bam &&
samtools view -h ${1}.bam Sc9M7eS_2_HRSCAF_41 > ${1}.Sc9M7eS_2_HRSCAF_41.sam &&
samtools view -bS ${1}.Sc9M7eS_2_HRSCAF_41.sam > ${1}.Sc9M7eS_2_HRSCAF_41.samtools.bam &&
rm ${1}.Sc9M7eS_2_HRSCAF_41.sam