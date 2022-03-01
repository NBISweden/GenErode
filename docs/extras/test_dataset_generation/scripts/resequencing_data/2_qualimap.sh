#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 4-00:00:00

### Usage: Move to the directory containing the bam files, and loop through the bam files to start the script  (adjust "*.bam" accordingly if only run on a specific bam file version). 
# for i in *.bam; do sbatch 2_qualimap.sh $i; done

module load bioinfo-tools QualiMap/2.2.1
unset DISPLAY

qualimap bamqc -bam ${1} --java-mem-size=24G -outdir ${1}.qualimap -outformat html
