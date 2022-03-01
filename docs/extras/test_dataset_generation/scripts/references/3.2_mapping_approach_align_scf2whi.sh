#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 2
#SBATCH -t 4-00:00:00


########## Input parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
whi="GCF_000283155.1_CerSimSim1.0_genomic" # enter reference fasta file prefix
sum="sumatran_rhino_22Jul2017_9M7eS_haploidified_headersFixed" # Sumatran rhino fasta file prefix
##########

cd ${refdir}

module load bioinfo-tools bwa/0.7.17 samtools/1.9 BEDTools/2.29.2

# Map the Sumatran rhino reads from each target scaffold to the full White rhino reference
bwa mem -R '@RG\tID:sumrhino\tSM:sumrhino\tPL:ILLUMINA\tPI:330' -t 8 ${whi}.fna ${sum}_Sc9M7eS_2_HRSCAF_41.fq.gz | \
    samtools view -@ 8 -h -q 1 -F 4 -F 256 | grep -v XA:Z | grep -v SA:Z | \
    samtools view -@ 8 -b - | samtools sort -@ 8 - > ${whi}_${sum}_Sc9M7eS_2_HRSCAF_41.bam &&

samtools index ${whi}_${sum}_Sc9M7eS_2_HRSCAF_41.bam &&

# Convert bam file to BED format, listing genomic locations where Sumatran rhino reads mapped to White rhino
bedtools bamtobed -i ${whi}_${sum}_Sc9M7eS_2_HRSCAF_41.bam > ${whi}_${sum}_Sc9M7eS_2_HRSCAF_41.bed &&

# Create a histogram of scaffolds, indicating how many reads mapped to each scaffold
awk '{counts[$1]++} END {for (c in counts) print c, counts[c]}' ${whi}_${sum}_Sc9M7eS_2_HRSCAF_41.bed  | sort -nk 2 > ${whi}_${sum}_Sc9M7eS_2_HRSCAF_41.bed.hist