#! /bin/bash -l
#SBATCH -A snic2021-5-47
#SBATCH -p core -n 1
#SBATCH -t 1-00:00:00

########## Input parameters:
refdir="../../testdata/reference" # enter correct path to reference genome (only the directory)
ref="GCF_000283155.1_CerSimSim1.0_genomic.fna" # enter correct reference genome file name (fasta format)
##########

cd ${refdir}
module load bioinfo-tools blast/2.11.0+
makeblastdb -in ${ref} -dbtype 'nucl' -input_type fasta