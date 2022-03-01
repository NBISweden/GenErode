# Workflow to generate a test dataset from whole-genome re-sequencing data

Whole genome re-sequencing data from three historical and three 
modern Sumatran rhinoceros samples from the now-extinct Malay 
Peninsula population were subset to paired-end reads
mapping to the Sumatran rhinoceros scaffold `Sc9M7eS_2_HRSCAF_41`
and a random selection of reads mapping elsewhere to the genome 
assembly or to the mitochondrial genome (GenBank accession number 
NC_012684.1).

All scripts were run on resources provided by the Swedish 
National Infrastructure for Computing (SNIC) at UPPMAX,
using the "module" system.


## Prepare the data for mapping

- Run fastQC on the untrimmed paired-end reads

Script: `scripts/resequencing_data/1_fastqc_raw.sh`

- Index the Sumatran rhinoceros mitochondrial genome (GenBank 
  accession number NC_012684.1)the Sumatran rhinoceros genome 
  assembly (GenBank accession number GCA_014189135.1) with bwa

Scripts: `scripts/resequencing_data/1_index_mt_bwa.sh` and
         `scripts/resequencing_data/1_index_sum_ref_bwa.sh`


## Map paired-end reads from the untrimmed whole genome re-sequencing data to the references

- Map paired-end reads from the untrimmed fastq files to the 
  reference fasta files (Sumatran rhinoceros genome, Sumatran 
  rhinoceros mitochondrial genome) with bwa mem

Scripts: `scripts/resequencing_data/2_map_historical_sum_ref_bwa.sh`,
         `scripts/resequencing_data/2_map_historical_sum_mt_bwa.sh`,
         `scripts/resequencing_data/2_map_modern_sum_ref_bwa.sh`,
         `scripts/resequencing_data/2_map_modern_sum_mt_bwa.sh`

- Run QualiMap on the BAM files and MultiQC on the QualiMap results

Scripts: `scripts/resequencing_data/2_qualimap.sh`,
         `scripts/resequencing_data/2_multiqc.sh`

- Extract scaffold `Sc9M7eS_2_HRSCAF_41` from the BAM files

Script: `scripts/resequencing_data/3_extract_sum_ref_scf_from_bams.sh`


## Convert to FASTQ format and add random mitochondrial and genomic reads

- Convert the BAM files (subset to `Sc9M7eS_2_HRSCAF_41` and those 
  based on the mitochondrial genome) to FASTQ format, removing 
  secondary alignments

Script: `scripts/resequencing_data/4_bam2fq.sh`

- Create FastQC reports for the FASTQ files

Scripts: `scripts/resequencing_data/4_fastqc_extracted_mt.sh`,
         `scripts/resequencing_data/4_fastqc_extracted_scf.sh`


## Extract random reads mapping to the mitochondrial or whole genome

- Extract randomly selected reads that had mapped to the mitochondrial 
  genome or the whole genome, adjust the numbers of extracted reads 
  according to the average read number per original, genome-wide
  FASTQ file

Scripts: `scripts/resequencing_data/5_sample_random_paired_reads_scf_historical_highcov.sh`,
         `scripts/resequencing_data/5_sample_random_paired_reads_scf_historical_lowcov.sh`,
         `scripts/resequencing_data/5_sample_random_paired_reads_scf_modern.sh`,
         `scripts/resequencing_data/5_sample_random_paired_reads_mt_historical_highcov.sh`,
         `scripts/resequencing_data/5_sample_random_paired_reads_mt_historical_lowcov.sh`,
         `scripts/resequencing_data/5_sample_random_paired_reads_mt_modern.sh`


## Concatenate all FASTQ files per sequencing library

- Move all FASTQ files that should be concatenated into a separate
  directory (`concatenation_input`), i.e. FASTQ files with reads
  mapped to scaffold `Sc9M7eS_2_HRSCAF_41` and FASTQ files with 
  randomly extracted reads that had mapped to the mitochondrial and 
  remaining genome, and run the following script to concatenate 
  the FASTQ files per sequencing library

Script: `scripts/resequencing_data/6_concatenate_per_sample.sh`


## Identify FASTQ files with duplicated reads

- The random reads extracted from the original, genome-wide FASTQ
  files contained some reads that later mapped to scaffold
  `Sc9M7eS_2_HRSCAF_41` and that were therefore duplicated in the
  concatenated FASTQ files

- Identify samples with duplicated reads and move their concatenated
  FASTQ files along with lists of duplicated read names to a
  separate directory

Script: `scripts/resequencing_data/7_check_dedup_read_pairs.sh`

- Remove the duplicated reads from the concatenated FASTQ files
  and move them back to the other concatenated FASTQ files

Script: `scripts/resequencing_data/8_extract_dedup_read_pairs.sh`
