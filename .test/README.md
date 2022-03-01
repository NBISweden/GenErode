# Testing with GitHub actions

This directory contains config files and data used to test 
the GenErode pipeline with GitHub actions. The dataset can 
also be used to try out the different options of the 
GenErode pipeline, but note that results are biologically 
not meaningful due to the small size of the dataset.

## config/

`config/` contains files that are used by the different 
test steps of GitHub actions.

## data/

`data/` contains two different datasets used to test the 
GenErode pipeline:

### FASTA and GTF files used as references for testing (`data/references`)

* `sumatran_rhino.fasta`: one scaffold extracted from the 
  full Sumatran rhinoceros reference genome (accession number 
  GCA_014189135.1) for testing of the pipeline. 
* `sumatran_rhino.gtf`: corresponding GFF file for testing of 
  the snpEff step of the pipeline. Conversion from GFF3 format 
  was done using the tool `agat_convert_sp_gff2gtf.pl` from 
  the AGAT toolkit (https://github.com/NBISweden/AGAT)
* To test features of the GERP step of the pipeline, the FASTA 
  and GTF files were split into three sequences using the tool 
  `SplitFastaAndGFF` from NBIS-UtilityCode 
  (https://github.com/NBISweden/NBIS-UtilityCode).
* `rhino_NC_012684.fasta`: mitochondrial genome of the Sumatran 
  rhinoceros in FASTA format for testing of the mitochondrial 
  mapping step of the pipeline.


### Short read data (`data/shortread_data`)

* `historical/`: Illumina short reads for one historical sample 
  in FASTQ format, from two independent library preparations and 
  sequenced on one lane.
* `modern/`: Illumina short reads for one modern sample in FASTQ 
  format.
* FASTQ files were generated the following: Paired-end reads from 
  the original FASTQ files were mapped to the full reference FASTA 
  file. The scaffold used as test reference 
  (`data/references/sumatran_rhino.fasta`) was extracted. The 
  resulting BAM files were converted to FASTQ files. The same 
  approach was used to extract reads mapping to the mitochondrial 
  genome, followed by an additional subsampling step to reduce 
  the number of mitochondrial reads.

### GERP++ test dataset (`data/gerp_data`)

* Time-calibrated phylogenetic tree and multiple FASTA files 
  (gzipped) from mammalian species.
* FASTA files from mammalian species were generated as following: 
  CDS of genes located on the Sumatran rhinoceros test scaffold were 
  blasted to each of the mammalian genomes using tblastx. The top 
  hit was extracted per CDS and blasted to the Sumatran rhinoceros 
  CDS sequences using blastn. Scaffolds with at least one top hit 
  from the reciprocal blast were extracted from each of the mammalian 
  genomes. Compressed FASTA files larger than 50 Mb were removed 
  from the test dataset.
* The time-calibrated phylogenetic tree was generated on 
  http://timetree.org from a list of all mammalian species and the 
  Sumatran rhinoceros.

#### Mammalian genomes and accession numbers

Suricata_suricatta GCF_006229205.1  
Hyaena_hyaena GCF_003009895.1  
Paradoxurus_hermaphroditus GCA_004024585.1  
Panthera_leo GCA_008795835.1  
Canis_lupus GCF_014441545.1  
Procyon_lotor GCA_015708975.1  
Ailurus_fulgens GCA_002007465.1  
Spilogale_gracilis GCA_004023965.1  
Leptonychotes_weddellii GCF_000349705.1  
Ursus_maritimus GCF_017311325.1  
Manis_javanica GCF_014570535.1  
Diceros_bicornis GCA_013634535.1  
Tapirus_indicus GCA_004024905.1  
Equus_asinus GCA_016077325.1  
Camelus_dromedarius GCF_000803125.2  
Hippopotamus_amphibius GCA_004027065.2  
Balaenoptera_acutorostrata GCF_000493695.1  
Physeter_catodon GCF_002837175.2  
Mesoplodon_bidens GCA_004027085.1  
Lipotes_vexillifer GCF_000442215.1  
Tursiops_truncatus GCF_011762595.1  
Tragulus_javanicus GCA_004024965.2  
Bubalus_bubalis GCF_003121395.1  
Ovis_aries GCF_002742125.1  
Antilocapra_americana GCA_007570785.1  
Cervus_elaphus GCA_002197005.1  
Catagonus_wagneri GCA_004024745.2  
