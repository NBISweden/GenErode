# Workflow to generate references for the test dataset

Genome assemblies of Sumatran rhinoceros (GenBank accession 
number GCA_014189135.1) and White rhinoceros (GenBank 
accession number GCF_000283155.1) were processed as 
follows to extract scaffolds to be used as references 
for the test dataset.

All scripts were run on resources provided by the Swedish 
National Infrastructure for Computing (SNIC) at UPPMAX,
using the "module" system.


## Extraction of Sumatran rhinoceros scaffold `Sc9M7eS_2_HRSCAF_41`

Script: `scripts/references/1.1_extract_sum_ref.sh`


## Identification and extraction of White rhinoceros scaffolds via reciprocal blast & mapping approach

### Reciprocal blast approach

#### Blast Sumatran rhinoceros scaffold proteins to White rhinoceros reference

- Extract proteins located on Sumatran rhinoceros scaffold 
  `Sc9M7eS_2_HRSCAF_41` from genome-wide protein FASTA file
  (protein sequences from annotation of GCA_014189135.1,
  as described in Lord et al. 2020)

Script: `scripts/references/1.2_extract_sum_prot.sh`

- Create a blast database from the White rhinoceros 
  reference genome

Script: `scripts/references/2.1_blast_db_whi_ref.sh`

- Blast the Sumatran rhinoceros proteins to the White 
  rhinoceros genome database

Script: `scripts/references/2.2_blast_prot2whi.sh`

- Extract the top blast hit per Sumatran rhinoceros protein

Script: `scripts/references/2.3_extract_prot2whi_blast_hits.sh`


#### Blast the top hits back to the Sumatran rhinoceros proteins

- Create a blast database from the Sumatran rhinoceros proteins

Script: `scripts/references/2.4_blast_db_sum_prot.sh`

- Blast the top hits to the Sumatran rhinoceros protein database

Script: `scripts/references/2.5_reciprocal_blasthits2sum.sh`


#### Filter the results and create a histogram to identify putatively orthologous scaffolds

- Filter the output to get the top hits per reciprocally blasted 
  protein & create a histogram from the scaffold names in the 
  list of top hits from the reciprocal blast

Script: `scripts/references/2.6_reciprocal_blast_evaluation.sh`

- Top scaffolds from histograms plus number of protein top hits: 
  `NW_004454182.1` (114), `NW_004454248.1` (17), `NW_004454260.1` (31)


### Mapping approach

#### Convert the Sumatran rhinoceros scaffold to FASTQ format

Script: `scripts/references/3.1_mapping_approach_fa2fq.sh`

- Requires the python script `fa2fq.py` from the GenErode pipeline


#### Map to White rhinoceros reference genome & analyze the results

- Map to White rhinoceros reference genome
- Convert BAM to BED format
- Create a histogram to identify putatively orthologous scaffolds

Script: `scripts/references/3.2_mapping_approach_align_scf2whi.sh`

- Top scaffolds from histograms: `NW_004454182.1`, `NW_004454248.1`,
  `NW_004454260.1`


### Conclusion

- Both approaches find the same scaffolds
- `Sc9M7eS_2_HRSCAF_41` is 40,842,778 bp, the White rhinoceros 
  scaffolds sum up to 41,195,616 bp (100.9 % of the Sumatran 
  rhinoceros scaffold)


### Extract the putatively orthologous White rhinoceros scaffolds in FASTA format

Script: `scripts/references/4_concatenate_whi_scf.sh`


## GTF file edits

- Extract the scaffold `Sc9M7eS_2_HRSCAF_41` from the Sumatran 
  rhinoceros annotation:

`grep "Sc9M7eS_2_HRSCAF_41" SR_mespa_WR90_jan2019.gtf > SR_mespa_WR90_jan2019_Sc9M7eS_2_HRSCAF_41.gtf`


- Extract the scaffolds `NW_004454248.1`, `NW_004454260.1`, and 
  `NW_004454182.1` from the White rhinoceros annotation:

```
grep "NW_004454182.1" GCF_000283155.1_CerSimSim1.0_genomic.gtf > GCF_000283155.1_CerSimSim1.0_genomic.Sc9M7eS_2_HRSCAF_41.gtf &&
grep "NW_004454248.1" GCF_000283155.1_CerSimSim1.0_genomic.gtf >> GCF_000283155.1_CerSimSim1.0_genomic.Sc9M7eS_2_HRSCAF_41.gtf &&
grep "NW_004454260.1" GCF_000283155.1_CerSimSim1.0_genomic.gtf >> GCF_000283155.1_CerSimSim1.0_genomic.Sc9M7eS_2_HRSCAF_41.gtf
```

## References

Lord E, Dussex N, Kierczak M, Díez-del-Molino D, Ryder OA, 
Stanton DWG, et al. Pre-extinction Demographic Stability 
and Genomic Signatures of Adaptation in the Woolly Rhinoceros.
Curr Biol 2020;30:3871-3879.