# Workflow to obtain outgroup references for GERP analyses with the test dataset

30 mammalian outgroup genomes were downloaded from NCBI and 
subset to genome regions with reciprocal blast hits 
with the Sumatran rhinoceros proteins from scaffold 
`Sc9M7eS_2_HRSCAF_41` (protein sequences from annotation 
of the Sumatran rhinoceros genome assembly with GenBank 
accession number GCA_014189135.1, as described in Lord 
et al. 2020).

List of species and accession numbers: 
`gerp_30_outgroups_accessions.list`

All scripts were run on resources provided by the Swedish 
National Infrastructure for Computing (SNIC) at UPPMAX,
using the "module" system or conda environments.


## Identification and extraction of putatively ortholog scaffolds via reciprocal blast

### Download, unzip, and move genome assemblies from 30 mammalian genomes

Scripts: `scripts/gerp_outgroups/1_download_gerp_30_outgroups.sh` 
         (requires conda environment from `scripts/gerp_outgroups/datasets_env.yml`), 
         `scripts/gerp_outgroups/1.1_unzip_gerp_30_outgroups.sh`,
         `scripts/gerp_outgroups/1.2_cat_move_gerp_30_outgroups.sh`


### Blast the Sumatran rhino scaffold proteins to each of the 30 outgroup genomes

- Requires the proteins from Sumatran rhinoceros scaffold
  `Sc9M7eS_2_HRSCAF_41` that were extracted as part of
  `workflows/references.md`

Script: `scripts/references/1.2_extract_sum_prot.sh`

- Create a blast database for each of the 30 mammalian 
  genomes

Script: `scripts/gerp_outgroups/2.1_blast_db_outgroup.sh`

- Blast the proteins from Sumatran rhinoceros scaffold 
  `Sc9M7eS_2_HRSCAF_41` to each of the 30 databases

Script: `scripts/gerp_outgroups/2.2_blast_prot2outgroup.sh`

- Extract top blast hit per protein and mammalian genome in FASTA format

Script: `scripts/gerp_outgroups/2.3_extract_prot2outgroup_blast_hits.sh`


### Blast the top hits back to the Sumatran rhinoceros proteins 

- Requires the blast database from the Sumatran rhinoceros 
  proteins that was created in `workflows/references.md`

Script: `scripts/references/2.4_blast_db_sum_prot.sh`

- Blast the top hits per protein and mammalian genome back 
  to the Sumatran rhinoceros protein database 

Script: `scripts/gerp_outgroups/2.4_reciprocal_blasthits2sum.sh`


### Filter the results and create a histogram to identify putatively orthologous scaffolds

- Filter the output to get the top hit per reciprocally 
  blasted protein per mammalian genome & create a histogram 
  from the scaffold names in the list of top hits

Script: `scripts/gerp_outgroups/2.5_reciprocal_blast_evaluation.sh`


### Extract all scaffolds with reciprocal blast hits from the 30 outgroup genomes in FASTA format

Script: `scripts/gerp_outgroups/3_extract_scaffolds.sh `


## References

Lord E, Dussex N, Kierczak M, Díez-del-Molino D, Ryder OA, 
Stanton DWG, et al. Pre-extinction Demographic Stability 
and Genomic Signatures of Adaptation in the Woolly Rhinoceros.
Curr Biol 2020;30:3871-3879.