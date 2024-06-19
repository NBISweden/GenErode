# Mutational load pipeline

Snakemake pipeline to process and analyse VCF files that 
were annotated with snpEff to calculate potential and 
realised mutational load.

Requirements:

- Conda/mamba. To run the pipeline, the GenErode conda 
environment needs to be activated.

- Singularity/apptainer. The pipeline executes each rule 
in a container.

- A slurm profile is provided with the pipeline for runs
on HPC clusters with the slurm workload manager. Please 
add your compute project to the file 
`utilities/mutational_load_snpeff/slurm/settings.json` 
by modifying the following line:

```
    "SBATCH_DEFAULTS": "account=XXXX-XX-XXXX",
```

For job- and cluster-specific adjustments, please edit the 
file `utilities/mutational_load_snpeff/slurm/config.yaml` 
accordingly.

The pipeline performs the following filtering and data 
processing steps:

- Identify sites that are are fixed for the ALT allele across 
all samples (1/1) from a merged and filtered VCF file 
containing all samples. When the reference genome of an 
outgroup species/population is used, these sites are not 
informative for population-level analyses. This filter 
currently only handles VCF files without any missing data. 
The pipeline produces additionally a merged VCF file from 
which these sites are removed.

- From SNPeff-annotated VCF files (one per individual), 
intergenic and intron variants are removed, along with 
the sites that were identified as being fixed for the ALT 
allele across all samples.

- Variants allocated to the following categories are extracted: 
a) Low impact: mostly harmless or unlikely to change protein 
behaviour; b) Moderate impact: non-disruptive variants that
might change protein effectiveness; c) High impact: variants 
assumed to have high (disruptive) impact on protein, probably 
causing  protein truncation, loss of function or triggering 
nonsense-mediated decay and including stop gained codons, 
splice donor variant and splice acceptor, start codon lost; 
d) Synonymous variants.

- Each individual's potential (total) load is calculated by 
taking the sum of the number of variants of each variant impact 
category i, correcting for potential mapping biases by dividing 
the number of each category by the total number of synonymous 
SNPs (Xue et al. 2015).

- Realised load is calculated by dividing the total number 
of homozygous variants of category i by twice the total 
number of variants for category i per individual (Mathur & 
DeWoody 2021).


References:

- Mathur, S. & DeWoody, J. A. Genetic load has potential in 
large populations but is realized in small inbred populations. 
Evol. Appl. 14, 1540–1557 (2021).

- Xue, Y. et al. Mountain gorilla genomes reveal the impact of 
long-term population decline and inbreeding. Science 348, 
242–245 (2015).
