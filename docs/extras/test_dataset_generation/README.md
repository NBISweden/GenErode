# Description of test dataset generation

We have generated a test dataset based on Sumatran rhinoceros 
(_Dicerorhinus sumatrensis_) whole-genome re-sequencing data 
(from von Seth et al. 2021) that we reduced in size so that 
users have the possibility to get familiar with the pipeline 
before analyzing their own genome-wide datasets. We published 
the test dataset along with the GenErode pipeline in the 
Scilifelab data repository (DOI: pending). Results presented 
in the pipeline article (Kutschera et al. 2022) were obtained 
from analyzing the test dataset with GenErode.


## Brief summary

The test dataset contains three references in FASTA format. 
We extracted scaffold ‘Sc9M7eS_2_HRSCAF_41’ of size 40,842,778 bp 
from the Sumatran rhinoceros genome assembly (_Dicerorhinus sumatrensis_ 
_harrissoni_; GenBank accession number GCA_014189135.1) to be 
used as reference genome for PCAs, mlRho and ROH analyses in 
GenErode. Some of the GenErode steps require the reference genome 
of a closely related species. We therefore identified scaffolds 
from the White rhinoceros genome assembly (_Ceratotherium simum_
_simum_; GenBank accession number GCF_000283155.1) as being 
putatively orthologs to the Sumatran rhinoceros scaffold through 
a reciprocal blast and a mapping approach, to be used as reference 
genome for PCAs, snpEff and GERP analyses. The three putatively 
orthologous White rhinoceros scaffolds ‘NW_004454182.1’, 
‘NW_004454248.1’, and ‘NW_004454260.1’ have a combined length 
of 41,195,616 bp. For snpEff analyses, we additionally provide 
gene predictions from the three White rhinoceros scaffolds in 
GTF format. The repository also contains a Sumatran rhinoceros 
mitochondrial genome (GenBank accession number NC_012684.1) to 
be used as reference for the optional mitochondrial mapping step 
in GenErode. A detailed workflow description is provided in 
`workflows/references.md` and the corresponding scripts can be
found in `scripts/references`.

The whole-genome re-sequencing data was processed as follows to 
be included into the test dataset. We extracted a subset of reads 
from whole-genome re-sequencing data from three historical and 
three modern Sumatran rhinoceros samples from the now-extinct 
Malay Peninsula population. These samples had been sequenced and 
analyzed for patterns of genome erosion by von Seth et al. (2021) 
(SRA identifiers: ERS4044060, ERS4044061, ERS4044063, ERS4042484, 
ERS4042485, ERS4042486). For two of the historical Sumatran 
rhinoceros samples, three sequencing libraries are available 
per sample that had been sequenced on two lanes each. For the 
third historical sample, 24 sequencing libraries are available 
that had been sequenced on two lanes (12 libraries per lane). 
For each of the three modern Sumatran rhinoceros samples, one 
sequencing library is available. Paired-end reads from each 
sequencing library were included into the test dataset that mapped 
to the Sumatran rhinoceros scaffold ‘Sc9M7eS_2_HRSCAF_41’, along 
with a small proportion of randomly selected reads that mapped 
to the Sumatran rhinoceros mitochondrial genome or elsewhere in 
the genome.  A detailed workflow description is provided in 
`workflows/resequencing_data.md` and the corresponding scripts can be
found in `scripts/resequencing_data`.

For GERP analyses, scaffolds from the genome assemblies of 30 
mammalian outgroup species were identified as putative orthologs 
to the Sumatran rhinoceros scaffold using reciprocal blast, and 
are provided in FASTA format (gzipped). Further, a phylogeny of 
the White rhinoceros and the 30 outgroup species including 
divergence time estimates (in billions of years) from 
timetree.org is provided in NEWICK format.  A detailed workflow 
description is provided in `workflows/gerp_outgroups.md` and the 
corresponding scripts can be found in `scripts/gerp_outgroups`.

Finally, the repository contains configuration and metadata 
files that were used for three separate runs of GenErode to 
generate the results presented in Kutschera et al. (2022).

Further details on the generation of this test dataset are 
provided in Kutschera et al. (2022) and in the README.md
files in each respective folder.

### References

von Seth J, Dussex N, Díez-Del-Molino D, van der Valk T, 
Kutschera VE, Kierczak M, et al. Genomic insights into the 
conservation status of the world’s last remaining Sumatran 
rhinoceros populations. Nat Commun. 2021;12:2393.

Kutschera VE, Kierczak M, van der Valk T, von Seth J, Dussex N, 
Lord E, et al. GenErode: a bioinformatics pipeline to investigate 
genome erosion in endangered and extinct species. bioRxiv. 2022.

