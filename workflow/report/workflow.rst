.. image:: https://github.com/NBISweden/GenErode/blob/main/docs/source/img/logga_viridis2.png?raw=true
   :width: 124.0px
   :height: 175.4px
   :alt: GenErode pipeline logo

This is an automatic report from an execution of GenErode, created by Snakemake. 

Each subsection in the Results section to the left contains a brief description of the output and links to some of the output files and the code used to generate them.
If you publish any research using results produced with GenErode, please cite our publication and refer to the source code.

Result sections are available for the following pipeline steps:

- BAM file processing
- mlRho
- PCA (plink)
- Runs of homozygosity (ROH; plink)
- snpEff
- GERP

Please set all steps to "True" in the config.yaml file that you would like to include into the report.

The rulegraph for GenErode is very large. If you wish to zoom into the image and inspect how the different rules of the pipeline are connected to each other, please save it as SVG or PNG (via the icon in the upper right corner of the rulegraph).