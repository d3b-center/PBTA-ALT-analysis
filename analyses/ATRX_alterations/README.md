## Description

In this analysis we find distribution of telmerehunter ratio with respect to ATRX alterations. Using PBTA v17 release files we identify SNV alterations in ATRX (no CNV ATRX alt were found) from strelka2,mutect2 or concensus method and plot boxplot for ratio distributions.

## Input
data/pbta-snv-mutect2.vep.maf.gz
data/pbta-snv-strelka2.vep.maf.gz
data/pbta-snv-consensus-mutation.maf.tsv.gz
data/consensus_seg_annotated_cn_autosomes.tsv.gz

## Ouptut
plots/ATRX_telomere_ratio_boxplot.pdf

## Code
01_get_alt_cohort3a.Rmd	: reads in strelka2,mutect2 and concensus files and looks for ATRX mutation and plots boxplots

