## Description

In this analysis we find distribution of TelmereHunter ratio with respect to ATRX alterations. 

 - D3b SNV consensus SNV (2 out of 4 callers) alterations in ATRX|DAXX|TERT
 - CNV alterations in ATRX|DAXX|TERT from consensus CNV method were also identified. 
 - Calls from manta are also added in files but not used in `defining_alt` columns yet since manta is included in CNV consensus call mentioned above.

TelomereHunter was run for telomre content in tumor and normal. We plot log2(ratio) in the plots for clarification.


## Input OpenPBTA v22, OpenPedCan v10
- consensus_wgs_plus_cnvkit_wxs.tsv.gz (OpenPedCan)
- pbta-fusion-putative-oncogenic.tsv (OpenPBTA)
- pbta-gene-expression-rsem-fpkm-collapsed.polya.rds (OpenPBTA)
- pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds (OpenPBTA)
- pbta-histologies.tsv (OpenPBTA)
- pbta-sv-manta.tsv.gz (OpenPBTA)
- snv-consensus-plus-hotspots.maf.tsv.gz (OpenPedCan)

## Input telhunter scores
 - telomere_940_ratio.tsv

## Ouptut

```
plots/
├── all_alterations_ATRX.png # ATRX boxplot for samples with alterations `_mut` and without `_WT`
├── all_alterations_ATRX_DAXX_TERT.png # ATRX,DAXX and TERT altrations boxplot for samples with alterations `_mut` and without `_WT`
├── all_alterations_ATRX_or_DAXX.png # ATRX or DAXX boxplot for samples with alterations `_mut` and without `_WT`
├── all_alterations_DAXX.png # DAXX boxplot for samples with alterations `_mut` and without `_WT`
├── all_alterations_TERT.png # TERT boxplot for samples with alterations `_mut` and without `_WT`
└── corr_extend_scores_telhunter_ratio.png # Telomrase score and telomere hunter log2(ratio) scatter plot
```

```
output/
├── PutativeDriver_ATRTX_DAXX_TERT_RNA_alt.tsv # ATRX TERT and DAXX merged rna fpkm and driver fusion
├── all_alterations_ATRX_DAXX_TERT.tsv # merge all alterations in ATRX,DAXX and TERT
├── colnames_all_alterations_ATRX_DAXX_TERT.tsv	# file describes column names in the output file 
└── PutativeDriver_ATRTX_DAXX_TERT_DNA_alt.tsv # ATRX TERT and DAXX merged SNV and CNV alterations 
	
```


## Code
01_get_alt_cohort3a.Rmd	: reads in SNV and CNV concensus files and looks for ATRX , DAXX and TERT DNA mutations
02_get_alt_rnaseq.R : reads in fpkm files and fusion calls to idenitify fusion in ATRX, DAXX and TERT as well as gather fpkm values			
03_merge_telomere_content.Rmd :	Merges telomere hunter ratio, extend scores and alterations status 	


Alterations in ATRX
```
    !is.na(ATRX_MOD_HIGH_mut_consensus) | !is.na(ATRX_loss_consensus) ~ "ATRX_loss",
           !is.na(ATRX_gain_consensus)  ~ "ATRX_gain",
    TRUE ~ NA_character_
```

Alterations in DAXX
```
    !is.na(DAXX_loss_consensus) ~ "DAXX_loss",
       !is.na(DAXX_gain_consensus) ~ "DAXX_gain",
      TRUE ~ NA_character_
```


Alterations in TERT
```
    !is.na(TERT_MOD_HIGH_mut_consensus) | !is.na(TERT_upstream_mut_consensus) |
      !is.na(TERT_amplification_consensus) | !is.na(TERT_gain_consensus) |
      !is.na(TERT_loss_consensus) ~ "TERT_mod",
    TRUE ~ NA_character_

```



