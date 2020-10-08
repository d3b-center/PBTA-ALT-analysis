## Description

In this analysis we find distribution of telmerehunter ratio with respect to ATRX alterations. Using PBTA v17 release files we identify SNV alterations in ATRX|DAXX|TERT from strela2,mutect2 and concensus snv method . CNV alterations in ATRX|DAXX|TERT from consensus CNV method were also identified. Calls from manta are also added in files but not used in `defining_alt` columns since almost all samples have SVs and might need more filtering/review.

TelomereHunter was run for telomre content in tumor and normal. We plot log2(ratio) in the plots for clear clarification.
TelomeraseScores were derived from [here](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/telomerase-activity-prediction)


## Input v17
 - data/pbta-snv-mutect2.vep.maf.gz
 - data/pbta-snv-strelka2.vep.maf.gz
 - data/pbta-snv-consensus-mutation.maf.tsv.gz
 - data/consensus_seg_annotated_cn_autosomes.tsv.gz data/consensus_seg_annotated_cn_x_and_y.tsv.gz
 - data/cnvkit_annotated_cn_autosomes.tsv.gz data/cnvkit_annotated_cn_x_and_y.tsv.gz
 - data/controlfreec_annotated_cn_autosomes.tsv.gz data/controlfreec_annotated_cn_x_and_y.tsv.gz
 - data/pbta-sv-manta.tsv.gz
 - data/pbta-gene-expression-rsem-fpkm.polya.rds
 - data/pbta-gene-expression-rsem-fpkm.stranded.rds
 - data/pbta-fusion-arriba.tsv.gz
 - data/pbta-fusion-starfusion.tsv.gz

## Input telhunter and extend scores
 - TelomeraseScores_PTBAStranded_FPKM.txt
 - TelomeraseScores_PTBAPolya_FPKM.txt
 - telomere_940_ratio.tsv

## Ouptut

plots/
├── all_alterations_ATRX.png # ATRX boxplot for samples with alterations `_mut` and without `_WT`
├── all_alterations_ATRX_DAXX_TERT.png # ATRX,DAXX and TERT altrations boxplot for samples with alterations `_mut` and without `_WT`
├── all_alterations_ATRX_or_DAXX.png # ATRX or DAXX boxplot for samples with alterations `_mut` and without `_WT`
├── all_alterations_DAXX.png # DAXX boxplot for samples with alterations `_mut` and without `_WT`
├── all_alterations_TERT.png # TERT boxplot for samples with alterations `_mut` and without `_WT`
└── corr_score_fpkm.png


output/
├── ATRTX_DAXX_TERT_FPKM_expression.tsv # ATRX,DAXX and TERT expression from stranded and polya 
├── FilteredFusionAnnoFuse.tsv # All PBTA filtered fusion in PBTA from annoFuse v0.90.0
├── PutativeDriver_ATRTX_DAXX_TERT_AnnoFuse.tsv # ATRX TERT and DAXX driver fusion in PBTA from annoFuse v0.90.0
├── PutativeDriver_ATRTX_DAXX_TERT_RNA_alt.tsv # ATRX TERT and DAXX merged rna fpkm and driver fusion
├── PutativeDriver_ATRX_DNA_alt.tsv # ATRX SNV, CNV and SV alterations
├── PutativeDriver_DAXX_DNA_alt.tsv # DAXX SNV, CNV and SV alterations
├── PutativeDriver_TERT_DNA_alt.tsv # TERT SNV, CNV and SV alterations
├── all_alterations_ATRX_DAXX_TERT.tsv # merge all alterations in ATRX,DAXX and TERT
├── pbta_telomere_ratio_merged_ATRX.tsv # merged alterations in ATRX with telomere ratio 
├── pbta_telomere_ratio_merged_DAXX.tsv # merged alterations in DAXX with telomere ratio
└── pbta_telomere_ratio_merged_TERT.tsv # merged alterations in TERT with telomere ratio




## Code
01_get_alt_cohort3a.Rmd	: reads in strelka2,mutect2 and concensus files and looks for ATRX , DAXX and TERT DNA mutations
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


04_score_per_alterations.Rmd : correlations within telomere ratio and extend scores as well as FPKM values



