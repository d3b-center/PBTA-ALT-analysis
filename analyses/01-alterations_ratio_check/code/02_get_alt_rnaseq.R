library(readr)
library(tidyverse)
library(reshape2)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

data_dir <- file.path(root_dir, "data")


# read in rnaseq data
expressionMatrixPolya<-readRDS(file.path(data_dir,"pbta-gene-expression-rsem-fpkm-collapsed.polya.rds"))
expressionMatrixStranded<-readRDS(file.path(data_dir,"pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"))

# read in histology
histology <-read_tsv(file.path(data_dir, "pbta-histologies.tsv"))

# get putative oncogene fusions 
putative_driver_fusions <- read_tsv(file.path(data_dir, "pbta-fusion-putative-oncogenic.tsv"))

# filter for ATRX,DAXX and TERT
fus_goi <- putative_driver_fusions %>% dplyr::filter(grepl("ATRX|TERT|DAXX",FusionName)) %>%
  dplyr::select(Sample,FusionName) %>%
  unique() %>%
  group_by(Sample) %>%
  summarise(FusionName=toString(FusionName))

# reshape fpkm values
exp_stranded <- expressionMatrixStranded[c("ATRX","DAXX", "TERT"),]
exp_stranded <- exp_stranded %>%
  t() %>% as.data.frame() %>%
  dplyr::rename(ATRX_fpkm=ATRX, TERT_fpkm=TERT,  DAXX_fpkm=DAXX)

exp_polya <- expressionMatrixPolya[c("ATRX","DAXX", "TERT"),]
exp_polya <- exp_polya %>%
  t() %>% as.data.frame() %>%
  dplyr::rename(ATRX_fpkm=ATRX, TERT_fpkm=TERT,  DAXX_fpkm=DAXX)

exp_goi <- rbind(exp_stranded,exp_polya) %>% rownames_to_column()

rna_alt <- exp_goi %>% left_join(fus_goi,by=c("rowname"="Sample"))  %>%
  dplyr::rename(Kids_First_Biospecimen_ID=rowname) %>%
  dplyr::mutate(ATRX_fusion = case_when(grepl("ATRX",FusionName) ~ FusionName,
                                        TRUE ~ NA_character_ ),
                DAXX_fusion = case_when(grepl("DAXX",FusionName) ~ FusionName,
                                        TRUE ~ NA_character_ ),
                TERT_fusion = case_when(grepl("TERT",FusionName) ~ FusionName,
                                        TRUE ~ NA_character_ )) %>%
  # remove column name FusionName
  dplyr::select(-FusionName)

histology %>%
  dplyr::filter(experimental_strategy %in% c("RNA-Seq")) %>%
  dplyr::select("Kids_First_Biospecimen_ID","Kids_First_Participant_ID","sample_id") %>%
  left_join(rna_alt) %>% 
  write_tsv(file.path(root_dir,"analyses/01-alterations_ratio_check","output","PutativeDriver_ATRTX_DAXX_TERT_RNA_alt.tsv"))

