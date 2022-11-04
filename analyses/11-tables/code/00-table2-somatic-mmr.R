# load libraries
library(tidyverse)
library(openxlsx)

# define directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
onco_dir <- file.path(root_dir, "analyses", "05-oncoplot", "input")
anno_maf_dir <- file.path(root_dir, "analyses", "09-lollipop", "output")
tel_dir <- file.path(root_dir, "analyses", "01-alterations_ratio_check", "input")
output_dir <- file.path(root_dir, "analyses", "11-tables", "output")
data_dir <- file.path(root_dir, "data")

# make output dir
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# histologies files
v11 <- read_tsv(file.path(data_dir, "histologies.tsv"), guess_max = 100000)

source(file.path(onco_dir, "mutation-colors.R"))
mut_of_interest <- c(names(colors), "3'UTR", "5'UTR", "Splice_Region")

goi <- read_table(file.path(onco_dir, "KEGG_MISMATCH_REPAIR.txt"),
                  col_names = "genes")

table2_ids <- c("7316-2594",
                "7316-515",
                "7316-2561",
                "7316-2308",
                "7316-2640",
                "7316-3225",
                "7316-2189",
                "7316-4215",
                "7316-4917")

pbta_maf <- read_tsv(file.path(anno_maf_dir, "snv-consensus-plus-hotspots-hgat-oncokb.maf.tsv")) %>%
  rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  filter(Hugo_Symbol %in% goi$genes,
         Variant_Classification %in% mut_of_interest,
         ONCOGENIC %in% c("Likely Oncogenic", "Oncogenic")) %>%
  mutate(HGVSp_Short = case_when(is.na(HGVSp_Short) & Variant_Classification == "Splice_Region" ~ "splice",
                                 TRUE ~ as.character(HGVSp_Short)),
         MMR = paste(Hugo_Symbol, HGVSp_Short, sep = " "))

pbta_mut <- pbta_maf %>%
  select(Kids_First_Biospecimen_ID, MMR) %>%
  left_join(v11[,c("Kids_First_Biospecimen_ID", "sample_id")]) %>%
  select(sample_id, MMR) 

#DGD MAF
dgd_maf <- read_tsv(file.path(anno_maf_dir, "dgd_maf-goi-oncokb.tsv")) %>%
  rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  filter(Hugo_Symbol %in% goi$genes,
         Variant_Classification %in% mut_of_interest,
         ONCOGENIC %in% c("Likely Oncogenic", "Oncogenic")) %>%
  mutate(HGVSp_Short = case_when(is.na(HGVSp_Short) & Variant_Classification == "Splice_Region" ~ "splice",
                                 TRUE ~ as.character(HGVSp_Short)),
         MMR = paste(Hugo_Symbol, HGVSp_Short, sep = " ")) %>%
  unique()


dgd_mut <- dgd_maf %>%
  select(Kids_First_Biospecimen_ID, MMR) %>%
  left_join(v11[,c("Kids_First_Biospecimen_ID", "sample_id")]) %>%
  #filter out known germline variants
  filter(!(sample_id == "7316-4215" & MMR == "MSH6 p.K1288_F1289ins*"),
         !(sample_id == "7316-4917" & MMR == "MSH6 p.R497*")) %>%
  select(sample_id, MMR) %>%
  unique()


all_mut <- pbta_mut %>%
  bind_rows(dgd_mut) %>%
  filter(!is.na(MMR),
         sample_id %in% table2_ids) %>%
  unique() %>%
  group_by(sample_id) %>%
  summarise(`Somatic MMR` = str_c(unique(MMR), collapse = ", "))

all_mut %>%
  write_tsv(file.path(output_dir, "Table2-somatic-MMR-oncogenic-mutations.tsv"))
  