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

# read in input
# histologies files
v11 <- read_tsv(file.path(data_dir, "histologies.tsv"), guess_max = 100000)
v22 <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"), guess_max = 3000)

# 85 hgat used in primary analysis
hgat_subset_ids <- read_tsv(file.path(onco_dir,"hgat_subset.tsv")) %>% 
  pull(sample_id)

# MAFs for PBTA + DGD
# pbta MAF - read in only oncogenic/likely oncogenic variants for ATRX
# source mutations used for oncoprint
source(file.path(onco_dir, "mutation-colors.R"))
mut_of_interest <- c(names(colors), "3'UTR", "5'UTR", "Splice_Region")

pbta_maf <- read_tsv(file.path(anno_maf_dir, "snv-consensus-plus-hotspots-hgat-oncokb.maf.tsv")) %>%
  rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  filter(Hugo_Symbol == "ATRX",
         Variant_Classification %in% mut_of_interest) %>%
  mutate(HGVSp_Short = case_when(is.na(HGVSp_Short) & Variant_Classification == "Splice_Region" ~
           "splice",
           HGVSp_Short = is.na(HGVSp_Short) & Variant_Classification == "3'UTR" ~
           "3'UTR",
         TRUE ~ as.character(HGVSp_Short))) %>%
  mutate(ATRXm = case_when(ONCOGENIC %in% c("Likely Oncogenic", "Oncogenic") ~ 
                           paste0(HGVSp_Short),
                           ONCOGENIC == "Unknown" ~ paste0(HGVSp_Short, " (VUS)")))
pbta_atrx <- pbta_maf %>%
  select(Kids_First_Biospecimen_ID, ATRXm) %>%
  left_join(v11[,c("Kids_First_Biospecimen_ID", "sample_id")]) %>%
  select(sample_id, ATRXm) 

#DGD MAF
dgd_maf <- read_tsv(file.path(data_dir, "snv-dgd.maf.tsv.gz")) %>%
  rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  filter(Hugo_Symbol == "ATRX",
         Variant_Classification %in% mut_of_interest)

dgd_atrx <- dgd_maf %>%
  select(Kids_First_Biospecimen_ID, HGVSp_Short) %>%
  left_join(v11[,c("Kids_First_Biospecimen_ID", "sample_id")]) %>%
  rename(ATRXm = HGVSp_Short) %>%
  select(sample_id, ATRXm) 

# remove samples duplicated in DGD
remove_atrx <- intersect(dgd_atrx$sample_id, pbta_atrx$sample_id)
dgd_atrx <- dgd_atrx %>%
  filter(sample_id != remove_atrx)

atrx_mut <- bind_rows(pbta_atrx, dgd_atrx)

samples_mut_profiled <- v11 %>%
  filter(cohort %in% c("PBTA", "DGD"),
         experimental_strategy %in% c("WGS", "Targeted Sequencing"),
         is.na(RNA_library),
         !is.na(pathology_diagnosis)) %>%
  pull(sample_id) %>%
  unique()

# IHC (TMA) results
ihc <- readxl::read_excel(file.path(onco_dir, "TMA table for HGAT paper_052722_kac.xlsx")) %>%
  rename(sample_id = ID,
         `ATRX IHC` = `ATRX IHC (Pathology)`,
         `UBTF` = `Presence of UBTF`) %>%
  mutate(`ATRX IHC` = case_when(`ATRX IHC` == 0 ~ "RETAINED",
                                `ATRX IHC` == 1 ~ "LOST",
                                TRUE ~ "Not done"),
         `C-Circle Assay` = case_when(`C-Circle Assay` == 0 ~ "NEG",
                                      `C-Circle Assay` == 1 ~ "POS",
                                      TRUE ~ "Not done"),
         UBTF = case_when(UBTF == 1 ~ "POS",
                          UBTF == 0 ~ "NEG",
                                      TRUE ~ "Not done")) %>%
  select(sample_id, `C-Circle Assay`, UBTF, `ATRX IHC`, `Research Subject ID`)

# telhunt results
telhunt <- read_tsv(file.path(tel_dir, "telomere_940_ratio.tsv")) %>%
  rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_tumor,
         `T/N TelHunt ratio` = ratio) %>%
  mutate(`T/N TelHunt ratio` = round(`T/N TelHunt ratio`, 3)) %>%
  select(Kids_First_Biospecimen_ID, `T/N TelHunt ratio`) %>%
  left_join(v22) %>%
  filter(composition != "Derived Cell Line",
         sample_id %in% ihc$sample_id) %>%
  select(sample_id, `T/N TelHunt ratio`) %>%
  unique()

# join everything together
all_pbta_ihc <- v11 %>% 
  filter(!is.na(pathology_diagnosis)) %>%
  select(cohort_participant_id, sample_id, tumor_descriptor) %>%
  unique() %>%
  right_join(ihc) %>%
  left_join(atrx_mut) %>%
  left_join(telhunt) %>%
  mutate(Cohort = case_when(sample_id %in% hgat_subset_ids ~ "Primary Analysis",
                            TRUE ~ "Validation"),
         cohort_participant_id = case_when(is.na(cohort_participant_id) ~ `Research Subject ID`,
                                           TRUE ~ as.character(cohort_participant_id)),
         tumor_descriptor = case_when(is.na(tumor_descriptor) ~ "Not Reported",
                                      TRUE ~ as.character(tumor_descriptor)),
         ATRXm = case_when(is.na(ATRXm) & sample_id %in% samples_mut_profiled ~ "Not mutated",
                           is.na(ATRXm) & !sample_id %in% samples_mut_profiled ~ "Not profiled",
                           TRUE ~ as.character(ATRXm)),
         # filter duplicate rows if VUS + oncogenic
         remove = case_when(sample_id == "7316-2189" & grepl("VUS", ATRXm) ~ "remove",
                            sample_id == "7316-2756" & grepl("VUS", ATRXm) ~ "remove",
                            TRUE ~ "keep")
         )  %>%
  filter(remove == "keep") %>%
  select(-c(`Research Subject ID`, remove)) %>%
  arrange(cohort_participant_id, Cohort) %>%
  filter(!is.na(cohort_participant_id)) %>%
  rename(`Patient ID` = cohort_participant_id,
         `Tumor ID` = sample_id,
         `Phase of Therapy` = tumor_descriptor) 

  # which samples are duplicated - choose oncogenic mutation only?
all_pbta_ihc[duplicated(all_pbta_ihc$`Tumor ID`)|duplicated(all_pbta_ihc$`Tumor ID`, fromLast=TRUE),]

all_pbta_ihc %>%
  openxlsx::write.xlsx(file.path(output_dir, "Table-S1.xlsx"), overwrite = T, keepNA=TRUE)
