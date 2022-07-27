# load libraries
library(tidyverse)
library(openxlsx)

# define directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
onco_dir <- file.path(root_dir, "analyses", "05-oncoplot", "input")
tel_dir <- file.path(root_dir, "analyses", "01-alterations_ratio_check", "input")
output_dir <- file.path(root_dir, "analyses", "11-tables", "output")
data_dir <- file.path(root_dir, "data")

# make output dir
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

##read in input
hgat_subset_ids <- read_tsv(file.path(onco_dir,"hgat_subset.tsv")) %>% 
  pull(sample_id)

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

v22 <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"), guess_max = 3000)

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

dups <- telhunt[duplicated(telhunt$sample_id),]


all_pbta_ihc <- v11 %>% 
  filter(!is.na(pathology_diagnosis)) %>%
  select(cohort_participant_id, sample_id, tumor_descriptor) %>%
  unique() %>%
  right_join(ihc) %>%
  left_join(telhunt) %>%
  mutate(Cohort = case_when(sample_id %in% hgat_subset_ids ~ "Primary Analysis",
                            TRUE ~ "Validation"),
         cohort_participant_id = case_when(is.na(cohort_participant_id) ~ `Research Subject ID`,
                                           TRUE ~ as.character(cohort_participant_id)),
         tumor_descriptor = case_when(is.na(tumor_descriptor) ~ "Not Reported",
                                      TRUE ~ as.character(tumor_descriptor))
         )  %>%
  select(-`Research Subject ID`) %>%
  arrange(cohort_participant_id, Cohort) %>%
  rename(`Patient ID` = cohort_participant_id,
         `Tumor ID` = sample_id,
         `Phase of Therapy` = tumor_descriptor) %>%
  openxlsx::write.xlsx(file.path(output_dir, "Table-S1.xlsx"), overwrite = T, keepNA=TRUE)
