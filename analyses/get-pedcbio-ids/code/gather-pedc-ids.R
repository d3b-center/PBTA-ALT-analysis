library(readr)
library(tidyr)
library(dplyr)

meta<- read_tsv("analyses/add-histologies/output/ALT PBTA oct 2021 (including all plates)-updated-hist-alt.tsv")

meta %>%
  dplyr::group_by(phenotype, group) %>%
  dplyr::select(Kids_First_Biospecimen_ID_DNA)

v21 <- read_tsv("analyses/add-histologies/input-v21/pbta-histologies.tsv")

rna <- v21 %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  pull(Kids_First_Biospecimen_ID)

pedc <- read_tsv("analyses/get-pedcbio-ids/input/ped_opentargets_2021_clinical_data.tsv") %>%
   # remove RNA-only samples
   filter(!SPECIMEN_ID %in% rna) %>%
   # separate DNA and RNA ids
   separate(col = SPECIMEN_ID, 
            into = c("Kids_First_Biospecimen_ID_DNA", "Kids_First_Biospecimen_ID_RNA"), 
            sep = ";") %>%
   mutate(Kids_First_Participant_ID = `Patient ID`)
 
# join with meta
meta_ids <- meta %>%
  dplyr::left_join(pedc, by = c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID_DNA", "Kids_First_Biospecimen_ID_RNA")) %>%
  dplyr::select(Kids_First_Participant_ID, Kids_First_Biospecimen_ID_DNA, Kids_First_Biospecimen_ID_RNA, 
         sample_id, phenotype, `CCA Sept 2021`, `Sample ID`, group, `TH T/TH N`) %>%
  dplyr::filter(!is.na(`Sample ID`))

# write case lists
alt_hgat <- meta_ids %>%
  dplyr::filter(group == "HGAT" & phenotype == "ALT") %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv("analyses/get-pedcbio-ids/output/alt_hgat_ids.txt", col_names = F)

nonalt_hgat <- meta_ids %>%
  dplyr::filter(group == "HGAT" & phenotype == "non-ALT") %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv("analyses/get-pedcbio-ids/output/nonalt_hgat_ids.txt", col_names = F)

alt_all <- meta_ids %>%
  dplyr::filter(phenotype == "ALT") %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv("analyses/get-pedcbio-ids/output/alt_ids.txt", col_names = F)

nonalt_all <- meta_ids %>%
  dplyr::filter(phenotype == "non-ALT") %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv("analyses/get-pedcbio-ids/output/nonalt_ids.txt", col_names = F)

# add additional two groups
hgat_big_telhunt <- meta_ids %>%
  dplyr::filter(group == "HGAT" & `TH T/TH N` > 1.07) %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv("analyses/get-pedcbio-ids/output/hgat_telhunt_over107_ids.txt", col_names = F)

bind_rows(alt_hgat, hgat_big_telhunt) %>% 
  distinct() %>%
  readr::write_tsv("analyses/get-pedcbio-ids/output/alt_hgat_w_hgat_telhunt_over107_ids.txt", col_names = F)

hgat_small_telhunt <- meta_ids %>%
  dplyr::filter(group == "HGAT" & `TH T/TH N` < 1.07) %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv("analyses/get-pedcbio-ids/output/hgat_telhunt_less107_ids.txt", col_names = F)

bind_rows(nonalt_hgat, hgat_small_telhunt) %>%
  distinct() %>% 
  readr::write_tsv("analyses/get-pedcbio-ids/output/alt_hgat_telhunt_less107_ids.txt", col_names = F)

 names(meta)
 