library(readr)
library(tidyr)
library(dplyr)

# define directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "03-get-pedcbio-ids")
input_dir <- file.path(root_dir, "analyses", "02-add-histologies")
output_dir <- file.path(analysis_dir, "output")
data_dir <- file.path(root_dir, "data")
  
meta<- read_tsv(file.path(input_dir, 
                         "output",
                         "stundon_hgat_updated_hist_alt.tsv"))

meta %>%
  dplyr::group_by(`alt final`, group) %>%
  dplyr::select(Kids_First_Biospecimen_ID_DNA)

hist <- read_tsv(file.path(data_dir,
                          "pbta-histologies.tsv"))
rna <- hist %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  pull(Kids_First_Biospecimen_ID)

pedc <- read_tsv(file.path(analysis_dir, 
                           "input",
                           "ped_opentargets_2021_clinical_data.tsv")) %>%
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
         sample_id, `alt final`, `CCA Sept 2021`, `Sample ID`, group, `TH T/TH N`) %>%
  dplyr::filter(!is.na(`Sample ID`))

# write case lists
alt_hgat <- meta_ids %>%
  dplyr::filter(group == "HGAT" & `alt final` == "POS") %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv(file.path(output_dir,
                             "alt_hgat_ids.txt"), col_names = F)

nonalt_hgat <- meta_ids %>%
  dplyr::filter(group == "HGAT" & `alt final` == "NEG") %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv(file.path(output_dir,
                             "nonalt_hgat_ids.txt"), col_names = F)

alt_all <- meta_ids %>%
  dplyr::filter(`alt final` == "POS") %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv(file.path(output_dir,"alt_ids.txt"), col_names = F)

nonalt_all <- meta_ids %>%
  dplyr::filter(`alt final` == "NEG") %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv(file.path(output_dir, "nonalt_ids.txt"), col_names = F)

# add additional two groups
hgat_big_telhunt <- meta_ids %>%
  dplyr::filter(group == "HGAT" & `TH T/TH N` > 1.07) %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv(file.path(output_dir, "hgat_telhunt_over107_ids.txt"), col_names = F)

bind_rows(alt_hgat, hgat_big_telhunt) %>% 
  distinct() %>%
  readr::write_tsv(file.path(output_dir, "alt_hgat_w_hgat_telhunt_over107_ids.txt"), col_names = F)

hgat_small_telhunt <- meta_ids %>%
  dplyr::filter(group == "HGAT" & `TH T/TH N` < 1.07) %>%
  dplyr::select(`Sample ID`) %>%
  readr::write_tsv(file.path(output_dir, "hgat_telhunt_less107_ids.txt"), col_names = F)

bind_rows(nonalt_hgat, hgat_small_telhunt) %>%
  distinct() %>% 
  readr::write_tsv(file.path(output_dir, "nonalt_hgat_w_hgat_telhunt_less107_ids.txt"), col_names = F)