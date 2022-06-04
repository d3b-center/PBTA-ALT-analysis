library(tidyverse)

# define directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "02-add-histologies")

# read in v21 histologies
hist <- read_tsv(file.path(analysis_dir, 
                           "input-v21",
                           "pbta-histologies.tsv"))

# read in file from jenny
jen_hgat <- readxl::read_excel(file.path(analysis_dir, 
                                    "input-jenny",
                                    "Stundon_HGAT_One_Sample_Per_Pt_adults_removed.xlsx"),
                          .name_repair = "unique") %>%
  mutate(age_last_update_days = as.numeric(age_last_update_days))

jen_nonhgat <- readxl::read_excel(file.path(analysis_dir, 
                                         "input-jenny",
                                         "non_HGAT_4.21.22.xlsx"),
                               .name_repair = "unique") %>%
  dplyr::rename(`alt final` = `FINAL ALT`)

jen <- jen_hgat %>%
  bind_rows(jen_nonhgat)

# read in alterations file
alt <- read_tsv(file.path(root_dir, 
                          "analyses",
                          "01-alterations_ratio_check",
                          "output",
                          "PutativeDriver_ATRX_DAXX_TERT_DNA_alt.tsv"), guess_max = 10000)

# remove columns from jenny's while which are going to be replaced by histologies file columns
rm.cols <- intersect(names(jen), names(hist))
jen.new <- jen[,!colnames(jen) %in% rm.cols]

#grab only appropriate rows from hist
hist_subset <- hist %>%
  dplyr::filter(sample_type == "Tumor" & experimental_strategy != "RNA-Seq") %>%
  dplyr::select(c(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, sample_id, all_of(rm.cols)))

jen_mer <- jen.new %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_DNA) %>%
  dplyr::left_join(hist_subset, by = c("Kids_First_Biospecimen_ID")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) %>%
  distinct() 
jen_mer$phenotype <- ifelse(jen_mer$`CCA Sept 2021` == "POS", "ALT", 
                            ifelse(jen_mer$`CCA Sept 2021` == "NEG", "non-ALT", NA))
jen_mer %>%
  select(phenotype, `CCA Sept 2021`)
jen_mer$group <- ifelse(jen_mer$short_histology == "HGAT", "HGAT", "non-HGAT")
jen_mer$telomere_ratio <- jen_mer$tel_content_tumor/jen_mer$tel_content_normal

# update tel analysis
names(alt)
rm_alt_cols <- names(alt[,c(4:ncol(alt))])

# columns which are not in alt
jen_mer_rm <- jen_mer[,!colnames(jen_mer) %in% rm_alt_cols]
names(jen_mer_rm)

# merge back cols
jen_mer_alt <- jen_mer_rm %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_DNA) %>%
  dplyr::left_join(alt, by = c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID","sample_id")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) %>%
  distinct() 

# take the base for parent aliquot id to only keep DNA-RNA match with the same parent aliquot ID
# get base parent aliquot ID for RNA samples
rna_aliquot_df <- hist %>% 
  dplyr::filter(sample_type == "Tumor" & experimental_strategy == "RNA-Seq") %>% 
  dplyr::select(Kids_First_Biospecimen_ID, sample_id, parent_aliquot_id) %>% 
  dplyr::mutate(match_aliquot_id = gsub("\\..*", "", parent_aliquot_id)) %>% 
  dplyr::select(-parent_aliquot_id) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_RNA=Kids_First_Biospecimen_ID)
  
# only keep entries where DNA-RNA match have the same parent aliquot ID or is NA
jen_mer_alt_match <- jen_mer_alt %>% 
  dplyr::select(-Kids_First_Biospecimen_ID_RNA) %>% 
  dplyr::mutate(match_aliquot_id = gsub("\\..*", "", parent_aliquot_id)) %>% 
  dplyr::left_join(rna_aliquot_df) %>%
  dplyr::select(-match_aliquot_id)

# additionally, for each sample ID we need to take independent DNA sample
# implement independent sample list module in OpenPedCan
set.seed(2021)
# first sample them
jen_mer_alt_match <- jen_mer_alt_match[sample(nrow(jen_mer_alt_match)), ]
# then take distinct
jen_mer_alt_match <- jen_mer_alt_match %>%
  dplyr::distinct(sample_id, .keep_all = TRUE)

jen_mer_alt_match %>%
  dplyr::select(Kids_First_Biospecimen_ID_DNA, Kids_First_Biospecimen_ID_RNA, colnames(jen_mer_alt_match)[2:114]) %>%
  arrange(Kids_First_Biospecimen_ID_DNA) %>%
  write_tsv(file.path(analysis_dir,
                      "output/stundon_hgat_updated_hist_alt.tsv"))


