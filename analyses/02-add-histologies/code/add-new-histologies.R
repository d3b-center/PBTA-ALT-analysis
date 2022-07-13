library(tidyverse)

# define directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "02-add-histologies")
data_dir <- file.path(root_dir, "data") 
  
# read in v22 histologies
hist <- read_tsv(file.path(data_dir,
                           "pbta-histologies.tsv"))

# read in file from jenny
hgat <- readxl::read_excel(file.path(analysis_dir, 
                                    "input-jenny",
                                    "Stundon_HGAT_One_Sample_Per_Pt_adults_removed.xlsx"),
                          .name_repair = "unique") %>%
  mutate(age_last_update_days = as.numeric(age_last_update_days))

nonhgat <- readxl::read_excel(file.path(analysis_dir, 
                                         "input-jenny",
                                         "non_HGAT_4.21.22.xlsx"),
                               .name_repair = "unique") %>%
  dplyr::rename(`alt final` = `FINAL ALT`)

full_set <- hgat %>%
  bind_rows(nonhgat)

# read in alterations file
alt <- read_tsv(file.path(root_dir, 
                          "analyses",
                          "01-alterations_ratio_check",
                          "output",
                          "PutativeDriver_ATRX_DAXX_TERT_DNA_alt.tsv"), guess_max = 10000)

# remove columns from jenny's while which are going to be replaced by histologies file columns
rm.cols <- intersect(names(jen), names(hist))
jen.new <- jen[,!colnames(jen) %in% rm.cols]

#grab only appropriate rows from hist based on DNA sample
hist_subset <- hist %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% full_set$Kids_First_Biospecimen_ID_DNA) %>%
  dplyr::select(c(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, sample_id, all_of(rm.cols)))

jen_mer <- jen.new %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_DNA) %>%
  dplyr::left_join(hist_subset, by = c("Kids_First_Biospecimen_ID")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) %>%
  distinct() 

# update phenotype
jen_mer$phenotype <- ifelse(jen_mer$`CCA Sept 2021` == "POS", "ALT", 
                            ifelse(jen_mer$`CCA Sept 2021` == "NEG", "non-ALT", NA))
# add grouping
jen_mer$group <- ifelse(jen_mer$short_histology == "HGAT", "HGAT", "non-HGAT")

# calculate tel ratio
jen_mer$telomere_ratio <- jen_mer$tel_content_tumor/jen_mer$tel_content_normal

# update analysis columns
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
  distinct() %>%
  arrange(Kids_First_Biospecimen_ID_DNA) %>%
  write_tsv(file.path(analysis_dir,
                      "output/stundon_hgat_updated_hist_alt.tsv"))
  

