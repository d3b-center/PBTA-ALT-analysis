library(tidyverse)

# read in v21 histologies
hist <- read_tsv("analyses/add-histologies/input-v21/pbta-histologies.tsv")

# read in file from jenny
jen <- readxl::read_excel("analyses/add-histologies/input-jenny/ALT PBTA oct 2021 (including all plates).xlsx", 
                          .name_repair = "unique") %>%
  rename(Sample_id = sample_id)

# read in alterations file
alt <- read_tsv("analyses/alterations_ratio_check/output/PutativeDriver_ATRX_DAXX_TERT_DNA_alt.tsv", guess_max = 10000)

# remove columns from jenny's while which are going to be replaced by histologies file columns
rm.cols <- intersect(names(jen), names(hist))
jen.new <- jen[,!colnames(jen) %in% rm.cols]

#grab only appropriate rows from hist
hist_subset <- hist %>%
  filter(sample_type == "Tumor" & experimental_strategy != "RNA-Seq") %>%
  select(c(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, sample_id, all_of(rm.cols)))

jen_mer <- jen.new %>%
  rename(sample_id = Sample_id,
        # Kids_First_Participant_ID = pt_id,
         Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_DNA) %>%
  left_join(hist_subset, by = c("Kids_First_Biospecimen_ID", "sample_id")) %>%
  rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) %>%
  distinct() 
jen_mer$phenotype <- ifelse(jen_mer$`CCA Sept 2021` == "POS", "ALT", "non-ALT")
jen_mer$group <- ifelse(jen_mer$short_histology == "HGAT", "HGAT", "non-HGAT")
jen_mer$telomere_ratio <- jen_mer$tel_content_tumor/jen_mer$tel_content_normal

# update tel analysis
rm_alt_cols <- intersect(names(jen_mer), names(alt))
rm_alt_cols <- rm_alt_cols[-c(1,29)]
rm_alt_cols

# columns which are not in alt
jen_mer_rm <- jen_mer[,!colnames(jen_mer) %in% rm_alt_cols]
names(jen_mer_rm)
names(jen_mer)
# merge back cols
jen_mer_alt <- jen_mer_rm %>%
  rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_DNA) %>%
  left_join(alt, by = c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID","sample_id")) %>%
  rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) %>%
  distinct() 

jen_mer_alt %>%
  write_tsv("analyses/add-histologies/output/ALT PBTA oct 2021 (including all plates)-updated-hist-alt.tsv")

# which are duplicate sample_ids?
dups <- jen_mer_alt[duplicated(jen_mer_alt$sample_id),]
dups <- dups %>%
  select(Kids_First_Biospecimen_ID_DNA, Kids_First_Biospecimen_ID_RNA, sample_id, 
         tumor_descriptor, sample_type, composition) %>%
  arrange(sample_id)

