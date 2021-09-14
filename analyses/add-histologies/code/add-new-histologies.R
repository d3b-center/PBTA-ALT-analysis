library(tidyverse)
v19 <- read_tsv("~/Documents/GitHub/OpenPBTA-analysis/data/release-v19-20210423/pbta-histologies.tsv")

jen <- readxl::read_excel("~/Box Sync/D3B-share/collaborations/Kristina/from-jenny/ALT May 2021 JS.xlsx", 
                          .name_repair = "unique") %>%
  rename(Sample_id = sample_id)
kc <- readxl::read_excel("~/Box Sync/D3B-share/collaborations/Kristina/from-jenny/ALT May 2021 KC.xlsx", 
                         .name_repair = "unique") %>%
  rename(Sample_id = sample_id)


jen %>%
  select(Sample_id, sample)


#kc %>%
 # select(Sample_id, sample)

rm.cols <- intersect(names(jen), names(v19))
jen.new <- jen[,!colnames(jen) %in% rm.cols]

rm.cols.kc <- intersect(names(kc), names(v19))
kc.new <- kc[,!colnames(kc) %in% rm.cols.kc]


#grab only appropriate columns from v19
v19_subset <- v19 %>%
  filter(sample_type == "Tumor" & experimental_strategy != "RNA-Seq") %>%
  select(c(Kids_First_Participant_ID, Kids_First_Biospecimen_ID, sample_id, rm.cols))

jen_mer <- jen.new %>%
  rename(sample_id = Sample_id,
        # Kids_First_Participant_ID = pt_id,
         Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_DNA) %>%
  left_join(v19_subset, by = c("Kids_First_Biospecimen_ID", "sample_id")) %>%
  rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) %>%
  distinct() 
jen_mer$phenotype <- ifelse(jen_mer$`CCA Binary` == 1, "ALT", "non-ALT")
jen_mer$group <- ifelse(jen_mer$short_histology == "HGAT", "HGAT", "non-HGAT")
jen_mer$telomere_ratio <- jen_mer$tel_content_tumor/jen_mer$tel_content_normal

jen_mer %>%
  write_tsv("~/Box Sync/D3B-share/collaborations/Kristina/from-jenny/ALT_May_2021_JS_plus_v19_histologies.tsv")

kc_mer <- kc.new %>%
  rename(sample_id = Sample_id,
         # Kids_First_Participant_ID = pt_id,
         Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_DNA) %>%
  left_join(v19_subset, by = c("Kids_First_Biospecimen_ID", "sample_id")) %>%
  rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID) %>%
  distinct() 
kc_mer$telomere_ratio <- kc_mer$tel_content_tumor/jen_mer$tel_content_normal

kc_mer %>%
  write_tsv("~/Box Sync/D3B-share/collaborations/Kristina/from-jenny/ALT_May_2021_KC_plus_v19_histologies.tsv")
  
hgat <- jen_mer %>%
  filter(group == "HGAT")
non_hgat <- jen_mer %>%
  filter(group == "non-HGAT")