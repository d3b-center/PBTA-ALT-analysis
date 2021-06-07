# Package names
cran_packages <- c("tidyverse","R.utils","rprojroot","devtools")

# Install packages not yet installed
installed_packages <- cran_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(cran_packages[!installed_packages],repos = "http://cran.us.r-project.org")
}

devtools::install_github("jokergoo/ComplexHeatmap")

# Packages loading
invisible(lapply(c(cran_packages,"ComplexHeatmap"), library, character.only = TRUE))

# Setup dir
input_dir <- "input"
data_dir <- "../../data"

##read in input
source(file.path(input_dir, "mutation-colors.R"))
goi.list <- read_tsv(file.path(input_dir, "goi-mutations"), col_names = "genes")

# https://chopri.box.com/s/pjfwokqko90onbc9zosujtlikdd3r9iq
hgat <- read_tsv(file.path(data_dir,"ALT_May_2021_JS_plus_v19_histologies.tsv")) %>%
  mutate(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID_DNA) %>%
  filter(short_histology == "HGAT")

tmb <- read_tsv(file.path(data_dir,"pbta-snv-consensus-TMB_intarget.txt"))

hgat <- hgat %>%
  left_join(tmb,by=c("Tumor_Sample_Barcode"))

# select 1 id from multiple RNA 
hgat_uniq_rna <- hgat %>%
  group_by(Kids_First_Participant_ID,tumor_descriptor) %>%
  summarise(Kids_First_Biospecimen_ID_RNA=sample(Kids_First_Biospecimen_ID_RNA,1),
            Kids_First_Biospecimen_ID_DNA=sample(Kids_First_Biospecimen_ID_DNA,1))

hgat <- arrange(hgat,as.numeric(telomere_ratio)) %>%
  filter(Kids_First_Biospecimen_ID_RNA %in% hgat_uniq$Kids_First_Biospecimen_ID_RNA) %>%
  as.data.frame() %>% remove_rownames() 

###read in maf
keep_cols <- c("Hugo_Symbol", 
               "Chromosome", 
               "Start_Position", 
               "End_Position", 
               "Reference_Allele", 
               "Tumor_Seq_Allele2", 
               "Variant_Classification", 
               "Variant_Type",
               "Tumor_Sample_Barcode",
               "VAF")

consensus <- data.table::fread(file.path(data_dir, "pbta-snv-consensus-mutation.maf.tsv.gz"),
                      select = keep_cols, data.table = FALSE)

# Hotspots were copied over from OpenPBTA-analysis/analyses/hotspots-detection/results
hotspots <- data.table::fread(file.path(data_dir,"pbta-snv-scavenged-hotspots.maf.tsv.gz"),
                     select = keep_cols, data.table = FALSE)

full_maf <- consensus %>%
  bind_rows(hotspots) %>%
  # filter for HGAT
  filter(Tumor_Sample_Barcode %in% hgat$Tumor_Sample_Barcode) %>%
  unique()

##read in cnv copied over from OpenPBTA-analysis/analyses/focal-cn-file-preparation/results
cnv_autosomes_df <- readr::read_tsv(file.path(data_dir, "consensus_seg_annotated_cn_autosomes.tsv.gz")) %>%
  left_join(select(hgat,c("Tumor_Sample_Barcode","germline_sex_estimate")),
            by=c("biospecimen_id"="Tumor_Sample_Barcode")
  )
cnv_xy_df <- readr::read_tsv(file.path(data_dir, "consensus_seg_annotated_cn_x_and_y.tsv.gz")) 
cnv_df <- rbind(cnv_autosomes_df,cnv_xy_df)
cnv_df <- cnv_df %>%
  mutate(Hugo_Symbol = gene_symbol,
         Tumor_Sample_Barcode = biospecimen_id,
         Variant_Classification = status) %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification) %>%
  # mutate loss and amplification to Del and Amp to fit Maftools format
  dplyr::mutate(Variant_Classification = dplyr::case_when(Variant_Classification == "loss" ~ "Del",
                                                          Variant_Classification == "amplification" ~ "Amp",
                                                          TRUE ~ as.character(Variant_Classification))) %>%
  # only keep Del and Amp calls
  filter(Variant_Classification %in% c("Del", "Amp")) %>%
  # filter for HGAT and gene of interest
  filter(Tumor_Sample_Barcode %in% hgat$Tumor_Sample_Barcode) %>%
  distinct()

cnv_df <- cnv_df %>%
  filter(!is.na(Variant_Classification))%>%
  distinct() %>%
  as.data.frame()

# complex heatmap
collapse_snv_dat <- select(full_maf,c(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification)) %>%
  unique() %>%
  filter(Variant_Classification %in% names(colors)) %>%
  group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  summarise(count = as.double(length(Variant_Classification[!is.na(Variant_Classification)])),
            Variant_Classification=paste(Variant_Classification,collapse = ",")) 

collapse_cnv_dat <- cnv_df %>%
  group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  summarise(count = as.double(length(Variant_Classification[!is.na(Variant_Classification)])),
            Variant_Classification=paste(Variant_Classification,collapse = ",")) 


merged_dat <- collapse_snv_dat %>%
  full_join(collapse_cnv_dat, by=c("Hugo_Symbol","Tumor_Sample_Barcode"),suffix=c("_SNV","_CNV")) %>%
  ungroup() %>%
  full_join(select(hgat,c("Tumor_Sample_Barcode","sample_id"))) %>%
  replace_na(list(Variant_Classification_SNV="",Variant_Classification_CNV="")) %>%
  mutate(count=if_else(Variant_Classification_CNV != "",count_SNV+count_CNV, count_SNV)) %>%
  unique()


merged_dat$Variant_Classification <- apply( merged_dat[ , c("Variant_Classification_SNV","Variant_Classification_CNV") ] , 1 , paste , collapse = "," ) 
merged_dat <- merged_dat %>%
  mutate(Variant_Classification = case_when(
    count>1~"Multi_Hit",
    TRUE ~ Variant_Classification))

gene_matrix<-reshape2::acast(merged_dat,
                             Hugo_Symbol~Tumor_Sample_Barcode,value.var = "Variant_Classification")

saveRDS(gene_matrix,"hgat_snv_cnv_alt_matrix.RDS")
write_tsv(hgat,"hgat_subset.tsv")

    
    
    