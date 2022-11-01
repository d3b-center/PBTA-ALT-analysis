# Package names
cran_packages <- c("tidyverse","R.utils","rprojroot","devtools")
# Packages loading
invisible(lapply(c(cran_packages,"ComplexHeatmap"), library, character.only = TRUE))

# define directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "05-oncoplot")
input_dir <- file.path(analysis_dir, "input")
output_dir <- file.path(analysis_dir, "output")

source(file.path(input_dir, "mutation-colors.R"))

##read in input
goi_list <- read_tsv(file.path(input_dir, "goi-mutations"), col_names = "genes")
kegg_list <- read_tsv(file.path(input_dir, "KEGG_MISMATCH_REPAIR.txt"), col_names = "genes")

hgat <- read_tsv(file.path(root_dir,
                           "analyses",
                           "02-add-histologies",
                           "output",
                           "stundon_hgat_updated_hist_alt.tsv")) %>%
  mutate(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID_DNA) %>%
  filter(short_histology == "HGAT")

# format for MAF

###read in maf
keep_cols <- c("Hugo_Symbol", 
               "Chromosome", 
               "Start_Position", 
               "End_Position", 
               "Reference_Allele", 
               "Tumor_Seq_Allele2", 
               "Variant_Classification", 
               "Variant_Type",
               "Tumor_Sample_Barcode")

consensus <- data.table::fread(file.path(data_dir, "snv-consensus-plus-hotspots.maf.tsv.gz"),
                      select = keep_cols, data.table = FALSE) %>% 
  filter(Tumor_Sample_Barcode %in% hgat$Tumor_Sample_Barcode) %>% # filter for HGAT
  unique()


# are there somatic MMR alterations?
consensus_mmr <- consensus %>%
  filter(Variant_Classification %in% names(colors),
         Hugo_Symbol %in% kegg_list$genes) %>%
  left_join(hgat[,c("Tumor_Sample_Barcode", "sample_id")])


## read in cnv (OpenPedCan v11)
cnv_df <- readr::read_tsv(file.path(data_dir, "consensus_wgs_plus_cnvkit_wxs.tsv.gz")) %>%
  # filter for HGAT and gene of interest
  filter(biospecimen_id %in% hgat$Tumor_Sample_Barcode) %>%
  mutate(Hugo_Symbol = gene_symbol,
         Tumor_Sample_Barcode = biospecimen_id,
         Variant_Classification = status) %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode, Variant_Classification) %>%
  # mutate loss and amplification to Del and Amp to fit Maftools format
  dplyr::mutate(Variant_Classification = dplyr::case_when(Variant_Classification == "deep deletion" ~ "Del",
                                                          Variant_Classification == "amplification" ~ "Amp",
                                                          TRUE ~ as.character(Variant_Classification))) %>%
  # only keep Del and Amp calls
  filter(Variant_Classification %in% c("Del", "Amp")) %>%
  distinct()

cnv_df <- cnv_df %>%
  filter(!is.na(Variant_Classification))%>%
  distinct() %>%
  as.data.frame()

# complex heatmap
collapse_snv_dat <- select(consensus,c(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification)) %>%
  unique() %>%
  filter(Variant_Classification %in% names(colors)) %>%
  group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  summarise(count = as.double(length(Variant_Classification[!is.na(Variant_Classification)])),
            Variant_Classification=paste(Variant_Classification,collapse = ",")) 

collapse_cnv_dat <- cnv_df %>%
  group_by(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  summarise(count = as.double(length(Variant_Classification[!is.na(Variant_Classification)])),
            Variant_Classification=paste(Variant_Classification,collapse = ",")) 

# No fusion calls to add for this cohort
fusion_df <- read_tsv(file.path(data_dir, "pbta-fusion-putative-oncogenic.tsv")) %>%
  filter(Sample %in% hgat$Tumor_Sample_Barcode)


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
merged_dat %>% 
  saveRDS(file.path(input_dir, "merged_mut_data.RDS"))

gene_matrix<-reshape2::acast(merged_dat,
                             Hugo_Symbol~Tumor_Sample_Barcode,value.var = "Variant_Classification")

saveRDS(gene_matrix, file.path(input_dir, "hgat_snv_cnv_alt_matrix.RDS"))
write_tsv(hgat, file.path(input_dir, "hgat_subset.tsv"))

    
    
    