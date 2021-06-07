# Package names
cran_packages <- c("tidyverse","BiocManager","R.utils","rprojroot")

# Install packages not yet installed
installed_packages <- cran_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(cran_packages[!installed_packages],repos = "http://cran.us.r-project.org")
}

bioc_packages <- "maftools"
# Install Bioconductor packages not yet installed
installed_packages <- bioc_packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(bioc_packages[!installed_packages])
  
}

# Packages loading
invisible(lapply(c(cran_packages,bioc_packages), library, character.only = TRUE))

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

expression <- readRDS(file.path("../../data/pbta-gene-expression-rsem-fpkm-collapsed.stranded.rds"))


maf_obj <- maftools::read.maf(maf = full_maf, clinicalData = hgat, cnTable=cnv_df,
                   vc_nonSyn = c(
                     "Frame_Shift_Del",
                     "Frame_Shift_Ins",
                     "Splice_Site",
                     "Nonsense_Mutation",
                     "Nonstop_Mutation",
                     "In_Frame_Del",
                     "In_Frame_Ins",
                     "Missense_Mutation",
                     "Fusion",
                     "Multi_Hit",
                     "Multi_Hit_Fusion",
                     "Amp",
                     "Del",
                     "Intron",
                     "5'Flank",
                     "3'Flank"))

maf_obj@clinical.data[,telomere_ratio := as.numeric(as.character(telomere_ratio))]
str(maf_obj@clinical.data$telomere_ratio)
maf_obj@clinical.data[,ATRX_fpkm := as.numeric(as.character(ATRX_fpkm))]
str(maf_obj@clinical.data$ATRX_fpkm)
maf_obj@clinical.data[,TERT_fpkm := as.numeric(as.character(TERT_fpkm))]
str(maf_obj@clinical.data$TERT_fpkm)
maf_obj@clinical.data[,DAXX_fpkm := as.numeric(as.character(DAXX_fpkm))]
str(maf_obj@clinical.data$DAXX_fpkm)

pdf("telomere_hgat.pdf", height = 15, width = 15)
    maftools::oncoplot(maf = maf_obj, genes = goi.list$genes, 
                           removeNonMutated = F, bgCol = "whitesmoke", 
                           showTumorSampleBarcodes = F, drawRowBar = T,
                           annotationFontSize = 10, 
                           sortByAnnotation = T, fontSize = 10, legendFontSize =20,
                           colors = colors, writeMatrix = F,
                           clinicalFeatures = c("CCA_Binary","telomere_ratio","ATRX_fpkm","TERT_fpkm","DAXX_fpkm"))
    dev.off()
    
    
    
    