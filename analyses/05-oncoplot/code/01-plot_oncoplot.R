# load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)

# define directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "05-oncoplot")
input_dir <- file.path(analysis_dir, "input")
output_dir <- file.path(analysis_dir, "output")
data_dir <- file.path(root_dir, "data")

##read in input
source(file.path(input_dir, "mutation-colors.R"))
goi.list <- read_tsv(file.path(input_dir, "goi-mutations"), col_names = "genes")
germline <- read_tsv(file.path(input_dir, "germline_variants_meta_format.tsv"))
me3 <- read_tsv(file.path(input_dir, "h3k27me3_tma.tsv")) %>%
  dplyr::rename(`H3K27me3 IHC` = H3K27me3)
ihc <- readxl::read_excel(file.path(input_dir, "TMA table for HGAT paper_052722_kac.xlsx")) %>%
  rename(sample_id = ID,
         `ATRX IHC` = `ATRX IHC (Pathology)`,
         `Telomeric foci` = `Presence of UBTF`) %>%
  mutate(`ATRX IHC` = case_when(`ATRX IHC` == 0 ~ "POS",
                                `ATRX IHC` == 1 ~ "NEG",
                                TRUE ~ "Not done"),
         `Telomeric foci` = case_when(`Telomeric foci` == 1 ~ "POS",
                                        `Telomeric foci` == 0 ~ "NEG",
                                TRUE ~ "Not done")) %>%
  full_join(me3) %>%
  mutate(`H3K27me3 IHC` = case_when(is.na(`H3K27me3 IHC`) ~ "Not done",
                                    TRUE ~ `H3K27me3 IHC`)) %>%
  # remove those without research ID
  filter(!is.na(`Research Subject ID`))


# read in telomerase scores
tel <- readxl::read_excel(file.path(data_dir, "TableS3-RNA-results-table.xlsx"), sheet = 2) %>%
  rename(`Telomerase score` = NormEXTENDScores_fpkm) %>%
  select(Kids_First_Biospecimen_ID_RNA, `Telomerase score`)


# read processed files
hgat <- read_tsv(file.path(input_dir,"hgat_subset.tsv")) %>% 
  select(-`ATRX IHC`) %>%
  left_join(germline, by = "sample_id") %>%
  mutate(`Germline MMR` = case_when(is.na(`Germline MMR`) ~ "no",
                                    TRUE ~ as.character(`Germline MMR`)),
         `Somatic MMR` = case_when(is.na(`Somatic MMR`) ~ "no",
                                   TRUE ~ as.character(`Somatic MMR`))
         ) %>%
  left_join(ihc, by = "sample_id") %>%
  left_join(tel, by = "Kids_First_Biospecimen_ID_RNA") %>%
  dplyr::rename(`C-circle` = `CCA Sept 2021`) %>%
  column_to_rownames("Tumor_Sample_Barcode") %>%
  arrange(telomere_ratio)

gene_matrix<- readRDS(file.path(input_dir,"hgat_snv_cnv_alt_matrix.RDS"))
gene_matrix <- gene_matrix[goi.list$genes,]

tmb <- read_tsv(file.path(data_dir,"pbta-snv-consensus-mutation-tmb-coding.tsv")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Tumor_Sample_Barcode) %>%
  select(Kids_First_Biospecimen_ID_DNA, tmb)

# mutate the hgat dataframe for plotting
hgat <- hgat %>%
  dplyr::left_join(tmb) %>%
  dplyr::mutate(`C-circle` = case_when(
    `C-circle` %in% c("POS", "NEG") ~ `C-circle`,
    TRUE ~ "Not done"
  )) %>%
  dplyr::mutate(
    TMB = case_when(
      tmb < 10 ~ "Normal",
      tmb >= 10 & tmb < 100 ~ "Hypermutant",
      tmb>= 100 ~ "Ultra-hypermutant")
  ) %>%
  dplyr::rename(`Phase of therapy` = tumor_descriptor,
                Sex = germline_sex_estimate,
                `Telomere ratio` = telomere_ratio )


# order columns for plotting
hgat$`C-circle` <- factor(hgat$`C-circle`, levels = c("POS", "NEG", "Not done"))
hgat$Sex <- factor(hgat$Sex, levels = c("Male", "Female"))
hgat$TMB <- factor(hgat$TMB, levels = c("Ultra-hypermutant", "Hypermutant", "Normal"))


## color for barplot
col = colors
df = hgat[,c("Kids_First_Biospecimen_ID_DNA", "Sex","Phase of therapy", "Telomere ratio", "C-circle", "ATRX IHC", "H3K27me3 IHC", "Telomeric foci", "TMB", "Germline MMR", "Somatic MMR")]

colorder <- df$Kids_First_Biospecimen_ID_DNA

#subset for what's in the meta file/order matrix
gene_matrix_df <- as.data.frame(gene_matrix)
gene_matrix_ordered <- gene_matrix_df %>%
  select(all_of(colorder)) %>%
  as.matrix()
  

# check if in same order                           
identical(colnames(gene_matrix_ordered), df$Kids_First_Biospecimen_ID_DNA)

palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


# match annotations and gene matrix by bs_id with 
ha = HeatmapAnnotation(name = "annotation", df = hgat[,c("Sex","Phase of therapy", "Telomere ratio", "C-circle", "Telomeric foci", "ATRX IHC", "H3K27me3 IHC", "TMB", "Germline MMR", "Somatic MMR")],
                       # "TMB"=anno_barplot(hgat$TMB, ylim = c(0,6), gp = gpar(fill = "#CCCCCC80")),
                        col=list(
                          "Sex" = c("Male" = "#56B4E9",
                                    "Female" = "#CC79A7"),
                          "Phase of therapy" = c("Initial CNS Tumor" = "#F0E442",
                                                 "Progressive" = "#56B4E9",
                                                 "Progressive Disease Post-Mortem" = "#009E73",
                                                 "Recurrence" = "#E69F00",
                                                 "Second Malignancy" = "#0072B2"),
                          #"Phase of therapy" = c("Initial CNS Tumor" = "#7FFFD4",
                          #                       "Progressive" = "#FFFFB5",
                          #                       "Progressive Disease Post-Mortem" = "#EEAEEE",
                          #                       "Recurrence" = "#ABDEE6",
                          #                       "Second Malignancy" = "#CBAACB"),
                          "Telomere ratio" = colorRamp2(c(0, 1.05, 1.06), c("whitesmoke", "#CAE1FF","#0072B2")),
                         # "Telomerase score" = colorRamp2(c(0, 0.5, 1.0), c("whitesmoke", "#E69F00", "#0072B2")),
                          "C-circle" = c("POS"="#0072B2",
                                         "NEG"="lightsteelblue1",
                                         "Not done" = "whitesmoke"),
                          "ATRX IHC" = c("POS"="#0072B2",
                                         "NEG"="lightsteelblue1",
                                         "Not done" = "whitesmoke"),
                         "H3K27me3 IHC" = c("POS"="#0072B2",
                                           "NEG"="lightsteelblue1",
                                           "Not done" = "whitesmoke"),
                          "Telomeric foci" = c("POS"="#0072B2",
                                         "NEG"="lightsteelblue1",
                                         "Not done" = "whitesmoke"),
                          "Germline MMR" = c("yes" = "#56B4E9",
                                             "no" = "whitesmoke"),
                          "Somatic MMR" = c("yes" = "#56B4E9",
                                            "no" = "whitesmoke"),
                          "TMB" = c("Ultra-hypermutant" = "#CC79A7", 
                                                "Hypermutant" = "#009E73", 
                                                "Normal" = "whitesmoke")),
                       annotation_name_side = "right", 
                       annotation_name_gp = gpar(fontsize = 9),
                       na_col = "whitesmoke")

#hgat_bt_anno = hgat[,c("ATRX_fpkm","DAXX_fpkm","TERT_fpkm")] %>%
 # mutate("zscore_ATRX_fpkm" = scale(ATRX_fpkm),
  #       "zscore_TERT_fpkm" = scale(TERT_fpkm),
   #      "zscore_DAXX_fpkm" = scale(DAXX_fpkm)) %>%
  #select("zscore_ATRX_fpkm","zscore_DAXX_fpkm","zscore_TERT_fpkm")

#ha1 = HeatmapAnnotation( df = hgat_bt_anno,
 #                       height = unit(3, "cm"),
  #                      col=list(
   #                       "zscore_ATRX_fpkm"= colorRamp2(c(-5, 0, 5), c("midnightblue", "white", "firebrick1")),
    #                      "zscore_DAXX_fpkm" = colorRamp2(c(-5, 0, 5), c("midnightblue", "white", "firebrick1")),
     #                     "zscore_TERT_fpkm" = colorRamp2(c(-5, 0, 5), c("midnightblue", "white", "firebrick1"))
      #                  ),
       #                 annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 9),
        #                na_col = "gainsboro")


#pdf(file.path(output_dir, "oncoprint_hgat.pdf"), height = 3, width = 15, onefile = FALSE)
tiff(file.path(output_dir, "oncoprint_hgat.tiff"), height = 900, width = 4500, units = "px", res = 300)
# global option to increase space between heatmap and annotations
ht_opt$ROW_ANNO_PADDING = unit(1, "cm")
oncoPrint(gene_matrix_ordered, get_type = function(x) strsplit(x, ",")[[1]],
          column_names_gp = gpar(fontsize = 9), show_column_names = F, #show_row_barplot = F,
          alter_fun = list(
            background = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "whitesmoke",col="whitesmoke")),
            Missense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Missense_Mutation"]),col = NA)),
            Nonsense_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Nonsense_Mutation"]),col = NA)),
            Frame_Shift_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Frame_Shift_Del"]), col = NA)),
            Frame_Shift_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Frame_Shift_Ins"]), col = NA)),
            Splice_Site = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Splice_Site"]), col = NA)),
            Translation_Start_Site = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Translation_Start_Site"]), col = NA)),
            Nonstop_Mutation = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Nonstop_Mutation"]),col = NA)),
            In_Frame_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["In_Frame_Del"]),col = NA)),
            In_Frame_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["In_Frame_Ins"]), col = NA)),
            Stop_Codon_Ins = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Stop_Codon_Ins"]), col = NA)),
            Start_Codon_Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Start_Codon_Del"]), col = NA)),
            Fusion = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Fusion"]),col = NA)),
            Multi_Hit_Fusion  = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Multi_Hit_Fusion"]),col = NA)),
            Multi_Hit = function(x, y, w, h) grid.rect(x, y, w*0.75, h*0.85, gp = gpar(fill = unname(col["Multi_Hit"]), col = NA)),
            Del = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Del"]), col = NA)),
            Amp = function(x, y, w, h) grid.rect(x, y, w*0.85, h*0.85, gp = gpar(fill = unname(col["Amp"]), col = NA))),
          col = col,
          top_annotation = ha,
          #bottom_annotation = ha1,
          column_order =  colnames(gene_matrix_ordered)
          )
dev.off()

