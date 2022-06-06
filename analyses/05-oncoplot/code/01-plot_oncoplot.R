# load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# define directories 
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "05-oncoplot")
input_dir <- file.path(analysis_dir, "input")
output_dir <- file.path(analysis_dir, "output")

##read in input
source(file.path(input_dir, "mutation-colors.R"))
goi.list <- read_tsv(file.path(input_dir, "goi-mutations"), col_names = "genes")

# read processed files
hgat <- read_tsv(file.path(input_dir,"hgat_subset.tsv")) %>% 
  arrange(telomere_ratio) %>% 
  column_to_rownames("Tumor_Sample_Barcode") %>%
  mutate(`C-circle` = `alt final`)
gene_matrix<- readRDS(file.path(input_dir,"hgat_snv_cnv_alt_matrix.RDS"))
tmb <- read_tsv(file.path(input_dir,"pbta-snv-consensus-mutation-tmb-coding.tsv")) %>%
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

#subset for what's in the meta file
gene_matrix<- gene_matrix[goi.list$genes, colnames(gene_matrix) %in% hgat$Kids_First_Biospecimen_ID_DNA]
setdiff(hgat$Kids_First_Biospecimen_ID_DNA, colnames(gene_matrix))

## color for barplot
col = colors
names(hgat)
df = hgat[,c("Sex","Phase of therapy", "Telomere ratio","C-circle", "TMB")]

ha = HeatmapAnnotation( name = "annotation", df = hgat[,c("Sex","Phase of therapy", "Telomere ratio","C-circle", "TMB")],
                       # "TMB"=anno_barplot(hgat$TMB, ylim = c(0,6), gp = gpar(fill = "#CCCCCC80")),
                        col=list(
                          "Sex" = c("Male" = "#CAE1FF",
                                                      "Female" = "#FFC1C1"),
                          "Phase of therapy" = c("Initial CNS Tumor" = "#7FFFD4",
                                                 "Progressive" = "#FFFFB5",
                                                 "Progressive Disease Post-Mortem" = "#EEAEEE",
                                                 "Recurrence" = "#ABDEE6",
                                                 "Second Malignancy" = "#CBAACB"),
                          "Telomere ratio" = colorRamp2(c(0, 1.05, 1.06), c("whitesmoke", "#CAE1FF","dodgerblue4")),
                          "C-circle" = c("POS"="dodgerblue4",
                                         "NEG"="whitesmoke",
                                         "Not done" = "gainsboro"),
                          "TMB" = c("Ultra-hypermutant" = "dodgerblue4", 
                                                "Hypermutant" = "darkorange1", 
                                                "Normal" = "whitesmoke")),
                      annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 9),
                      na_col = "gainsboro")

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

pdf(file.path(output_dir, "oncoprint_hgat.pdf"), height = 2, width = 15, onefile = FALSE)
oncoPrint(gene_matrix, get_type = function(x) strsplit(x, ",")[[1]],
          column_names_gp = gpar(fontsize = 9), show_column_names = F,#show_row_barplot = F,
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
          column_order =  colnames(gene_matrix)
          )

dev.off()

