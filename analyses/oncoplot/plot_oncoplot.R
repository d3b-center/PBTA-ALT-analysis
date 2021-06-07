# Setup dir
input_dir <- "input"
data_dir <- "../../data"

##read in input
source(file.path(input_dir, "mutation-colors.R"))
goi.list <- read_tsv(file.path(input_dir, "goi-mutations"), col_names = "genes")

# read processed files
hgat <- read_tsv("hgat_subset.tsv") %>% 
  arrange(telomere_ratio) %>% 
  column_to_rownames("Tumor_Sample_Barcode")
gene_matrix<-readRDS("hgat_snv_cnv_alt_matrix.RDS")

gene_matrix<- gene_matrix[goi.list$genes,hgat$Kids_First_Biospecimen_ID_DNA]


## color for barplot
col = colors

ha = HeatmapAnnotation( name = "annotation", df = hgat[,c("telomere_ratio","CCA Binary","tumor_descriptor")],
                        "TMB"=anno_barplot(hgat$TMB,ylim = c(0,10)),
                        col=list(
                          "telomere_ratio" = colorRamp2(c(0, 1,11), c("white", "blue","darkblue")),
                          "CCA Binary" = c("1"="#CD96CD",
                                           "0"="#000000",
                                           "NA" = "#f0efef")),
                      annotation_name_side = "right",annotation_name_gp = gpar(fontsize = 9))

hgat_bt_anno = hgat[,c("ATRX_fpkm","DAXX_fpkm","TERT_fpkm")] %>%
  mutate("zscore_ATRX_fpkm" = scale(ATRX_fpkm),
         "zscore_TERT_fpkm" = scale(TERT_fpkm),
         "zscore_DAXX_fpkm" = scale(DAXX_fpkm)) %>%
  select("zscore_ATRX_fpkm","zscore_DAXX_fpkm","zscore_TERT_fpkm")

ha1 = HeatmapAnnotation( df = hgat_bt_anno,
                        height = unit(3, "cm"),
                        col=list(
                          "zscore_ATRX_fpkm"= colorRamp2(c(-5, 0, 5), c("green", "white", "red")),
                          "zscore_DAXX_fpkm" = colorRamp2(c(-5, 0, 5), c("green", "white", "red")),
                          "zscore_TERT_fpkm" = colorRamp2(c(-5, 0, 5), c("green", "white", "red"))
                        ),
                        annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 9))

pdf("telomere_hgat.pdf", height = 9, width = 20)
oncoPrint(gene_matrix, get_type = function(x) strsplit(x, ",")[[1]],
          column_names_gp = gpar(fontsize = 9),show_column_names = F,show_row_barplot = F,
          alter_fun = list(
            background = function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "white",col="white")),
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
            Del = function(x, y, w, h) grid.rect(x, y, w*0.75, h*0.5, gp = gpar(fill = unname(col["Del"]), col = NA)),
            Amp = function(x, y, w, h) grid.rect(x, y, w*0.75, h*0.5, gp = gpar(fill = unname(col["Amp"]), col = NA))),
          col = col,
          top_annotation = ha,
          bottom_annotation = ha1,
          column_order = NULL
)

dev.off()

