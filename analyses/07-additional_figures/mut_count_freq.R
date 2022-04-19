library(readr)
library(tidyverse)
library(ggpubr)

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "07-additional_figures")
plots_dir <- file.path(analysis_dir, "plots")
results_dir <- file.path(analysis_dir, "results", "count_freq_4_genes_stats")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}

input_dir <- file.path(root_dir, "analyses", "05-oncoplot", "input")

# metadata read in
metadata <- read_tsv(file.path(root_dir, "analyses", "02-add-histologies", "output",
                               "stundon_hgat_03312022_updated_hist_alt.tsv")) %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID_DNA) %>%
  dplyr::select(Tumor_Sample_Barcode, `alt final`)

# merged mutation matrix read in
merged_dat <- readRDS(file.path(input_dir, "merged_mut_data.RDS")) %>%
  dplyr::mutate(Variant_Classification = gsub(",", "", Variant_Classification)) %>%
  dplyr::filter(Variant_Classification %in% c("Missense_Mutation",
                                              "Nonsense_Mutation",
                                              "Frame_Shift_Del",
                                              "Frame_Shift_Ins",
                                              "Splice_Site",
                                              "Translation_Start_Site",
                                              "Nonstop_Mutation",
                                              "In_Frame_Del",
                                              "In_Frame_Ins",
                                              "Stop_Codon_Ins",
                                              "Start_Codon_Del",
                                              "Fusion",
                                              "Multi_Hit_Fusion",
                                              "Multi_Hit")) %>%
  dplyr::left_join(metadata) 

# count how many samples have a particular mutation
count_df <- merged_dat %>%
  dplyr::group_by(Hugo_Symbol,
                  `alt final`) %>%
  dplyr::mutate(sample_count = length(Tumor_Sample_Barcode)) %>%
  dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode, sample_count, `alt final`) %>%
  ungroup()

# get number of pos and neg
num_pos <- count_df %>%
  dplyr::filter(`alt final` == "POS") %>% 
  pull(Tumor_Sample_Barcode) %>% unique() %>% length()
num_neg <- count_df %>%
  dplyr::filter(`alt final` == "NEG") %>% 
  pull(Tumor_Sample_Barcode) %>% unique() %>% length()

# calculate frequency
count_df <- count_df %>%
  dplyr::mutate(freq = case_when(
    `alt final` == "POS" ~ sample_count / num_pos,
    `alt final` == "NEG" ~ sample_count / num_neg)) %>%
  distinct()

# filter to only 4 genes
count_df_4_genes <- count_df %>%
  dplyr::filter(Hugo_Symbol %in% c("TP53", 
                                   "ATRX", 
                                   "H3F3A", 
                                   "NF1")) %>%
  dplyr::select(Hugo_Symbol, sample_count, `alt final`, freq) %>%
  distinct()
count_df_4_genes$`alt final` <- factor(count_df_4_genes$`alt final`, levels = c("POS", "NEG"))
count_df_4_genes$Hugo_Symbol <- factor(count_df_4_genes$Hugo_Symbol, levels = c("TP53", "H3F3A", "ATRX", "NF1"))

# output plots
pdf(file.path(plots_dir, "freq_mut_barplot.pdf"))
p <- ggplot(data = count_df_4_genes, aes(x = Hugo_Symbol, y = freq, fill = `alt final`)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.75)
print(p)
dev.off()

# output the results 
for(gene in c("TP53", "ATRX", "H3F3A", "NF1")){
  count_each_gene <- count_df_4_genes %>%
    dplyr::filter(Hugo_Symbol == gene)
  
  sink(file = file.path(results_dir, paste0("wilcoxon_rank_sum_freq_", gene, ".txt")))
  print(pairwise.wilcox.test(count_each_gene$freq, 
                             count_each_gene$`alt final`, 
                             p.adjust.method = "bonf",
                             paired = TRUE))
  sink()
  
}
# output all combined
sink(file = file.path(results_dir, "wilcoxon_rank_sum_freq_all.txt"))
print(pairwise.wilcox.test(count_df_4_genes$freq, 
                           count_df_4_genes$`alt final`, 
                           p.adjust.method = "bonf",
                           paired = TRUE))
sink()

