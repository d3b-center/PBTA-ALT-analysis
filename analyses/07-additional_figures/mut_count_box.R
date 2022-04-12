
library(readr)
library(tidyverse)
library(ggpubr)

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "07-additional_figures")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(root_dir, "analyses", "05-oncoplot", "input")

# metadata read in
metadata <- read_tsv(file.path(root_dir, "analyses", "02-add-histologies", "output",
                               "stundon_hgat_03312022_updated_hist_alt.tsv")) %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID_DNA) %>%
  dplyr::select(Tumor_Sample_Barcode, `alt final`)

# get TMB 
tmb <- read_tsv(file.path(input_dir,"pbta-snv-mutation-tmb-coding.txt")) %>%
  select(Tumor_Sample_Barcode, tmb)

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
  dplyr::left_join(tmb) %>%
  dplyr::filter(tmb < 10)

################# Generate figures with combined mutation counts 
# generate count dataframe
count_df <- merged_dat %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::mutate(mut_count = sum(count)) %>%
  dplyr::select(Tumor_Sample_Barcode, mut_count) %>%
  distinct() %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(log2_mut_count = log2(mut_count))
count_df$`alt final` <- factor(count_df$`alt final`, levels = c("POS", "NEG"))

# output plots
pdf(file.path(plots_dir, "mut_count_alt_all_genes.pdf"))
p <- ggplot(count_df, aes(x =`alt final`, y = log2_mut_count)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test') 
  
print(p)
dev.off()

# filter to genes of interest
merged_dat_4_genes <- merged_dat %>%
  dplyr::filter(Hugo_Symbol %in% c("TP53",
                                   "H3F3A",
                                   "ATRX",
                                   "NF1"))

# generate count dataframe
count_df_4_genes <- merged_dat_4_genes %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::mutate(mut_count = sum(count)) %>%
  dplyr::select(Tumor_Sample_Barcode, mut_count) %>%
  distinct() %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(log2_mut_count = log2(mut_count))
count_df_4_genes$`alt final` <- factor(count_df_4_genes$`alt final`, levels = c("POS", "NEG"))

# output plots
pdf(file.path(plots_dir, "mut_count_alt_4_genes.pdf"))
p <- ggplot(count_df_4_genes, aes(x =`alt final`, y = log2_mut_count)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(method='t.test')

print(p)
dev.off()

################# Generate figures mutation counts faceted by type
# generate count dataframe
count_df_facet <- merged_dat %>%
  dplyr::group_by(Tumor_Sample_Barcode,
                  Variant_Classification) %>%
  dplyr::mutate(mut_count = sum(count)) %>%
  dplyr::select(Tumor_Sample_Barcode, mut_count, Variant_Classification) %>%
  distinct() %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(log2_mut_count = log2(mut_count))
count_df_facet$`alt final` <- factor(count_df_facet$`alt final`, levels = c("POS", "NEG"))

# output plots
pdf(file.path(plots_dir, "mut_count_alt_all_genes_faceted.pdf"))
p <- ggplot(count_df_facet, aes(x =`alt final`, y = log2_mut_count)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test', label.y = 7) + 
  facet_wrap( ~ Variant_Classification)

print(p)
dev.off()

# generate count dataframe
count_df_4_genes_faceted <- merged_dat_4_genes %>%
  dplyr::group_by(Tumor_Sample_Barcode,
                  Variant_Classification) %>%
  dplyr::mutate(mut_count = sum(count)) %>%
  dplyr::select(Tumor_Sample_Barcode, mut_count, Variant_Classification) %>%
  distinct() %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(log2_mut_count = log2(mut_count))
count_df_4_genes_faceted$`alt final` <- factor(count_df_4_genes_faceted$`alt final`, levels = c("POS", "NEG"))

# output plots
pdf(file.path(plots_dir, "mut_count_alt_4_genes_faceted.pdf"))
p <- ggplot(count_df_4_genes_faceted, aes(x =`alt final`, y = log2_mut_count)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(method='t.test') +
  facet_wrap( ~ Variant_Classification)

print(p)
dev.off()

