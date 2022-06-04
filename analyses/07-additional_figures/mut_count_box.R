
library(readr)
library(tidyverse)
library(ggpubr)
library(openxlsx)


# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "07-additional_figures")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(root_dir, "analyses", "05-oncoplot", "input")
output_dir <- file.path(analysis_dir, "output")
# metadata read in
metadata <- read_tsv(file.path(root_dir, "analyses", "02-add-histologies", "output",
                               "stundon_hgat_updated_hist_alt.tsv")) %>%
  filter(short_histology == "HGAT") %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID_DNA) %>%
  dplyr::mutate(atrx_mut = case_when(
    !is.na(`ATRX Mutation`) ~ "ATRX_mut",
    TRUE ~ "non_ATRX_mut"),
    `alt final` = case_when(`alt final` == "pos" ~ "POS",
                            `alt final` == "neg" ~ "NEG",
                            TRUE ~ as.character(`alt final`)) 
  ) %>%
  dplyr::select(Tumor_Sample_Barcode, `alt final`, atrx_mut)

# get TMB 
tmb_coding <- read_tsv(file.path(input_dir,"pbta-snv-consensus-mutation-tmb-coding.tsv")) %>%
  dplyr::select(Tumor_Sample_Barcode, tmb)


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
  dplyr::left_join(tmb_coding) %>%
  dplyr::filter(tmb < 10) %>%
  dplyr::rename(`TMB Coding` = tmb) 

################# Generate figures with combined mutation counts 
# generate count dataframe
count_df <- merged_dat %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::mutate(mut_count = sum(count)) %>%
  dplyr::select(Tumor_Sample_Barcode, mut_count, `TMB Coding`) %>%
  distinct() %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(log2_mut_count = log2(mut_count))
count_df$`alt final` <- factor(count_df$`alt final`, levels = c("POS", "NEG")) 

write.xlsx(count_df, 
           file.path(output_dir, "hgat_tmb_mut_counts_df.xlsx"), 
             overwrite=TRUE, 
             keepNA=TRUE)


# output plots for all mutation coutns
pdf(file.path(plots_dir, "mut_count_alt_all_genes.pdf"))
p <- ggplot(count_df, aes(x =`alt final`, y = log2_mut_count)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test') +
  theme_bw() + 
  ylab("Log2 Mutation Count")
  
print(p)
dev.off()

# output plots for TMB 
pdf(file.path(plots_dir, "tmb_alt_all_genes.pdf"))
p <- ggplot(count_df, aes(x =`alt final`, y = `TMB Coding`)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test') +
  theme_bw()

print(p)
dev.off()

pdf(file.path(plots_dir, "tmb_atrx_all_genes.pdf"))
p <- ggplot(count_df, aes(x = atrx_mut, y = `TMB Coding`)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test') +
  theme_bw()

print(p)
dev.off()

# output plots for TMB 
pdf(file.path(plots_dir, "tmb_alt_all_genes_atrx.pdf"))
p <- ggplot(count_df, aes(x =`alt final`, y = `TMB Coding`, color = atrx_mut)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test') +
  theme_bw()

print(p)
dev.off()


################# Generate figures mutation counts faceted by type
# generate count dataframe
count_df_facet <- merged_dat %>%
  dplyr::group_by(Tumor_Sample_Barcode,
                  Variant_Classification) %>%
  dplyr::mutate(mut_count = sum(count)) %>%
  dplyr::select(Tumor_Sample_Barcode, mut_count, Variant_Classification, `TMB Coding`) %>%
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
  facet_wrap( ~ Variant_Classification) +
  theme_bw() + 
  ylab("Log2 Mutation Count")

print(p)
dev.off()

pdf(file.path(plots_dir, "tmb_alt_all_genes_faceted.pdf"))
p <- ggplot(count_df_facet, aes(x =`alt final`, y = `TMB Coding`)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test', label.y = 7) + 
  facet_wrap( ~ Variant_Classification) +
  theme_bw()

print(p)
dev.off()



