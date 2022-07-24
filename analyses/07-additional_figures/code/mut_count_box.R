
library(readr)
library(tidyverse)
library(ggpubr)
library(openxlsx)


# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "07-additional_figures")
plots_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(root_dir, "analyses", "05-oncoplot", "input")
output_dir <- file.path(analysis_dir, "output")
maf_dir <- file.path(root_dir, "analyses", "09-lollipop", "output")

# source mutations used for oncoprint
source(file.path(input_dir, "mutation-colors.R"))
mut_of_interest <- c(names(colors), "3'UTR", "5'UTR", "Splice_Region")

# metadata read in
goi <- read_table(file.path(root_dir, "analyses", "05-oncoplot", "input", "goi-mutations"), col_names = F)
  
metadata <- read_tsv(file.path(root_dir, "analyses", "02-add-histologies", "output",
                               "stundon_hgat_updated_hist_alt.tsv")) %>%
  filter(short_histology == "HGAT") %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID_DNA) %>%
  dplyr::mutate(`alt final` = case_when(`alt final` == "pos" ~ "POS",
                            `alt final` == "neg" ~ "NEG",
                            TRUE ~ as.character(`alt final`)) 
  ) %>%
  dplyr::select(Tumor_Sample_Barcode, Kids_First_Biospecimen_ID_RNA, sample_id, `alt final`, atrx_mut, telomere_ratio)

# get TMB 
tmb_coding <- read_tsv(file.path(data_dir,"pbta-snv-consensus-mutation-tmb-coding.tsv")) %>%
  dplyr::select(Tumor_Sample_Barcode, tmb)

# read in annotated MAF
maf <- read_tsv(file.path(maf_dir, "snv-consensus-plus-hotspots-hgat-oncokb.maf.tsv"))

maf %>%
  filter(ONCOGENIC == "Unknown") %>%
  count(Variant_Classification)

maf %>%
  filter(ONCOGENIC != "Unknown") %>%
  count(Variant_Classification)

# read in mut sigs
sigs <- readxl::read_excel(file.path(data_dir, "TableS2-DNA-results-table.xlsx"), sheet = 2) %>%
  rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID)

# read in telomerase scores
tel <- readxl::read_excel(file.path(data_dir, "TableS3-RNA-results-table.xlsx"), sheet = 2) %>%
  rename(telomerase_score = NormEXTENDScores_fpkm) %>%
  select(Kids_First_Biospecimen_ID_RNA, telomerase_score)

names(tel)

maf_reanno <- maf %>%
  filter(Variant_Classification %in% mut_of_interest) %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, ONCOGENIC) %>%
  unique() %>%
  # pull together if multiple annotations per TSB
  group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  summarise(ONCOGENIC = str_c(unique(ONCOGENIC), collapse=";")) %>%
  # recode as likely oncogenic for those with both
  mutate(ONCOGENIC = case_when(grepl(";", ONCOGENIC) ~ "Oncogenic or Likely Oncogenic",
                               ONCOGENIC == "Unknown" ~ "VUS",
                               ONCOGENIC == "Likely Oncogenic" ~ "Oncogenic or Likely Oncogenic",
                               ONCOGENIC == "Oncogenic" ~ "Oncogenic or Likely Oncogenic",
                               TRUE ~ as.character(ONCOGENIC))) %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, ONCOGENIC) %>%
  unique()


alt_pos <- metadata %>% filter(`alt final` == "POS") %>% nrow()
alt_neg <- metadata %>% filter(`alt final` == "NEG") %>% nrow()

for (gene in goi$X1) {
  maf_reanno_goi <- maf_reanno %>%
  filter(Hugo_Symbol == gene)
  
  # add non-mutated
  if (gene == "ATRX"){
  metadata_atrx <- metadata %>%
    left_join(maf_reanno_goi) %>%
    mutate(ONCOGENIC = case_when(is.na(ONCOGENIC) ~ "Not mutated",
                                 TRUE ~ as.character(ONCOGENIC)),
           atrx_mut = case_when(ONCOGENIC == "Not mutated" ~ "nonATRX_mut",
             TRUE ~ "ATRX_mut"),
    )
  
  table_df_atrx <- as.data.frame(table(metadata_atrx$`alt final`, metadata_atrx$ONCOGENIC))
  colnames(table_df_atrx) <- c("ALT", "Mutation Type", "Count")
  
  # add freq
  table_df_atrx <- table_df_atrx %>%
    mutate(Frequency = case_when(ALT == "POS" ~ Count/alt_pos,
                                 ALT == "NEG" ~ Count/alt_neg))
  
  write_tsv(table_df_atrx, file.path(output_dir, paste0("alt_ATRX_counts_by_oncogenicity.tsv")))
  
  }
  if (gene != "ATRX"){
    metadata_maf <- metadata %>%
      left_join(maf_reanno_goi) %>%
      mutate(ONCOGENIC = case_when(is.na(ONCOGENIC) ~ "Not mutated",
                                   TRUE ~ as.character(ONCOGENIC))) 
    table_df <- as.data.frame(table(metadata_maf$`alt final`, metadata_maf$ONCOGENIC)) 
    colnames(table_df) <- c("ALT", "Mutation Type", "Count")
  
    # add freq
    table_df <- table_df %>%
      mutate(Frequency = case_when(ALT == "POS" ~ Count/alt_pos,
                                 ALT == "NEG" ~ Count/alt_neg))
    
    write_tsv(table_df, file.path(output_dir, paste0("alt_", gene, "_counts_by_oncogenicity.tsv")))
  }
}



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
  dplyr::rename(`TMB Coding` = tmb) 


################# Generate figures with combined mutation counts 
# generate count dataframe
count_df <- merged_dat %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::mutate(mut_count = sum(count)) %>%
  dplyr::select(Tumor_Sample_Barcode, mut_count, `TMB Coding`) %>%
  distinct() %>%
  dplyr::left_join(metadata_atrx) %>%
  dplyr::left_join(tel) %>%
  dplyr::mutate(log2_mut_count = log2(mut_count))

count_df$`alt final` <- factor(count_df$`alt final`, levels = c("POS", "NEG")) 

# add mutational sigs
count_df_sigs <- count_df %>%
  left_join(sigs) %>%
  rename(Kids_First_Biospecimen_ID_DNA = Tumor_Sample_Barcode) %>%
  select(Kids_First_Biospecimen_ID_DNA, Kids_First_Biospecimen_ID_RNA, sample_id, telomere_ratio, `alt final`,
         atrx_mut, ONCOGENIC, `TMB Coding`, log2_mut_count, everything(sigs))

write.xlsx(count_df_sigs, 
           file.path(output_dir, "hgat_tmb_mut_counts_sigs_tel_df.xlsx"), 
             overwrite=TRUE, 
             keepNA=TRUE)

# Remove hypermutant for plots
count_df_for_plots <- count_df %>%
  filter(`TMB Coding` < 10)
# output plots for all mutation coutns
pdf(file.path(plots_dir, "mut_count_alt_all_genes.pdf"))
p <- ggplot(count_df_for_plots, aes(x =`alt final`, y = log2_mut_count)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test') +
  theme_bw() + 
  ylab("Log2 Mutation Count")
  
print(p)
dev.off()

# output plots for TMB 
pdf(file.path(plots_dir, "tmb_alt_all_genes.pdf"))
p <- ggplot(count_df_for_plots, aes(x =`alt final`, y = `TMB Coding`)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test') +
  theme_bw()

print(p)
dev.off()

pdf(file.path(plots_dir, "tmb_atrx_all_genes.pdf"))
p <- ggplot(count_df_for_plots, aes(x = atrx_mut, y = `TMB Coding`)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test') +
  theme_bw()

print(p)
dev.off()


pdf(file.path(plots_dir, "atrx_alt_stacked_by_oncogenicity.pdf"))
p <- ggplot(table_df_atrx, aes(x = ALT, y = Count, fill = `Mutation Type`)) +
  geom_bar(stat = "identity") + 
  scale_fill_brewer(palette="Paired") +
  theme_bw() +
  ylab("Number of tumors") +
  xlab("ALT status")
print(p)
dev.off()

pdf(file.path(plots_dir, "atrx_alt_faceted_by_oncogenicity.pdf"))
p <- ggplot(table_df_atrx, aes(x = ALT, y = Count, fill = `Mutation Type`)) +
  geom_bar(stat = "identity") +
  facet_wrap(~`Mutation Type`)+
  scale_fill_brewer(palette="Paired")+
  theme_bw() +
  ylab("Number of tumors") +
  xlab("ALT status")
print(p)
dev.off()



# output plots for TMB 
pdf(file.path(plots_dir, "tmb_alt_all_genes_atrx.pdf"))
p <- ggplot(count_df_for_plots, aes(x =`alt final`, y = `TMB Coding`, color = atrx_mut)) +
  geom_boxplot() + 
  geom_jitter() + 
  stat_compare_means(method='t.test') +
  theme_bw()

print(p)
dev.off()


################# Generate figures mutation counts faceted by type
# generate count dataframe
count_df_facet <- merged_dat %>%
  filter(`TMB Coding` < 10) %>%
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

cor.test(count_df$telomere_ratio, count_df$telomerase_score)

pdf(file.path(plots_dir, "alt_tel_cor.pdf"))
p <- ggplot(count_df, aes(x = telomere_ratio, y = telomerase_score)) +
  geom_jitter() + 
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x)  +
  theme_bw()

print(p)
dev.off()

