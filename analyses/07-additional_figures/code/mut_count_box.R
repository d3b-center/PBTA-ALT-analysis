library(readr)
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(stats)


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
# source publication theme
source(file.path(root_dir, "analyses", "04-cutpoint-analysis", "code", "theme.R"))
mut_of_interest <- c(names(colors), "3'UTR", "5'UTR", "Splice_Region")

# metadata read in
goi <- read_table(file.path(root_dir, "analyses", "05-oncoplot", "input", "goi-mutations"), col_names = "genes")
  
metadata <- read_tsv(file.path(root_dir, "analyses", "02-add-histologies", "output",
                               "stundon_hgat_updated_hist_alt.tsv")) %>%
  filter(short_histology == "HGAT") %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID_DNA) %>%
  dplyr::mutate(`alt final` = case_when(`alt final` == "pos" ~ "POS",
                            `alt final` == "neg" ~ "NEG",
                            TRUE ~ as.character(`alt final`)) 
  ) %>%
  dplyr::select(Tumor_Sample_Barcode, Kids_First_Biospecimen_ID_RNA, sample_id, `alt final`, telomere_ratio)

# get TMB 
tmb_coding <- read_tsv(file.path(data_dir,"pbta-snv-consensus-mutation-tmb-coding.tsv")) %>%
  dplyr::select(Tumor_Sample_Barcode, tmb)

# read in annotated MAF
maf <- read_tsv(file.path(maf_dir, "snv-consensus-plus-hotspots-hgat-oncokb.maf.tsv")) %>%
  filter(Variant_Classification %in% mut_of_interest) %>%
  filter(Hugo_Symbol %in% goi$genes) %>%
  mutate(gene_var = paste(Hugo_Symbol, Variant_Classification, sep = "_"))
table(maf$gene_var, maf$ONCOGENIC)

tert <- maf %>%
  filter(Hugo_Symbol == "TERT")

maf %>%
  filter(ONCOGENIC == "Unknown") %>%
  count(Variant_Classification)

maf %>%
  filter(ONCOGENIC != "Unknown") %>%
  count(Variant_Classification)

maf %>%
  count(ONCOGENIC)

# read in mut sigs
sigs <- readxl::read_excel(file.path(data_dir, "TableS2-DNA-results-table.xlsx"), sheet = 2) %>%
  dplyr::rename(Tumor_Sample_Barcode = Kids_First_Biospecimen_ID)

# read in telomerase scores
tel <- readxl::read_excel(file.path(data_dir, "TableS3-RNA-results-table.xlsx"), sheet = 2) %>%
  rename(telomerase_score = NormEXTENDScores_fpkm) %>%
  select(Kids_First_Biospecimen_ID_RNA, telomerase_score)

maf_reanno <- maf %>%
  # remove 5'Flank from genes other than TERT, as they were not used initially
  filter(Variant_Classification != "5'Flank") %>%
  bind_rows(tert) %>%
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

for (gene in goi$genes) {
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
  rename(Kids_First_Biospecimen_ID_DNA = Tumor_Sample_Barcode)
count_df_sigs <- count_df_sigs %>%
  select(Kids_First_Biospecimen_ID_DNA, Kids_First_Biospecimen_ID_RNA, sample_id, telomere_ratio, `alt final`,
         atrx_mut, ONCOGENIC, `TMB Coding`, log2_mut_count, everything(count_df_sigs), -mut_count)

write.xlsx(count_df_sigs, 
           file.path(output_dir, "hgat_tmb_mut_counts_sigs_tel_df.xlsx"), 
             overwrite=TRUE, 
             keepNA=TRUE)

# Remove hypermutant for plots
count_df_for_plots <- count_df_sigs %>%
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
p <- ggplot(table_df_atrx, aes(x = ALT, y = Frequency, fill = `Mutation Type`)) +
  geom_bar(stat = "identity") + 
  scale_fill_brewer(palette="Paired") +
  theme_bw() +
  ylab("Proportion of tumors") +
  xlab("ALT status")
print(p)
dev.off()

pdf(file.path(plots_dir, "atrx_alt_faceted_by_oncogenicity.pdf"))
p <- ggplot(table_df_atrx, aes(x = ALT, y = Frequency, fill = `Mutation Type`)) +
  geom_bar(stat = "identity") +
  facet_wrap(~`Mutation Type`)+
  scale_fill_brewer(palette="Paired")+
  theme_bw() +
  ylab("Proportion of tumors") +
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


# Do telhunt scores inversely correlate with telomerase scores? No
cor.test(count_df$telomere_ratio, count_df$telomerase_score)

pdf(file.path(plots_dir, "alt_tel_cor.pdf"))
p <- ggplot(count_df, aes(x = telomere_ratio, y = telomerase_score)) +
  geom_jitter() + 
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x)  +
  theme_bw()

print(p)
dev.off()

########### Assess ATRX clonality ############

# Create density plots for ATRX VAF by ALT status
metadata_atrx_vaf <- metadata_atrx %>%
  dplyr::rename(ALT_status = `alt final`) %>%
  left_join(tmb_coding) %>%
  left_join(tel) %>%
  # remove hypermutant samples and those without mutations
  filter(tmb < 10,
         !is.na(VAF)) %>%
  mutate(variable = case_when(ALT_status == "POS" & ONCOGENIC == "Oncogenic or Likely Oncogenic" ~ "ALT+ Oncogenic",
                              ALT_status == "POS" & ONCOGENIC == "VUS" ~ "ALT+ VUS",
                              ALT_status == "NEG" & ONCOGENIC == "Oncogenic or Likely Oncogenic" ~ "ALT- Oncogenic",
                              ALT_status == "NEG" & ONCOGENIC == "VUS" ~ "ALT- VUS"))
table(metadata_atrx_vaf$variable)

# get group means
mu <- plyr::ddply(metadata_atrx_vaf, "ALT_status", summarise, grp.mean=mean(VAF))
mu

# get group means
mu2 <- plyr::ddply(metadata_atrx_vaf, "variable", summarise, grp.mean=mean(VAF))
mu2

# plot
p1 <- ggplot(data=metadata_atrx_vaf, aes(x=VAF, group=ALT_status, color = ALT_status, after_stat(count))) +
  geom_density(adjust=1.5, alpha=.4) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=ALT_status),
             linetype="dashed") +
  xlim(0,1) +
  xlab("Somatic ATRX mutation VAF") +
  ylab("Count") +
  theme_Publication()

p2 <- ggplot(data=metadata_atrx_vaf, aes(x=VAF, group=variable, color = variable, after_stat(count))) +
  geom_density(adjust=1.5, alpha=.4) +
  geom_vline(data=mu2, aes(xintercept=grp.mean, color=variable),
             linetype="dashed") +
  xlim(0,1) +
  xlab("Somatic ATRX mutation VAF") +
  ylab("Count") +
  theme_Publication()


# correlate ATRXm VAF with telomere ratio
result <- cor.test(metadata_atrx_vaf$VAF, metadata_atrx_vaf$telomere_ratio, method = "pearson")
result

# Fit linear regression and create plot label
corr <- result$estimate
pval <-  result$p.value
label <- paste0("Adjusted Pearson's \n R = ", signif(corr, 3),
                ", p = ", signif(pval, 3))
label_grob <- grobTree(textGrob(label, x=0.05, y=0.9, hjust=0))

# plot correlation
p3 <- ggplot(metadata_atrx_vaf, aes(y = telomere_ratio, x = VAF)) +
  stat_smooth(fill="lightgray", method = "lm", col = "black", show.legend = FALSE) + 
  geom_point(size=4, color="black", stroke = 0.25) + 
  annotation_custom(label_grob) +
  ylim(c(0,5)) +
  xlim(c(0,1.2)) +
  xlab("Somatic ATRX mutation VAF") +
  ylab("T/N telomere content ratio") +
  theme_Publication()


tiff(file.path(plots_dir, "atrx_clonality.tiff"), height = 1400, width = 5000, res = 300)
p <- ggarrange(p1, p2, p3, widths = c(450, 500, 400), heights = c(150, 150, 150), align = "h", labels = c("B", "C", "D"), ncol = 3, nrow = 1)
print(p)
dev.off()

# Are the distributions significantly different? No

pos <- metadata_atrx_vaf %>%
  filter(ALT_status == "POS")
neg <- metadata_atrx_vaf %>%
  filter(ALT_status == "NEG")

ks.test(pos$VAF, neg$VAF)
#	Exact two-sample Kolmogorov-Smirnov test

# data:  pos$VAF and neg$VAF
# D = 0.35, p-value = 0.7283
# alternative hypothesis: two-sided

