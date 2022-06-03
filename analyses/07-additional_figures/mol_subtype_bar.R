library(readr)
library(tidyverse)
library(ggpubr)

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "07-additional_figures")
plots_dir <- file.path(analysis_dir, "plots")
results_dir <- file.path(analysis_dir, "output")
if(!dir.exists(results_dir)){
  dir.create(results_dir)
}
data_dir <- file.path(root_dir, "analyses", "02-add-histologies")

# metadata read in
metadata <- read_tsv(file.path(data_dir,
                               "output",
                               "stundon_hgat_03312022_updated_hist_alt.tsv")) %>%
  dplyr::select(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_DNA, `alt final`)

histology_pbta <- read_tsv(file.path(data_dir,
                                     "input-v21",
                                     "pbta-histologies.tsv")) %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype)

# calculate percentage of mutations for different subtypes 
num_pos <- metadata %>%
  dplyr::filter(`alt final` == "POS") %>% 
  nrow()
num_neg <- metadata %>%
  dplyr::filter(`alt final` == "NEG") %>% 
  nrow()

# calculate the percentage
combined_df <- metadata %>%
  left_join(histology_pbta) %>%
  dplyr::group_by(molecular_subtype) %>%
  dplyr::mutate(freq = case_when(
    `alt final` == "POS" ~ length(molecular_subtype) / num_pos,
    `alt final` == "NEG" ~ length(molecular_subtype) / num_neg)) %>%
  dplyr::mutate(mol_sub_group = case_when(
    molecular_subtype %in% c("DMG, H3 K28", 
                             "DMG, H3 K28, TP53 loss", 
                             "DMG, H3 K28, TP53 activated") ~ "H3 K28-mutant",
    molecular_subtype %in% c("HGG, H3 G35",
                             "HGG, H3 G35, TP53 loss") ~ "H3 G35-mutant",
    molecular_subtype %in% c("HGG, H3 wildtype",
                             "HGG, H3 wildtype, TP53 loss",
                             "HGG, H3 wildtype, TP53 activated") ~ "H3 Wild-type",
    molecular_subtype %in% c("HGG, IDH, TP53 loss",
                             "HGG, IDH, TP53 activated") ~ "IDH-mutant",
    TRUE ~ NA_character_
  ))
# reorder 
combined_df$`alt final` <- factor(combined_df$`alt final`, levels = c("POS", "NEG"))

# output the plot
pdf(file.path(plots_dir, "mol_subtype_stacked.pdf"))
p <- ggplot(combined_df, 
            aes(x=`alt final`, y=freq, fill=molecular_subtype)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()
print(p)
dev.off()

# generate chi square results
sink(file = file.path(results_dir, "chi-sqaure-sub-results.tsv"))
table(combined_df$molecular_subtype, combined_df$`alt final`)
chisq.test(table(combined_df$molecular_subtype, combined_df$`alt final`))
sink()

# output the plot
pdf(file.path(plots_dir, "mol_subtype_grouped_stacked.pdf"))
p <- ggplot(combined_df, 
            aes(x=`alt final`, y=freq, fill=mol_sub_group)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()
print(p)
dev.off()

# generate chi square results
sink(file = file.path(results_dir, "chi-sqaure-group-results.tsv"))
table(combined_df$mol_sub_group, combined_df$`alt final`)
chisq.test(table(combined_df$mol_sub_group, combined_df$`alt final`))
sink()


