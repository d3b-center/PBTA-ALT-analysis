library(tidyverse)
library(ggplot2)
library(ggforce)
library(forcats)

source(file.path(root_dir, "analyses", "04-cutpoint-analysis", "code", "theme.R"))
# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "10-distribution")
input_dir <- file.path(root_dir, "analyses", "01-alterations_ratio_check", "input")
output_dir <- file.path(analysis_dir, "output")
hist_dir <- file.path(analysis_dir, "v22-input")


# telhunter results
telhunt <- read_tsv(file.path(input_dir, "telomere_940_ratio.tsv")) %>%
  rename(Kids_First_Biospecimen_ID = Kids_First_Biospecimen_ID_tumor)

palette <- read_tsv(file.path(hist_dir, "broad_histology_cancer_group_palette.tsv")) %>%
  select(broad_histology, cancer_group, cancer_group_display) %>%
  mutate(cancer_group_display = case_when(is.na(cancer_group_display) ~ "Benign Tumor", 
                                          TRUE ~ as.character(cancer_group_display)))

# v22 OpenPBTA histologies
v22 <- read_tsv(file.path(hist_dir, "pbta-histologies.tsv"), guess_max = 3000) %>%
  filter(sample_type == "Tumor") %>%
  mutate(sample_id = case_when(sample_id == "7316-3214-A07083" ~ "7316-3214",
                               sample_id == "7316-3214-A07082" ~ "7316-3214",
                               TRUE ~ as.character(sample_id))) 
tel_df <- v22 %>%
  filter(Kids_First_Biospecimen_ID %in% telhunt$Kids_First_Biospecimen_ID) %>%
  left_join(telhunt, by = "Kids_First_Biospecimen_ID") %>%
  select(sample_id, broad_histology, cancer_group, ratio) %>%
  unique() %>%
  left_join(palette, by = c("broad_histology", "cancer_group"))

# We need to create a new scheme for labeling that shows wrapped cancer groups with `(n=X)`
tel_df <- tel_df %>%
  dplyr::count(cancer_group_display) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(
  dplyr::select(palette, dplyr::contains("cancer_group"))) %>%
  dplyr::select(cancer_group_display, n) %>%
  # Create wrapped with (n=X) factor column for cancer groups
  dplyr::mutate(cancer_group_display_n = glue::glue("{cancer_group_display} (N={n})")) %>%
  dplyr::inner_join(tel_df) %>%
  unique() %>%
  # reverse order levels to alphabetize when flipping coordinates
  mutate(cancer_group_display_n = fct_rev(cancer_group_display_n))

# plot
tiff(file.path(output_dir, "telhunt-distributions.tiff"), height = 1800, width = 4200, res = 300)
p <- ggplot(tel_df, aes(x = cancer_group_display_n, y = ratio)) +
  geom_violin() +
  geom_sina(alpha = 0.4) +
theme_Publication() +
  xlab("Histology") +
  ylab("T/N telomere content ratio") +
  coord_flip() +
  geom_hline(aes(yintercept = 1.0679, linetype = "HGAT")) +
  geom_hline(aes(yintercept = 0.9963, linetype = "non-HGAT"), colour = "gray") +
scale_linetype_manual(name = "ALT+ cutpoint", values = c(2, 2), 
                      guide = guide_legend(override.aes = list(color = c("black", "gray")))) 

p
dev.off()
