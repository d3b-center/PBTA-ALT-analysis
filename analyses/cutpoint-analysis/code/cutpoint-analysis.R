
library(cutpointr)
#all
#set.seed(2021)
#opt_cut <- cutpointr(jen_mer, telomere_ratio, phenotype, group, 
#         method = maximize_boot_metric,
#        boot_cut = 200, summary_func = mean,
#       metric = accuracy, silent = TRUE, na.rm = T)
opt_cut <- cutpointr(jen_mer, telomere_ratio, phenotype, group, metric = sum_sens_spec, 
                     tol_metric = 0.05, break_ties = c, na.rm = T)
summary(opt_cut)
pdf("~/Box Sync/D3B-share/collaborations/Kristina/cutoff-analysis/all-pbta-by-subgroup.pdf", height = 4, width = 6)
plot_metric(opt_cut)
dev.off()
opt_cut %>%
  filter(subgroup == "HGAT") %>%
  select(optimal_cutpoint, sum_sens_spec) %>% 
  unnest(cols = c(optimal_cutpoint, sum_sens_spec)) %>%
  arrange(-sum_sens_spec)

opt_cut %>%
  filter(subgroup == "non-HGAT") %>%
  select(optimal_cutpoint, sum_sens_spec) %>% 
  unnest(cols = c(optimal_cutpoint, sum_sens_spec)) %>%
  arrange(-sum_sens_spec)

#all
cp <- cutpointr(jen_mer, telomere_ratio, phenotype, 
                method = maximize_metric, metric = sum_sens_spec, na.rm = T)
summary(cp)
pdf("~/Box Sync/D3B-share/collaborations/Kristina/cutoff-analysis/all-pbta.pdf", height = 4, width = 8)
plot(cp)
dev.off()
#hgat
cp <- cutpointr(hgat, telomere_ratio, phenotype,
                method = maximize_metric, metric = sum_sens_spec, na.rm = T)
summary(cp)
pdf("~/Box Sync/D3B-share/collaborations/Kristina/cutoff-analysis/hgat.pdf", height = 4, width = 8)
plot(cp)
dev.off()
#non-hgat
cp <- cutpointr(non_hgat, telomere_ratio, phenotype, 
                method = maximize_metric, metric = sum_sens_spec, na.rm = T)
summary(cp)
pdf("~/Box Sync/D3B-share/collaborations/Kristina/cutoff-analysis/non-hgat.pdf", height = 4, width = 8)
plot(cp)
dev.off()

#cell lines
opt_cut <- cutpointr(hgat, telomere_ratio, phenotype, composition, metric = sum_sens_spec, 
                     tol_metric = 0.05, break_ties = c, na.rm = T)

opt_cut %>%
  filter(subgroup == "Derived Cell Line") %>%
  select(optimal_cutpoint, sum_sens_spec) %>% 
  unnest(cols = c(optimal_cutpoint, sum_sens_spec)) %>%
  arrange(-sum_sens_spec)

opt_cut %>%
  filter(subgroup == "Solid Tissue") %>%
  select(optimal_cutpoint, sum_sens_spec) %>% 
  unnest(cols = c(optimal_cutpoint, sum_sens_spec)) %>%
  arrange(-sum_sens_spec)
plot(opt_cut)


hgat_tumor <- hgat %>%
  filter(composition == "Solid Tissue")


cp <- cutpointr(hgat_tumor, telomere_ratio, phenotype,
                method = maximize_metric, metric = sum_sens_spec, na.rm = T)
summary(cp)
pdf("~/Box Sync/D3B-share/collaborations/Kristina/cutoff-analysis/hgat-tumor-only.pdf", height = 4, width = 8)
plot(cp)
dev.off()

