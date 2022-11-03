# Libraries
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(plyr)

source()

vaf_alt <- atrx_mut %>%
  left_join(telhunt) %>%
  left_join(ihc) %>%
  filter(!is.na(ATRXm_VAF)) %>%
  filter(!(is.na(`T/N TelHunt ratio`) & is.na(`C-Circle Assay`) & is.na(UBTF))) %>%
  # filter duplicate rows if VUS + oncogenic
  mutate(remove = case_when(sample_id == "7316-2189" & grepl("VUS", ATRXm) ~ "remove",
                     sample_id == "7316-2756" & grepl("VUS", ATRXm) ~ "remove",
                     sample_id == "7316-4723" & grepl("VUS", ATRXm) ~ "remove",
                     sample_id == "7316-5159" & grepl("VUS", ATRXm) ~ "remove",
                     ATRXm == "p.E929Qfs*8",
                     TRUE ~ "keep")) %>%
  filter(remove == "keep") %>%
  mutate(ALT_status = case_when(`T/N TelHunt ratio` > 1.699 | `C-Circle Assay` == "POS" | UBTF == "POS" ~ "POS",
                                `T/N TelHunt ratio` < 1.699 | `C-Circle Assay` == "NEG" & UBTF != "POS" ~ "NEG",
                                `T/N TelHunt ratio` < 1.699 | `C-Circle Assay` != "POS" & UBTF == "NEG" ~ "NEG",
                                TRUE ~ NA_character_)) %>%
  select(sample_id, ATRXm, ATRXm_VAF, `T/N TelHunt ratio`, `C-Circle Assay`, `UBTF`, ALT_status)

mu <- plyr::ddply(vaf_alt, "ALT_status", summarise, grp.mean=mean(ATRXm_VAF))


# With transparency (right)
ggplot(data=vaf_alt, aes(x=ATRXm_VAF, group=ALT_status, fill=ALT_status, after_stat(count))) +
  geom_density(adjust=1.5, alpha=.4) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=ALT_status),
             linetype="dashed") +
  theme_minimal()
