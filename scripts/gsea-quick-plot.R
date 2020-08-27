#' Tumor mutation burden vs ratio
#' Scatter plot tumor mutation burden vs the ratio score
#' @author Daniel Miller <millerd15@@email.chop.edu

library(ggplot2)

workDir <- "/Users/millerd15/CHOP/PBTA"
resultsDir <- file.path(workDir,"results")
dataDir <- file.path(workDir,"OpenPBTA-analysis/data/release-v12-20191217")
gsvaPolyAFile <- "OpenPBTA-analysis/analyses/gene-set-enrichment-analysis/results/gsva_scores_polya.tsv"
gsvaStrandedFile <- "OpenPBTA-analysis/analyses/gene-set-enrichment-analysis/results/gsva_scores_stranded.tsv"
telomereFile <- "telomere_940_ratio.tsv"
demographicFile <- "2020-01-13-pbta-histologies-v13.tsv"

rmin <- 1

#' Load theme
pThemeFile <- file.path(workDir,"Scripts/theme_Publication.R")
source(pThemeFile)

#' Load TMB files
gsvaP <- read.delim(file.path(workDir,gsvaPolyAFile), header = T, sep = "\t")
gsvaS <- read.delim(file.path(workDir,gsvaStrandedFile), header = T, sep = "\t")
telo <- read.delim(file.path(workDir,telomereFile), header = T, sep = "\t")
clin <- read.delim(file.path(workDir,demographicFile), header = T, sep = "\t")

gsvaP.clin <- merge(gsvaP,
                    clin[c("Kids_First_Biospecimen_ID","Kids_First_Participant_ID")],
                    by="Kids_First_Biospecimen_ID")
gsvaP.telo <- merge(gsvaP.clin,
                    telo[,c("Kids_First_Participant_ID","ratio")],
                    by="Kids_First_Participant_ID")
gsvaP.telo$alt <- gsvaP.telo$ratio >= rmin
gsvaS.clin <- merge(gsvaS,
                    clin[c("Kids_First_Biospecimen_ID","Kids_First_Participant_ID")],
                    by="Kids_First_Biospecimen_ID")
gsvaS.telo <- merge(gsvaS.clin,
                    telo[,c("Kids_First_Participant_ID","ratio")],
                    by="Kids_First_Participant_ID")
gsvaS.telo$alt <- gsvaS.telo$ratio >= rmin

ggplot(gsvaP.telo,aes(x=alt,y=gsea_score,fill=alt)) + 
  geom_violin() + 
  ggtitle("GSEA PolyA Score Distribution in ALT+/-") +
  theme_Publication()

ggplot(gsvaS.telo,aes(x=alt,y=gsea_score,fill=alt)) + 
  geom_violin() + 
  ggtitle("GSEA Stranded Score Distribution in ALT+/-") +
  theme_Publication()




