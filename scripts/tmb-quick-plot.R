#' Tumor mutation burden vs ratio
#' Scatter plot tumor mutation burden vs the ratio score
#' @author Daniel Miller <millerd15@@email.chop.edu

library(ggplot2)

workDir <- "/Users/millerd15/CHOP/PBTA"
resultsDir <- file.path(workDir,"results")
dataDir <- file.path(workDir,"OpenPBTA-analysis/data/release-v12-20191217")
tmbAllFileName <- "pbta-snv-consensus-mutation-tmb-all.tsv"
tmbCodingFileName <- "pbta-snv-consensus-mutation-tmb-coding.tsv"
telomereFile <- "telomere_940_ratio.tsv"

#' Load theme
pThemeFile <- file.path(workDir,"Scripts/theme_Publication.R")
source(pThemeFile)

#' Load TMB files
tmbAll <- read.delim(file.path(dataDir,tmbAllFileName), header = T, sep = "\t")
tmbCod <- read.delim(file.path(dataDir,tmbCodingFileName), header = T, sep = "\t")
telo <- read.delim(file.path(workDir,telomereFile), header = T, sep = "\t")

tmbAll.telo <- merge(tmbAll,
                     telo[,c("Kids_First_Biospecimen_ID_tumor","ratio")],
                     by.x="Tumor_Sample_Barcode",
                     by.y="Kids_First_Biospecimen_ID_tumor")
tmbCod.telo <- merge(tmbCod,
                     telo[,c("Kids_First_Biospecimen_ID_tumor","ratio")],
                     by.x="Tumor_Sample_Barcode",
                     by.y="Kids_First_Biospecimen_ID_tumor")

ggplot(tmbAll.telo,aes(x=log(tmb),y=log(ratio))) +
  geom_point() +
  ggtitle("Overall TMB vs Ratio") +
  theme_Publication()

ggplot(tmbCod.telo,aes(x=log(tmb),y=log(ratio))) +
  geom_point() +
  ggtitle("Coding TMB vs Ratio") +
  theme_Publication()




