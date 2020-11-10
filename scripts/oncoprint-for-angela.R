#install.packages("usethis")
#usethis::edit_r_environ()
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",dep = TRUE)

if (!requireNamespace("maftools", quietly = TRUE))
  BiocManager::install("maftools", ask = FALSE)

usePackage("tidyr")
usePackage("dplyr")
usePackage("circlize")
usePackage("gplots")


# Global Variables
workDir <- "/Users/millerd15/CHOP/PBTA"
scriptsDir <- file.path(workDir,"Scripts")
dataDir <- file.path(workDir,"OpenPBTA-analysis/data/release-v12-20191217")
resultsDir <- file.path(workDir,"results")
fusionsFile <- "pbta-fusion-putative-oncogenic.tsv"
mafFile <- "pbta-snv-consensus-mutation.maf.tsv.gz"
demographicFile <- "pbta-histologies.tsv"
goiFile <- "brain-goi-list-new.txt"
mutColorFile <- "mutation-color-function.R"
demColorFile <- "demog-color-function.R"

###source colors
source(file.path(scriptsDir,mutColorFile))
source(file.path(scriptsDir,demColorFile))

###reformat fusions for MAF
fus <-read.delim(file.path(dataDir,fusionsFile), header = T, sep = "\t")
names(fus)
colnames(fus)[colnames(fus)== "Sample"] <- "Kids_First_Biospecimen_ID"
fus.sep <- fus %>% tidyr::separate(FusionName, c("5'-gene", "3'-gene"), sep = "--")
names(fus.sep)

###determine which samples have multi-fused fusions (5')
multi <- fus.sep[,c("Kids_First_Biospecimen_ID", "5'-gene", "3'-gene")]
multi$ID <- paste0(multi$Kids_First_Biospecimen_ID, ";", multi$`5'-gene`)
reformat.multi <- multi %>% 
  group_by(ID) %>% 
  dplyr::summarise(Kids_First_Biospecimen_ID = n()) %>%
  as.data.frame()

##3' multi-fused genes:
multi3 <- fus.sep[,c("Kids_First_Biospecimen_ID", "5'-gene", "3'-gene")]
multi3$ID <- paste0(multi$Kids_First_Biospecimen_ID, ";", multi$`3'-gene`)
reformat.multi3 <- multi3 %>% 
  group_by(ID) %>% 
  dplyr::summarise(Kids_First_Biospecimen_ID = n()) %>%
  as.data.frame()

###add all fused genes together
new.fus <- rbind(reformat.multi, reformat.multi3)
new.fus$Variant_Classification <- ifelse(new.fus$Kids_First_Biospecimen_ID == 1, "Fusion", "Multi_Hit_Fusion")
new.fus <- new.fus %>% separate(ID, c("Kids_First_Biospecimen_ID", "Hugo_Symbol"), sep = ";")
new.fus$Variant_Type <- "OTHER"


###read in maf
print(paste(Sys.time(),"Begin Reading MAF..."))
strelka <- read.delim(file.path(dataDir,mafFile), skip = 1, sep = "\t", header = T)
print(paste(Sys.time(),"MAF Reading Complete."))

colnames(strelka)[colnames(strelka)== "Tumor_Sample_Barcode"] <- "Kids_First_Biospecimen_ID"

###add fusion to maf
###merge with MAF
maf.fus <- bind_rows(strelka, new.fus)
rm(strelka)

###add patient ID as TSB
###read in clinical
clin <- read.delim(file.path(dataDir,demographicFile), header = T, sep = "\t")
names(clin)
##create TSB from PT_id and sample_id
clin$Tumor_Sample_Barcode <- paste0(clin$Kids_First_Participant_ID, ";", clin$sample_id)
###merge with maf
maf.tsb <- merge(maf.fus, clin[,c("Kids_First_Biospecimen_ID", "Tumor_Sample_Barcode")], all.x = T)
rm(maf.fus)

###genelist
#genes <- read.delim("~/Box Sync/D3B-share/data/genelists/brain-goi-list-new.txt", header = F)
genes <- read.delim(file.path(workDir,goiFile), header = F)

###subset only LGG
#as.data.frame(table(clin$disease_type_new)) #use Low-grade astrocytic tumor
lgg <- subset(clin, disease_type_new == "Low-grade astrocytic tumor")
length(unique(lgg$Tumor_Sample_Barcode))
lgg.maf <- maf.tsb[maf.tsb$Kids_First_Biospecimen_ID %in% lgg$Kids_First_Biospecimen_ID,]
names(lgg.maf)
#colnames(lgg.maf)[colnames(lgg.maf)== "Tumor_Sample_Barcode"] <- "Kids_First_Participant_ID"
#lgg.maf <- merge(lgg.maf, clin[,c("Kids_First_Biospecimen_ID", "Tumor_Sample_Barcode")], all.x = T)
length(unique(lgg.maf$Tumor_Sample_Barcode))

table(lgg.maf$Variant_Classification)

##fix assay
rna <- subset(lgg, experimental_strategy == "RNA-Seq")
dna <- subset(lgg, experimental_strategy == "WGS")
rnaonly <- setdiff(rna$Tumor_Sample_Barcode, dna$Tumor_Sample_Barcode)
dnaonly <- setdiff(dna$Tumor_Sample_Barcode, rna$Tumor_Sample_Barcode)

###add column for experimental_strategy
lgg$assay <- ifelse(lgg$Tumor_Sample_Barcode %in% rnaonly, "RNA-Seq",
                    ifelse(lgg$Tumor_Sample_Barcode %in% dnaonly, "WGS", "Both"))

###collapse clinical by TSB
newclin <- unique(lgg[,c("Tumor_Sample_Barcode", "glioma_brain_region", "germline_sex_estimate", "ethnicity", "assay", "tumor_descriptor")])
names(newclin) <- c("Tumor_Sample_Barcode", "brain_region", "sex", "ethnicity", "assay", "phase_of_therapy")


###read maf
maf = read.maf(maf = lgg.maf, clinicalData = newclin, removeDuplicatedVariants = F,
                                                    vc_nonSyn = c("Frame_Shift_Del","Frame_Shift_Ins", "Splice_Site", "Splice_Region", "Translation_Start_Site","Nonsense_Mutation", 
                                                                        "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation",  
                                                                        "Stop_Codon_Ins", "Start_Codon_Del", "Fusion", "Multi_Hit_Fusion"))


sub.maf.goi = maftools::subsetMaf(maf = maf, genes = genes$V1, mafObj = TRUE)
gene.sum <- mafSummary(sub.maf.goi)$gene.summary
###remove silent/intronic variants and recalculate
gene.sum$Intron <- NULL
gene.sum$Silent <- NULL
gene.sum$RNA <- NULL
gene.sum$`5'Flank` <- NULL
gene.sum$`3'Flank` <- NULL
gene.sum$`3'UTR` <- NULL
gene.sum$`5'UTR` <- NULL
##new total sum
gene.sum$total <- rowSums(gene.sum[,2:12])
names(gene.sum)

###order by alteration totals and choose top N as noted below
goi.ordered <- gene.sum[order(gene.sum$total, decreasing = T),]

###select N top genes
N <- ifelse(nrow(gene.sum) < 25, nrow(gene.sum), 25)

goi.ordered.N <- goi.ordered[1:N,]  
###generate oncoprint
prefix_name <- "lgg"

pdf(file.path(resultsDir,"OpenPBTA-lgg-oncoprint-top25-logTMB.pdf"), height = 13, width = 15)
print(maftools::oncoplot(maf = maf,  genes = goi.ordered.N$Hugo_Symbol, 
                         #top = 20, 
                         removeNonMutated = F, bgCol = "whitesmoke", 
                         showTumorSampleBarcodes = F, drawRowBar = T, sepwd_genes = 1, sepwd_samples = 1,
                         annotationFontSize = 1, gene_mar = 7, barcode_mar = 6,
                         sortByAnnotation = T, fontSize = 1, #legendFontSize = 2,
                         colors = colores, writeMatrix = T, showTitle = F, logColBar = T,
                         clinicalFeatures = c( "phase_of_therapy","brain_region", "sex", "ethnicity", "assay"),
                         annotationColor = list(phase_of_therapy = phasecol,
                                                brain_region = braincol,
                                                sex = sexcol, 
                                                ethnicity = ethcol,
                                                assay = assaycol)))
                                              #age = colorRamp2(c(0, 4500), 
                                              #c("#4F94CD","#171717")))))
                         #                       tp53_score_discrete = tp53_score_discrete,
                          #                      nf1_score_discrete = nf1_score_discrete,
                           #                     Molecular.Subtype = subcol)))

dev.off()

  getwd()
