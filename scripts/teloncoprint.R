#' @title Teloncoprinter
#' @description Print custom oncoplots that have telomere ratios in them
#' @author Daniel Miller <millerd15@@mail.chop.edu>

#' Download a package if it's not installed
#' @param p the package name
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

#' Get your packages
usePackage("BiocManager")
usePackage("tidyr")
usePackage("dplyr")
usePackage("circlize")
usePackage("gplots")
usePackage("RColorBrewer")

#' Get maftools from BiocManager if not installed
if (!requireNamespace("maftools", quietly = TRUE))
  BiocManager::install("maftools", ask = FALSE)

library("maftools")

#' *Global Variables Begin*
#' These global variables should be overwritten with inputs
#' Directory structure and file paths
#workDir <- "/home/data"
workDir <- "/Users/millerd15/CHOP/PBTA"
scriptsDir <- file.path(workDir,"Scripts")
dataDir <- file.path(workDir,"OpenPBTA-analysis/data/release-v12-20191217")
resultsDir <- file.path(workDir,"results")
fusionsFile <- "pbta-fusion-putative-oncogenic.tsv"
mafFile <- "pbta-snv-consensus-mutation.maf.tsv.gz"
demographicFile <- "pbta-histologies.tsv"
telomereFile <- "telomere_940_ratio.tsv"
goiFile <- "brain-goi-list-new.txt"
mutColorFile <- "mutation-color-function.R"
demColorFile <- "demog-color-function.R"
gisticAllLesionsPath <- file.path(dataDir,"2019-12-10-gistic-results-cnvkit/all_lesions.conf_90.txt")
gisticAmpGenesPath <- file.path(dataDir,"2019-12-10-gistic-results-cnvkit/amp_genes.conf_90.txt")
gisticDelGenesPath <- file.path(dataDir,"2019-12-10-gistic-results-cnvkit/del_genes.conf_90.txt")
gisticScoresPath <- file.path(dataDir,"2019-12-10-gistic-results-cnvkit/scores.gistic")
outputName <- "OpenPBTA-gistic-telo-oncoprint-top25-logTMB.pdf"
#' Numeric values
nskip <- 0
maximumGenes <- 25
minimumRatio <- 1
#' The values for what we will be using to annotate the oncoplot
clinAnnoValues <- c(
  "Kids_First_Biospecimen_ID",
  "germline_sex_estimate",
  "ethnicity",
  "assay",
  "tumor_descriptor",
  "disease_type_new",
  "alt"
  )
clinAnnoNick <- c(
  "Tumor_Sample_Barcode", 
  "sex",
  "ethnicity",
  "assay",
  "phase_of_therapy",
  "disease_type_new",
  "alt"
  )
oncoAnnotations <- c(
  "alt",
  "disease_type_new",
  "sex",
  "ethnicity",
  "phase_of_therapy"
  )
#' *End global variables*

#' source colors
source(file.path(scriptsDir,mutColorFile))
colors <- c("Amp"="red","Del"="blue")
source(file.path(scriptsDir,demColorFile))

#' Reformat a fusions file to work with a MAF
#' Not really safe to change
#' @param df The fusions file data frame
#' @return The transformed data frame
fusionsForMaf <- function(df) {
  # First rename the Sample column to match MAF label
  colnames(df)[colnames(df)=="Sample"] <- "Kids_First_Biospecimen_ID"
  # Next split the FusionName column into a 5' and 3' gene column
  df.sep <- df %>% tidyr::separate(FusionName, c("5'-gene", "3'-gene"), sep = "--")
  # Make a smaller table of just the IDs and genes
  df.short5 <- df.sep[,c("Kids_First_Biospecimen_ID", "5'-gene", "3'-gene")]
  # Add a fourth column that has the first and second col values concat
  df.short5$ID <- paste0(df.short5$Kids_First_Biospecimen_ID, ";", df.short5$`5'-gene`)
  # Make a summary table of the ID and the total counts in the table
  reformat5 <- df.short5 %>% 
    group_by(ID) %>% 
    dplyr::summarise(Kids_First_Biospecimen_ID = n()) %>%
    as.data.frame()
  # Repeat the last three steps for the 3 prime genes
  df.short3 <- df.sep[,c("Kids_First_Biospecimen_ID", "5'-gene", "3'-gene")]
  df.short3$ID <- paste0(df.short3$Kids_First_Biospecimen_ID, ";", df.short3$`3'-gene`)
  reformat3 <- df.short3 %>% 
    group_by(ID) %>% 
    dplyr::summarise(Kids_First_Biospecimen_ID = n()) %>%
    as.data.frame()
  # Merge the tables
  df.new <- rbind(reformat5,reformat3)
  # Add a new column to denote Multi_Hit_Fusion if total over 1
  df.new$Variant_Classification <- ifelse(df.new$Kids_First_Biospecimen_ID == 1, "Fusion", "Multi_Hit_Fusion")
  # Now split back out the ID into its components
  df.new <- df.new %>% separate(ID, c("Kids_First_Biospecimen_ID", "Hugo_Symbol"), sep = ";")
  # Add a placeholder column for variant type
  df.new$Variant_Type <- "OTHER"
  # return the new data frame
  return(df.new)
}

#' Generate a clinical frame to accompany the raw MAF in generating the oncoplot
#' Make changes here to alter the stuff the ends up in the clinical tract in the oncoplot
#' @param clin The raw clinical file
#' @param anno List of column names from the clinical frame to be used as oncoplot annotations
#' @param nick List of nicknames for the anno list
#' @return newclin The new clinical data frame
clinForMaf <- function(clin, anno, nick) {
  # Add a TSB to aggregate values
  clin$Tumor_Sample_Barcode <- paste0(clin$Kids_First_Participant_ID, ";", clin$sample_id)
  # add assay value to clin
  rna <- subset(clin, experimental_strategy == "RNA-Seq")
  dna <- subset(clin, experimental_strategy == "WGS")
  # Get those barcodes that are only in the RNA table
  rnaonly <- setdiff(rna$Tumor_Sample_Barcode, dna$Tumor_Sample_Barcode)
  # Do the same for the DNA table
  dnaonly <- setdiff(dna$Tumor_Sample_Barcode, rna$Tumor_Sample_Barcode)
  # Add a new column to clin to denote if the assay used
  clin$assay <- ifelse(clin$Tumor_Sample_Barcode %in% rnaonly, "RNA-Seq",
                       ifelse(clin$Tumor_Sample_Barcode %in% dnaonly, "WGS", "Both"))
  # Grab the desired columns in the list
  newclin <- unique(clin[,anno])
  # Rename the columns
  names(newclin) <- nick
  # And return
  return(newclin)
}

#' Add telomere ratios to the clinical file
#' @param clin The raw clinical df
#' @param telo The raw telomere df
#' @param rmin The minimum ratio value to qualify the sample as ALT
#' @return the modified clinical df including telomere ratios
makeClinTelo <- function(clin,telo,rmin) {
  clin.telo <- merge(clin,
                     telo[,c("Kids_First_Biospecimen_ID_tumor","ratio")],
                     by.x="Kids_First_Biospecimen_ID",
                     by.y="Kids_First_Biospecimen_ID_tumor",
                     all.x = TRUE)
  clin.telo$alt <- clin.telo$ratio >= rmin
  return(clin.telo)
}

#' Build the MAF for plotting
#' @param fuse Fusions file formatted for MAF
#' @param rmaf Raw MAF df
#' @param clin Raw clinical df
#' @param cloi Modified clinical file for MAF
#' @return a maf with all the accompanying files integrated
buildOncoMaf <- function(rmaf,fuse,clin,cloi) {
  # Change the column name so we can match the maf to the fusion frame
  colnames(rmaf)[colnames(rmaf)=="Tumor_Sample_Barcode"] <- "Kids_First_Biospecimen_ID"
  # Add the fusions to the MAF file
  maf.fus <- bind_rows(rmaf, fuse)
  # Cut down on space by killing the old maf
  rm(rmaf)
  # Create TSB from PT_id and sample_id
  #clin$Tumor_Sample_Barcode <- paste0(clin$Kids_First_Participant_ID, ";", clin$sample_id)
  # Add the TSB to fusion maf where matching the KF ID for matching later
  #maf.tsb <- merge(maf.fus,
  #                 clin[,c("Kids_First_Biospecimen_ID", "Tumor_Sample_Barcode")], 
  #                 all.x = T)
  # Add back the TSB to the fram
  maf.fus$Tumor_Sample_Barcode <- maf.fus$Kids_First_Biospecimen_ID
  maf.tsb <- maf.fus
  # again delete the old MAF
  rm(maf.fus)
  #read maf
  maf = read.maf(maf = maf.tsb,
                 clinicalData = cloi,
                 removeDuplicatedVariants = F,
                 gisticAllLesionsFile = gisticAllLesionsPath,
                 gisticAmpGenesFile = gisticAmpGenesPath,
                 gisticDelGenesFile = gisticDelGenesPath,
                 gisticScoresFile = gisticScoresPath,
                 vc_nonSyn = c("Frame_Shift_Del","Frame_Shift_Ins",
                               "Splice_Site", "Splice_Region",
                               "Translation_Start_Site","Nonsense_Mutation", 
                               "Nonstop_Mutation", "In_Frame_Del",
                               "In_Frame_Ins", "Missense_Mutation",  
                               "Stop_Codon_Ins", "Start_Codon_Del",
                               "Fusion", "Multi_Hit_Fusion"))
  # return the oncomaf
  return(maf)
}

#' Create list of top genes in the final MAF
#' @param fmaf Final oncomaf
#' @param gene Raw gene data fram3
#' @param maxg Maximum number of genes to return
#' @return Data frame containing the top hit genes
topGenesForMaf <- function(fmaf,gene,maxg) {
  #' Subset the MAF by genes of interest
  sub.maf.goi = maftools::subsetMaf(maf = fmaf, genes = gene$V1, mafObj = TRUE)
  gene.sum <- mafSummary(sub.maf.goi)$gene.summary
  #' remove silent/intronic variants and recalculate
  gene.sum$Intron <- NULL
  gene.sum$Silent <- NULL
  gene.sum$RNA <- NULL
  gene.sum$`5'Flank` <- NULL
  gene.sum$`3'Flank` <- NULL
  gene.sum$`3'UTR` <- NULL
  gene.sum$`5'UTR` <- NULL
  #' new total sum
  gene.sum$total <- rowSums(gene.sum[,2:12])
  #' order by alteration totals and choose top N as noted below
  goi.ordered <- gene.sum[order(gene.sum$total, decreasing = T),]
  #' Determine which is larger, the max gene value or the number of rows
  N <- ifelse(nrow(gene.sum) < maxg, nrow(gene.sum), maxg)
  #' Take the smaller number maxg or the number of rows
  goi.ordered.N <- goi.ordered[1:N,]
  return(goi.ordered.N)
}

#' Make colors for your unique variables in a frame
#' @param data The data frame containing the variable
#' @param var The varaible to which we will assign colors
#' @return The color list
colorVars <- function(data,var) {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  values <- unique(data[var])[[1]]
  lenvar <- length(values)
  ncols <- sample(col_vector,lenvar)
  names(ncols) <- values
  return(ncols)
}

#' **MAIN SCRIPT**
### read in fusions
print(paste(Sys.time(),"Reading Fusion File..."))
fuse <-read.delim(file.path(dataDir,fusionsFile), header = T, sep = "\t")
### read in clinical
print(paste(Sys.time(),"Fusion File Loaded. Begin Reading Clinical File..."))
clin <- read.delim(file.path(dataDir,demographicFile), header = T, sep = "\t")
### read in telomere
print(paste(Sys.time(),"Clinical File Loaded. Begin Reading Telomere File..."))
telo <- read.delim(file.path(workDir,telomereFile), header = T, sep = "\t")
### read in MAF
print(paste(Sys.time(),"Telomere File Loaded. Begin Reading MAF..."))
#' Check to see how many lines you need to skip off the top of your MAF 
rmaf <- read.delim(file.path(dataDir,mafFile), skip = nskip, sep = "\t", header = T)
### read in genes of interst
print(paste(Sys.time(),"MAF File Loaded. Reading Genes of Interest..."))
gene <- read.delim(file.path(workDir,goiFile), header = F)
print(paste(Sys.time(),"All Files Loaded"))

# Make a new fusion frame to add to the MAF
new.fus <- fusionsForMaf(fuse)
# Make a new clinical frame with telomere values
clin.telo <- makeClinTelo(clin,telo,minimumRatio)
# Make a new clinical frame to accompany the MAF for oncoploting
new.clin <- clinForMaf(clin.telo,clinAnnoValues,clinAnnoNick)
# Create the new MAF
oncomaf <- buildOncoMaf(rmaf,new.fus,clin,new.clin)
# Get the top genes for the final MAF
topGenes <- topGenesForMaf(oncomaf,gene,maximumGenes)

#' Get your colors for your variables
altcol <- colorVars(new.clin,"alt")
newcol <- colorVars(new.clin,"disease_type_new")

# Make the plot
pdf(file.path(resultsDir,outputName), height = 17, width = 22)
print(maftools::oncoplot(maf = oncomaf,
                         genes = topGenes$Hugo_Symbol, 
                         removeNonMutated = T, 
                         showTumorSampleBarcodes = F,
                         drawRowBar = T,
                         bgCol = "whitesmoke",
                         sepwd_genes = 1, sepwd_samples = 1,
                         annotationFontSize = 0.8,
                         legendFontSize = 1,
                         colors = c(colors,colores),
                         gene_mar = 7, barcode_mar = 6,
                         sortByAnnotation = T, fontSize = 1,
                         writeMatrix = T, showTitle = F, logColBar = T,
                         clinicalFeatures = oncoAnnotations,
                         annotationColor = list(alt = altcol,
                                                disease_type_new = newcol, 
                                                phase_of_therapy = phasecol,
                                                sex = sexcol, 
                                                ethnicity = ethcol)))
#legendFontSize = 2,

#age = colorRamp2(c(0, 4500), 
#c("#4F94CD","#171717")))))
#                       tp53_score_discrete = tp53_score_discrete,
#                      nf1_score_discrete = nf1_score_discrete,
#                     Molecular.Subtype = subcol)))

dev.off()

getwd()

#' *Old Code no longer in use but may be useful*
###subset only LGG
#as.data.frame(table(clin$disease_type_new)) #use Low-grade astrocytic tumor
#lgg <- subset(clin, disease_type_new == "Low-grade astrocytic tumor")
#length(unique(lgg$Tumor_Sample_Barcode))
#lgg.maf <- maf.tsb[maf.tsb$Kids_First_Biospecimen_ID %in% lgg$Kids_First_Biospecimen_ID,]
#names(lgg.maf)
#colnames(lgg.maf)[colnames(lgg.maf)== "Tumor_Sample_Barcode"] <- "Kids_First_Participant_ID"
#lgg.maf <- merge(lgg.maf, clin[,c("Kids_First_Biospecimen_ID", "Tumor_Sample_Barcode")], all.x = T)
#length(unique(lgg.maf$Tumor_Sample_Barcode))
#table(lgg.maf$Variant_Classification)
#




