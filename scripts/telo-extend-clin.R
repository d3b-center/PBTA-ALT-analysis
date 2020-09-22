

#' *Global Variables Begin*
#' These global variables should be overwritten with inputs
#' Directory structure and file paths
#workDir <- "/home/data"
workDir <- "/Users/millerd15/CHOP/PBTA"
scriptsDir <- file.path(workDir,"Scripts")
telomereFile <- "telomere_940_ratio.tsv"
demographicFile <- "2020-01-13-pbta-histologies-v13.tsv"
extendStandardFile <- "pbta_clinical_extend_standard.tsv"
extendPolyAFile <- "pbta_clinical_extend_polya.tsv"
outputName <- "pbta_clinical_extend_telomere.tsv"

### read in clinical
print(paste(Sys.time(),"Reading Clinical File..."))
clin <- read.delim(file.path(workDir,demographicFile), header = T, sep = "\t")
### read in telomere
print(paste(Sys.time(),"Clinical File Loaded. Begin Reading Telomere File..."))
telo <- read.delim(file.path(workDir,telomereFile), header = T, sep = "\t")
### read in standard extend
print(paste(Sys.time(),"Telomere File Loaded. Begin Reading Standard EXTEND..."))
exts <- read.delim(file.path(workDir,extendStandardFile), header = T, sep = "\t")
### read in standard extend
print(paste(Sys.time(),"Telomere File Loaded. Begin Reading PolyA EXTEND..."))
extp <- read.delim(file.path(workDir,extendPolyAFile), header = T, sep = "\t")
print(paste(Sys.time(),"All Files Loaded."))

clin.telo <- merge(clin,
                   telo[,c("Kids_First_Biospecimen_ID_tumor","ratio")],
                   by.x="Kids_First_Biospecimen_ID",
                   by.y="Kids_First_Biospecimen_ID_tumor",
                   all=TRUE)
clin.telo.exts <- merge(clin.telo,
                        exts[,c("SampleID","NormEXTENDScores_FPKM")],
                        by.x="Kids_First_Biospecimen_ID",
                        by.y="SampleID",
                        all=TRUE)
clin.telo.exts.extp <- merge(clin.telo.exts,
                             extp[,c("SampleID","EXTENDScores_FPKM")],
                             by.x="Kids_First_Biospecimen_ID",
                             by.y="SampleID",
                             all=TRUE)

write.table(clin.telo.exts.extp,
            file.path(workDir,outputName),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)


