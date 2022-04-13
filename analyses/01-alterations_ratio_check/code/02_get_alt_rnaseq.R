#devtools::install_github("d3b-center/annoFuse",force=TRUE)

library("annoFuse")
library(readr)
library(tidyverse)
library(reshape2)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

#read in caller results
STARFusioninputfile<-read_tsv(file.path(root_dir,"data","pbta-fusion-starfusion.tsv.gz")) %>%
  dplyr::rename(`#FusionName`= FusionName)
Arribainputfile<-read_tsv(file.path(root_dir,"data","pbta-fusion-arriba.tsv.gz")) %>%
  dplyr::rename(`#gene1`= gene1)

outputfile <- "PBTA_v17"
expressionMatrixPolya<-readRDS(file.path(root_dir,"data","pbta-gene-expression-rsem-fpkm.polya.rds"))
expressionMatrixStranded<-readRDS(file.path(root_dir,"data","pbta-gene-expression-rsem-fpkm.stranded.rds"))

cohort_plot_labels <-read_tsv(file.path(root_dir,"data","pbta-histologies.tsv"))

standardizedSTARFusion<-fusion_standardization(fusion_calls = STARFusioninputfile,caller = "STARFUSION",tumorID = STARFusioninputfile$tumor_id)
standardizedArriba<-fusion_standardization(fusion_calls = Arribainputfile,caller = "ARRIBA",tumorID = Arribainputfile$tumor_id)
head(standardizedArriba[,"annots"])

#merge standardized fusion calls
standardFusioncalls<-rbind(standardizedSTARFusion,standardizedArriba) %>% as.data.frame()
#fix IGH@ issue
standardFusioncalls$FusionName[grep("IGH@",standardFusioncalls$FusionName)]<-gsub("IGH@","IGH",standardFusioncalls$FusionName[grep("IGH@",standardFusioncalls$FusionName)])
standardFusioncalls$FusionName[grep("IGH-@",standardFusioncalls$FusionName)]<-gsub("IGH-@","IGH",standardFusioncalls$FusionName[grep("IGH-@",standardFusioncalls$FusionName)])



# General fusion QC for read support and red flags
# artifact filter using GTEx_Recurrent/DGD_paralog/Normal/BodyMap/Conjoin 
# edit to retain mitelman annotated fusion in readthrough filter
fusionQCFiltered<-fusion_filtering_QC(standardFusioncalls=standardFusioncalls,readingFrameFilter="in-frame|frameshift|other",artifactFilter="GTEx_Recurrent|DGD_PARALOGS|Normal|BodyMap",junctionReadCountFilter=1,spanningFragCountFilter=100,readthroughFilter = TRUE)


# get gene symbol and ensemble ID
expressionMatrixPolya <- cbind(expressionMatrixPolya, colsplit(expressionMatrixPolya$gene_id, pattern = '_', names = c("EnsembleID","GeneSymbol")))
expressionMatrixStranded <- cbind(expressionMatrixStranded, colsplit(expressionMatrixStranded$gene_id, pattern = '_', names = c("EnsembleID","GeneSymbol")))

#Filter Expression for Polya
fusionExpressionFilteredPolya<-expression_filter_fusion(standardFusioncalls=fusionQCFiltered,expressionMatrix=expressionMatrixPolya,expressionFilter=1)
#Filter Expression Stranded
fusionExpressionFilteredStranded<-expression_filter_fusion(standardFusioncalls=fusionQCFiltered,expressionMatrix=expressionMatrixStranded,expressionFilter=1)

# Add reference gene list containing known oncogenes, tumor suppressors, kinases, and transcription factors
geneListReferenceDataTab <- read.delim(
  system.file("extdata", "genelistreference.txt", package = "annoFuse"),
  stringsAsFactors = FALSE
)
# Add fusion list containing previously reported oncogenic fusions.
fusionReferenceDataTab <- read.delim(
  system.file("extdata", "fusionreference.txt", package = "annoFuse"),
  stringsAsFactors = FALSE
)

# bioMartDataPfam <- readRDS(system.file("extdata", "pfamDataBioMart.RDS", package = "annoFuse"))
# kinaseid<-unique(bioMartDataPfam$pfam_id[grep("kinase",bioMartDataPfam$NAME)] )

# Annotate polya
filteredFusionAnnotatedPolya<-annotate_fusion_calls(standardFusioncalls=fusionExpressionFilteredPolya,geneListReferenceDataTab=geneListReferenceDataTab,fusionReferenceDataTab=fusionReferenceDataTab,checkReciprocal = TRUE)

# Annotate stranded
filteredFusionAnnotatedStranded<-annotate_fusion_calls(standardFusioncalls=fusionExpressionFilteredStranded,geneListReferenceDataTab=geneListReferenceDataTab,fusionReferenceDataTab=fusionReferenceDataTab,checkReciprocal = TRUE)


QCGeneFiltered_filtFusion<-rbind(filteredFusionAnnotatedStranded,filteredFusionAnnotatedPolya)

# subset for recurrent fusion detection and multifused genes QC
fusion_calls<-unique(QCGeneFiltered_filtFusion)

# get Telomere-length maintenance and DNA damage repair pfam id
bioMartDataPfam<-readRDS("~/Downloads/pfamDataBioMart.RDS")
id <- "PF11640"

# get putative oncogene fusions 
putative_driver_fusions <-fusion_driver(standardFusioncalls=fusion_calls,annotated=TRUE,geneListReferenceDataTab = geneListReferenceDataTab ,fusionReferenceDataTab = fusionReferenceDataTab,checkDomainStatus = TRUE,domainsToCheck = id) 

# filter for ATRX,DAXX and TERT
fus_goi <- putative_driver_fusions %>% dplyr::filter(grepl("ATRX|TERT|DAXX",FusionName)) %>%
  dplyr::select(Sample,FusionName) %>%
  unique() %>%
  group_by(Sample) %>%
  summarise(FusionName=toString(FusionName))


exp_stranded <-expressionMatrixStranded %>% 
  # Keep the gene symbols and the samples themselves
  dplyr::select(-one_of("gene_id", "EnsembleID")) %>%
  dplyr::filter(grepl("ATRX|DAXX|TERT",GeneSymbol)) %>%
  column_to_rownames("GeneSymbol") %>%
  t() %>% as.data.frame() %>%
  dplyr::rename(ATRX_fpkm=ATRX, TERT_fpkm=TERT,  DAXX_fpkm=DAXX)

exp_polya <-expressionMatrixPolya %>% 
  # Keep the gene symbols and the samples themselves
  dplyr::select(-one_of("gene_id", "EnsembleID")) %>%
  dplyr::filter(grepl("ATRX|DAXX|TERT",GeneSymbol)) %>%
  column_to_rownames("GeneSymbol") %>%
  t() %>% as.data.frame() %>%
  dplyr::rename(ATRX_fpkm=ATRX, TERT_fpkm=TERT,  DAXX_fpkm=DAXX) 

exp_goi <- rbind(exp_stranded,exp_polya) %>% rownames_to_column()

rna_alt <- exp_goi %>% left_join(fus_goi,by=c("rowname"="Sample"))  %>%
  dplyr::rename(Kids_First_Biospecimen_ID=rowname) %>%
  dplyr::mutate(ATRX_fusion = case_when(grepl("ATRX",FusionName) ~ FusionName,
                                        TRUE ~ NA_character_ ),
                DAXX_fusion = case_when(grepl("DAXX",FusionName) ~ FusionName,
                                        TRUE ~ NA_character_ ),
                TERT_fusion = case_when(grepl("TERT",FusionName) ~ FusionName,
                                        TRUE ~ NA_character_ )) %>%
  # remove column name FusionName
  dplyr::select(-FusionName)

cohort_plot_labels %>%
  dplyr::filter(experimental_strategy %in% c("RNA-Seq")) %>%
  dplyr::select("Kids_First_Biospecimen_ID","Kids_First_Participant_ID","sample_id") %>%
  left_join(rna_alt) %>% 
  write_tsv(file.path(root_dir,"analyses/01-alterations_ratio_check","output","PutativeDriver_ATRTX_DAXX_TERT_RNA_alt.tsv"))

