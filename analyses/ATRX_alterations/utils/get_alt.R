set.seed(12345)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

library(data.table)
get_alt<-function(gene){
  
  strelka2 <- strelka2 %>% unique() %>%
    # select only MODERATE/HIGH IMPACT mutations
    dplyr::filter((IMPACT %in% c("MODERATE","HIGH")))  %>%
    dplyr::select (c("Tumor_Sample_Barcode",
                     "Hugo_Symbol","HGVSp_Short")) %>%
    # select only gene
    dplyr::filter(Hugo_Symbol %in% gene) %>% 
    # unique so that 1 = mutation exists 0 mutation doesn't exist
    unique() %>% 
    dplyr::mutate(data_source="strelka2")
  
  mutect2 <- mutect2 %>% unique() %>%
    # select only MODERATE/HIGH IMPACT mutations
    dplyr::filter((IMPACT %in% c("MODERATE","HIGH")))  %>%
    dplyr::select (c("Tumor_Sample_Barcode",
                     "Hugo_Symbol","HGVSp_Short")) %>%
    # select only gene
    dplyr::filter(Hugo_Symbol %in%  gene) %>% 
    # unique so that 1 = mutation exists 0 mutation doesn't exist
    unique() %>% 
    dplyr::mutate(data_source="mutect2")
  
  # consensus maf
  consensus_snv <- consensus_snv %>%
    # select only MODERATE/HIGH IMPACT mutations
    dplyr::filter((IMPACT %in% c("MODERATE","HIGH")))  %>%
    dplyr::select (c("Tumor_Sample_Barcode",
                     "Hugo_Symbol","HGVSp_Short")) %>%
    # select only gene
    dplyr::filter(Hugo_Symbol %in%  gene ) %>% 
    # unique so that 1 = mutation exists 0 mutation doesn't exist
    unique() %>% 
    dplyr::mutate(data_source="consensus_snv")
 
  # all_cnv<- read_tsv(file.path(root_dir,
  #   "data",
  #   "pbta-cnv-consensus-gistic",
  #   "all_data_by_genes.txt"),col_names = TRUE) %>% 
  #   # select TP53 in all_data_by_genes.txt
  #   dplyr::filter(`Gene Symbol`=="TP53" ) %>%
  #   dplyr::select(-`Gene ID`,-Cytoband ) %>% 
  #   column_to_rownames("Gene Symbol") %>% t() %>% 
  #   as.data.frame() %>% rownames_to_column() %>% 
  #   # rename
  #   dplyr::rename("TP53_all_gistic"="TP53")
  
  
  consensus_cnv <- read_tsv(file.path(root_dir,
                                      "data",
                                      "consensus_seg_annotated_cn_autosomes.tsv.gz")) %>%
    # select TP53
    dplyr::filter(gene_symbol %in% gene) %>% 
    dplyr::mutate(data_source="consensus_cnv")
  
  gene_alt_list<-list("strelka2"=strelka2,"mutect2"=mutect2,"consensus_snv"=consensus_snv,"consensus_cnv"=consensus_cnv)
  
}
