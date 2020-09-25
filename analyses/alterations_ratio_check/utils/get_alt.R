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
  
  
  consensus_cnv <- consensus_cnv %>%
    # select TP53
    dplyr::filter(gene_symbol %in% gene) %>% 
    dplyr::mutate(data_source="consensus_cnv")
  
  manta_sv <- manta_sv %>% 
    dplyr::filter(grepl(gene,Gene.name),FILTER=="PASS") 
  
  cnvkit <- cnvkit %>%
    # select TP53
    dplyr::filter(gene_symbol %in% gene) %>% 
    dplyr::mutate(data_source="consensus_cnv")
    
  controlfreec <- controlfreec %>%
    # select TP53
    dplyr::filter(gene_symbol %in% gene) %>% 
    dplyr::mutate(data_source="consensus_cnv")
  
  gene_alt_list<-list("strelka2"=strelka2,"mutect2"=mutect2,"consensus_snv"=consensus_snv,"consensus_cnv"=consensus_cnv,"controlfreec"=controlfreec,"cnvkit"=cnvkit,"manta_sv"=manta_sv)
  
}


format_snv_binary<-function(gene,alt_caller,caller) {
  colname<-paste0(gene,"_MOD_HIGH_mut_",caller)
  alt_caller <- alt_caller %>% 
    dplyr::select(c("Tumor_Sample_Barcode",
                    "Hugo_Symbol")) %>%
    # select only TP53
    dplyr::filter(Hugo_Symbol %in% gene ) %>% 
    # unique so that 1 = mutation exists 0 mutation doesn't exist
    unique()
  if (nrow(alt_caller)>0) {
    alt_caller <- alt_caller %>%
      # reshape to make TP53 row
      reshape2::dcast(Tumor_Sample_Barcode ~ Hugo_Symbol,fun.aggregate = length) %>%
      # rename 
      dplyr::rename(!!colname := gene )
  }
}

format_snv_hgvs<-function(gene,alt_caller,caller) {
  colname<-paste0(gene,"_MOD_HIGH_mut_",caller)
  alt_caller <- alt_caller %>%
    dplyr::select (c("Tumor_Sample_Barcode",
                     "Hugo_Symbol","HGVSp_Short")) %>%
    # select only gene 
    dplyr::filter(Hugo_Symbol== gene) %>% 
    # unique so that 1 = mutation exists 0 mutation doesn't exist
    unique() 
  if (nrow(alt_caller)>0) {
    alt_caller <- alt_caller %>%
      dplyr::group_by(Tumor_Sample_Barcode ,Hugo_Symbol) %>%
      # aggregate mutations
      summarise( !! colname := toString(HGVSp_Short))
  }
}

format_cnv<-function(gene,alt_caller,caller) {
  colname_gain <- paste0(gene,"_gain_consensus_cnv_",caller)
  colname_loss <- paste0(gene,"_loss_consensus_cnv_",caller)
  alt_caller <- alt_caller %>% 
    # select only gene 
    dplyr::filter(gene_symbol== gene) %>% 
    # unique 
    unique() 
    
  if (nrow(alt_caller)>0) {
    alt_caller <- alt_caller %>%
      # reshape to make gain /loss columns
      reshape2::dcast(biospecimen_id ~ status ,value.var = "status",fun.aggregate = length) %>%
      # rename columns
      dplyr::rename(!! colname_gain :="gain",!! colname_loss :="loss") %>%
      # add column Hugo_Symbol to use for plotting
      dplyr::mutate(Hugo_Symbol = gene)
  }
}
