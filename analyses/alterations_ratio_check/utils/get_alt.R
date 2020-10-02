set.seed(12345)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
setwd(root_dir)

library(data.table)
get_alt<-function(gene,format_snv="hgvs",format_cnv="length",format_sv="length"){
  
  if (format_snv=="hgvs"){
  # strelka2
  strelka2 <- format_snv_hgvs(gene,strelka2,"strelka2")
  
  # mutect2
  mutect2 <- format_snv_hgvs(gene,mutect2,"mutect2")
  
  # consensus maf
  consensus_snv <- format_snv_hgvs(gene,consensus_snv,"consensus")
  }
  
  if (format_snv=="binary"){
    # strelka2
    strelka2 <- format_snv_binary(gene,strelka2,"strelka2")
    
    # mutect2
    mutect2 <- format_snv_binary(gene,mutect2,"mutect2")
    
    # consensus maf
    consensus_snv <- format_snv_binary(gene,consensus_snv,"consensus")
  }
  

  #consensus
  consensus_cnv <- format_cnv(gene,consensus_cnv,"consensus")
  
  # cnvkit
  cnvkit <- format_cnv(gene,cnvkit,"cnvkit")
  
  # controlfreec  
  controlfreec <- format_cnv(gene,controlfreec,"controlfreec")
  
  # manta
  manta_sv <- format_sv(gene,manta_sv,caller = "manta")

  
  gene_alt_list<-list("strelka2"=strelka2,"mutect2"=mutect2,"consensus_snv"=consensus_snv,"consensus_cnv"=consensus_cnv,"controlfreec"=controlfreec,"cnvkit"=cnvkit,"manta_sv"=manta_sv)
  
}


format_snv_binary<-function(gene,alt_caller,caller) {
  colname<-paste0(gene,"_MOD_HIGH_mut_",caller)
  alt_caller <- alt_caller <- alt_caller %>% unique() %>%
    # select only MODERATE/HIGH IMPACT mutations
    dplyr::filter((IMPACT %in% c("MODERATE","HIGH")))  %>%
    dplyr::select (c("Tumor_Sample_Barcode",
                     "Hugo_Symbol","HGVSp_Short")) %>%
    # select only gene
    dplyr::filter(Hugo_Symbol == gene) %>% 
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
  alt_caller <- alt_caller %>% unique() %>%
    # select only MODERATE/HIGH IMPACT mutations
    dplyr::filter((IMPACT %in% c("MODERATE","HIGH")))  %>%
    dplyr::select (c("Tumor_Sample_Barcode",
                     "Hugo_Symbol","HGVSp_Short")) %>%
    # select only gene
    dplyr::filter(Hugo_Symbol == gene) %>% 
    # unique so that 1 = mutation exists 0 mutation doesn't exist
    unique() %>% 
    dplyr::mutate(data_source = caller ) %>%
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
  return(alt_caller)
}

format_cnv<-function(gene,alt_caller,caller) {
  colname_gain <- paste0(gene,"_gain_",caller)
  colname_loss <- paste0(gene,"_loss_",caller)
  colname_amp <- paste0(gene,"_amplification_",caller)
  colname_del <- paste0(gene,"_deletion_",caller)
  alt_caller <- alt_caller %>% 
    # select only gene 
    dplyr::filter(gene_symbol == gene) %>% 
    # unique 
    unique() 
    
  if (nrow(alt_caller)>0) {
    alt_caller <- alt_caller %>%
      # reshape to make gain /loss columns
      reshape2::dcast(biospecimen_id ~ status ,value.var = "status",fun.aggregate = length) %>%
      # add column Hugo_Symbol to use for plotting
      dplyr::mutate(Hugo_Symbol = gene)
    if (!is.null(alt_caller$gain)){
      # rename columns
      alt_caller<-dplyr::rename(alt_caller,!! colname_gain :="gain")
    }
    if (!is.null(alt_caller$loss)){
      alt_caller<-dplyr::rename(alt_caller,!! colname_loss :="loss")
    }
    if(!is.null(alt_caller$amplification)){
      alt_caller<-dplyr::rename(alt_caller,!! colname_amp :="amplification")
    }
    if(!is.null(alt_caller$deletion)){
      alt_caller<-dplyr::rename(alt_caller,!! colname_del :="deletion")
    }
  }
  return(alt_caller)
}


format_sv <- function(gene,alt_caller,caller){
  colname_BND <- paste0(gene,"_BND_",caller)
  colname_DEL <- paste0(gene,"_DEL_",caller)
  colname_DUP <- paste0(gene,"_DUP_",caller)
  colname_INV <- paste0(gene,"_INV_",caller)
  
  alt_caller <- alt_caller %>%
    reshape2::dcast(Kids.First.Biospecimen.ID.Tumor ~ SV.type,fun.aggregate = length) %>%
    # add Hugo_Symbol column for plotting
    dplyr::mutate(Hugo_Symbol = gene) 
  
  if (!is.null(alt_caller$BND)){
    # rename columns
    alt_caller<-dplyr::rename(alt_caller,!! colname_BND :="BND")
  }
  if (!is.null(alt_caller$DEL)){
    alt_caller<-dplyr::rename(alt_caller,!! colname_DEL :="DEL")
  }
  if(!is.null(alt_caller$DUP)){
    alt_caller<-dplyr::rename(alt_caller,!! colname_DUP :="DUP")
  }
  if(!is.null(alt_caller$INV)){
    alt_caller<-dplyr::rename(alt_caller,!! colname_INV :="INV")
  }
  return(alt_caller)

  
}


