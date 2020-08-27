#' Create survival plots for ALT data
#' Overall and within histologies
#' @author Daniel Miller <millerd15@@email.chop.edu>

workDir <- "/Users/millerd15/CHOP/PBTA"
resultsDir <- file.path(workDir,"results")
analysisFile <- file.path(workDir,"OpenPBTA-analysis/analyses/survival-analysis/util/survival_models.R")
clinteloFile <- file.path(workDir,"Scripts/telo-extend-clin.R")
pThemeFile <- file.path(workDir,"Scripts/theme_Publication.R")
rmin = 1

# Get the survival analysis script and its dependencies
source(analysisFile)
# Get your clin telo file made
source(clinteloFile)
source(pThemeFile)
library("stringr")

clin.telo$alt <- clin.telo$ratio >= rmin

for (i in unique(na.omit(clin.telo$disease_type_new))) {
  print(i)
  tryCatch({
    kap_fit <- survival_analysis(clin.telo[clin.telo$disease_type_new == i,],
                                 ind_var = "alt"
    )
    surv_plot <- survminer::ggsurvplot(kap_fit$model,
                                       pval = TRUE,
                                       data = kap_fit$original_data,
                                       risk.table = TRUE,
                                       xlim = c(0, 2000),
                                       break.time.by = 500,
                                       ggtheme = theme_Publication(),
                                       title = paste(i,"Survival"),
                                       risk.table.y.text.col = TRUE,
                                       risk.table.y.text = FALSE
    )
    surv_plot <- cowplot::plot_grid(surv_plot[[1]], surv_plot[[2]], nrow = 2, 
                                    rel_heights = c(2.5, 1))
    outputName <- paste0("survival-alt-",str_replace_all(tolower(i)," ","-"),".pdf")
    pdf(file.path(resultsDir,outputName),height = 8.5,width = 11)
    print(surv_plot)
    dev.off()
  }, error=function(e) {
    cat("ERROR :",conditionMessage(e), "\n")
  })
}

  
  
  
  
  