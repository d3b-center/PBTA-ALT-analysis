#' @description  Generate graphs for telomere ratio data based on demographic variables
#' @author Daniel Miller <millerd15@@email.chop.edu> (D3b)
#' @note Date: January 2020

suppressPackageStartupMessages(library("csv"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("forcats"))

#' **Global Variables**
#' Global Variables for dealing with continuous variables
contBreaks = 30
contDigs = 5
ratioCutoff = 1
#' Gloval Variables for minimum samples to be graphed
demoMin = 5
#' Global Variables that are overwritten when variables are passed
workDir = "/Users/millerd15/CHOP/PBTA/"
outputDir = "/Users/millerd15/CHOP/PBTA/results/"
teloPath = file.path(workDir,"telomere_940_ratio.tsv")
demoPath = file.path(workDir,"2020-01-13-pbta-histologies-v13.tsv")
themePath = file.path(workDir,"Scripts/theme_Publication.R")
colorPath = file.path(workDir,"Scripts/demog-color-function.R")

source(themePath)
source(colorPath)

#' **Demographic Metadata of Interest**
#' *List of demographic metadata to check*
#' Provide a list of metadata items to great graphs for
demos <- c(
  "broad_histology",
  "disease_type_new",
  "ethnicity",
  "germline_sex_estimate",
  "age_at_diagnosis_days"
)
#' *List of demographic metadata that are continuous*
#' If a metadata item in the list above is a continuous variable
#' As in a discrete grouping does not already exist (such as age)
#' Add that item to this list
continuousDOI <- c(
  "age_at_diagnosis_days"
)

#' Generate histogram for pre-binned data
#' @param df data frame with combined telomere ratio and demographic data
#' @param demo binned demographic item of interested
#' @param minSample minimum samples for the demo bin to be displayed
#' @param continuous T/F to declare if variable is continuous
#' @returns the generated plot
demographicBar<- function(df, demo, minSample=5, minRatio=1, continuous=FALSE) {
  # If the variable is continuous we need to bin it ourselves
  # The breaks and digits are set globally
  if (continuous) {
    iBin = cut(as.numeric(levels(telo.demo[,demo]))[telo.demo[,demo]],breaks=contBreaks,dig.lab = contDigs)
    df[,demo] = iBin
  }
  # get summary data of the data frame
  mu = ddply(df, demo, summarise, total=length(ratio),above_one=sum(ratio>=1),percent_above_one=100*sum(ratio>=1)/length(ratio))
  # maintain numerical order if variable is continuous
  if (continuous) {
    xOrder = demo
  } else {
    xOrder = paste0("reorder(",demo,",-percent_above_one)")
  }
  # make a plot of the values with at least the minSample
  plot <- ggplot(mu[mu$total>=minSample,],
         aes_string(
           # order the x axist from highest to lowest percentage
           x=xOrder,
           y="percent_above_one")) +
    # make it a histogram
    geom_bar(stat="identity") +
    # relabel the x and y axis
    labs(x=demo,y="% Samples T/N Telomere Conent Ratio >= 1") +
    # adjust the size and justification of the x axis values
    theme(axis.text.x= element_text(angle=90,hjust=1,vjust=.5)) +
    # also wrap the x axis values in case they are too long
    scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
    # add the total number of samples to each bar
    geom_text(aes(label=total), vjust=1.6, color="green", size=3.0) + 
    # add a horizontal dashed line for the overall average percentage
    geom_hline(yintercept = 100*sum(mu$above_one)/sum(mu$total),color="red",linetype="dashed")
  return(plot)
}

#' Generate Box plot for given demographic value
#' @param df The data frame with combined telomere ratio and demographic data
#' @param xaxis Primary x axis graphic variable
#' @param fill Split the primary xaxis variable by this group
#' @param facet Variable to facet the primary x axis variable
#' @param violin TRUE/FALSE Should this be a violin plot
#' @param color Colors for the primary variable
#' @return the generated plot
demographicBox <- function(df, xaxis, fill, facet, order=FALSE, violin=FALSE, color) {
  #' If a fill value is not provided set it to the xaxis value
  if (missing(fill)) {
    fill <- xaxis
  }
  if (order) {
    xlabel <- xaxis
    xaxis <- paste0("fct_reorder(",xaxis,",ratio,.desc=TRUE)")
  }
  #' Set up the basic plot with xaxis, log of the ratio, and fill
  plot <- ggplot(df, aes_string(x=xaxis,y="log(ratio)",fill=fill)) +
    scale_fill_manual(values = color)
  #' If violin true plot a violin plot otherwise do a boxplot
  if (violin) {
    plot <- plot + geom_violin()
  } else {
    plot <- plot + geom_boxplot()
  }
  #' Fix the xlab
  if (order) {
    plot <- plot + xlab(xlabel)
  }
  #' If a facet was provided facet on that
  if (!missing(facet)) {
    plot <- plot + facet_wrap(as.formula(paste("~", facet)))
  }
  #' Apply the theme
  plot <- plot + theme_Publication()
  #' Get xaxis in order
  plot <- plot +
    # adjust the size and justification of the x axis values
    theme(axis.text.x= element_text(angle=45,hjust=1,vjust=1)) +
    # also wrap the x axis values in case they are too long
    scale_x_discrete(labels = function(x) str_wrap(x, width = 25))
  return(plot)
}

#' **Main Function**
#' Create the option list
optionList <- list(
  make_option(
    opt_str = c("-d","--demographic_file"), type = "character",
    default = NULL, help = "Path to the demographic file."
  ),
  make_option(
    opt_str = c("-f","--telomere_file"), type = "character",
    default = NULL, help = "Path to the file containing T/N telomere content ratio"
  ),
  make_option(
    opt_str = c("-o", "--output_path"), type = "character",
    default = ".", help = "Path for output files."
  )
)

#' Parse the options
opt <- parse_args(OptionParser(option_list = optionList))

#' Check that the files exist
if (!file.exists(opt$demographic_file)) {
  demo_missing <- paste("Error:", opt$demographic_file, "does not exist")
  if (!file.exists(opt$telomere_file)) {
    stop(paste(demo_missing, "\nError:", opt$telomere_file, "does not exist"))
  } else
    stop(demo_missing)
} else if (!file.exists(opt$telomere_file)) {
  stop(paste("Error:", opt$telomere_file, "does not exist"))
}

#' Reset the global variables
teloPath <- opt$telomere_file
demoPath <- opt$demographic_file
outputDir <- opt$output_path

#' Check to make sure the output directory is there
#' Create it if it isn't
if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

# Load the files
telo = read.csv(teloPath,sep="\t")
demo = read.csv(demoPath,sep="\t")                                                  
telo.demo = merge(telo,demo,by.x="Kids_First_Biospecimen_ID_tumor",by.y="Kids_First_Biospecimen_ID")

#' Loop through the demographic variables
#' Open a PDF document
#' Generate the graph
#' Print it to the graph
#' Close the device
for (i in demos) {
  fileName = paste0("pbta-telomere-ratio-",i,".pdf")
  pdf(file.path(outputDir,fileName),width=11,height=7)
  if (i %in% continuousDOI) {
    demoPlot <- demographicBar(telo.demo,i,demoMin,TRUE)
  } else {
    demoPlot <- demographicBar(telo.demo,i,demoMin,FALSE)
  }
  print(demoPlot)
  dev.off()
}


