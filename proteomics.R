#Proteomics analysis through DEPv1.21.0

##Setup and data loading
#First setup run
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DEP")

suppressPackageStartupMessages({
  library("BiocStyle")
  library("DEP")
  library("dplyr")
  library("readr")
  library("tidyverse")
})
require(DEP) #Citation 2018 publication Vermeulen; online vignette #https://bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html

require(ggplot2)
require(dplyr)
require(readr)
require(ggpubr)

##Reading the LFQ matrices
###We generated proteomics data for MResistant LNCaP in following setup:
###4x R0 (wt LNCaP:RFP-H2B)
###4x R9 clone A
###4x R9 clone B
###4x R9 clone C

#defining paths and loading data
path = ""
path_out = paste0(path, "R/")
path_plot = paste0(path_out, "plot/")

#loading the data
options(readr.num_columns = 0)
r9a <- data.frame(read.delim(paste0(path, "221124_A9_vs_WT_LFQ.txt"), 
                               header = T, sep = "\t", strip.white = TRUE))
r9b <- data.frame(read.delim(paste0(path, "221124_B9_vs_WT_LFQ.txt"), 
                               header = T, sep = "\t", strip.white = TRUE))
r9c <- data.frame(read.delim(paste0(path, "221124_C9_vs_WT_LFQ.txt"), 
                               header = T, sep = "\t", strip.white = TRUE))

expdes <- data.frame(read.delim(paste0(path_out, "221124_index.txt"), 
                               header = T, sep = "\t", strip.white = TRUE))
                               
