
#Make a directory, install all required packages 
#(clusterProfiler and enrichplot might need separate dependencies)
setwd("~/RNAseq/2020_P9_LNMph_v3/")

require(readr)
require(dplyr)
require(DESeq2)
require(ggplot2)
require(ggpubr)
require(tidyverse)
require(NMF)
require(EnhancedVolcano)
require(msigdbr)
require(org.Hs.eg.db)
require(magrittr)
require(readr)
require(clusterProfiler)
require(DOSE)
require(enrichplot)

#mkdir folders as follows, setting variables
path = "~/RNAseq/2020_P9_LNMph_v3/"
path_gcf = paste0(path,"GCF/")
path_plot = paste0(path, "R/Rplot/")
path_gsea = paste0(path,"GSEA/")
path_out = paste0(path, "GSEA/out/")
path_xena = paste0(path, "Xena/")

#loading raw data
options(readr.num_columns = 0)
index <- data.frame(read.delim(paste0(path, "R/221016_index.txt"), 
                               header = T, sep = "\t", strip.white = TRUE))
norm_10m <- data.frame(read.delim(paste0(path, "GCF/genecounts_normalized_10mil.txt"), 
                                  header = T, sep = "\t", strip.white = TRUE))
counts <- data.frame(read.delim(paste0(path, "GCF/genecounts.txt"), 
                                header = T, sep = "\t", strip.white = TRUE))

#Function definitions
#DESEq2 normalization function
normalizeWithDESeq2 = function(mat){ #Credits: Gergana Bounova    
  require(DESeq2) || stop("Need package \"DESeq2\" to perform this normalization.")
  ys = estimateSizeFactorsForMatrix(mat)
  yy = mat
  
  # divide by size factor, explained on page 4 in DESeq manual
  # for (i in 1:ncol(yy)){ yy[,i]=yy[,i]/ys[i] }
  yy = sweep(yy,2,ys,'/')
  
  return(yy)
}

#Genomic Core Facility (GCF) readcount parser for quick formatting
parse_readcounts_gcf = function(counts,p=c("no","yes")){
  p <- match.arg(p)
  if (p == "yes"){
    png("normalizationplots.png",res=300,units='in',width=14,height=7)
  }
  #Separate counts and gene annotation data:
  countdata<-counts[,2:(ncol(counts))]
  ncols <- ncol(countdata)
  annot<-counts[,1]  #add more columns if annotation is there
  sums<-apply(countdata,2,sum) #calculate number of reads per sample
  sums_zs<-(sums-mean(sums))/sd(sums) #calculate the average readcount per sample, and for each sample how many standard deviations they differ from the average
  print("Consider to remove samples with low readcounts:")
  if (length(sums[which(sums_zs < -1.5)]) > 0){
    print(as.data.frame(sums[which(sums_zs < -1.5)]))
  } else {
    print("None.")
  }
  cols<-rep("black",ncols)
  cols[which(sums_zs < -1.5)]<-"orange"
  cols[which(sums_zs < -2)]<-"red"
  par(mfcol=c(1,2))
  #Plot readcount distribution before normalization:
  countdata_exprs<-countdata[which(apply(countdata,1,max) >= 10),]
  if (ncols < 20){
    boxplot(log(countdata_exprs),las=2,ylab="log(readcounts per gene)",main="raw data",border=cols)
  } else {
    plot(density(log(countdata_exprs[,1])),main="raw data",xlab="log(readcounts per gene)",ylim=c(0,0.25),col=cols[1])
    for (s in 2:ncols){
      lines(density(log(countdata_exprs[,s])),col=cols[s])
    }
  }
  legend('topleft',lwd=1,col=c('red','orange'),legend=c('total reads > 2 st dev below mean','total reads > 1.5 st dev below mean'),cex=0.5)
  countdata_norm<-normalizeWithDESeq2(as.matrix(countdata)) #Normalize using DESeq2
  countdata_norm_exprs<-countdata_norm[which(apply(countdata_norm,1,max) >= 10),] #Remove non-expressed genes
  countdata_norm_log<-log(countdata_norm_exprs+1) #Log transform
  #Plot readcount distribution normalized & transformed data:
  if (ncols < 20){
    boxplot(countdata_norm_log,las=2,ylab="log(readcounts per gene)",main="normalized data",border=cols)
  } else {
    plot(density(countdata_norm_log[,1]),main="normalized data",xlab="log(readcounts per gene)",ylim=c(0,0.25),col=cols[1])
    for (s in 2:ncols){
      lines(density(countdata_norm_log[,s]),col=cols[s])
    }
  }
  #Re-attach gene information to transformed counts:
  countdata_transformed<-cbind(annot[which(apply(countdata_norm,1,max) >= 10)],countdata_norm_log) #deleted , in >= 10),]
  rownames(countdata_transformed)<-countdata_transformed[,1]
  par(mfrow=c(1,1))
  if (p == "yes"){
    dev.off()
  }
  return(countdata_transformed)
}

#Quick functions for rankplots and genelist exporting
rankplot = function(resinput,cfpadj, cf_l2fc,top){
  require("ggplot2")
  #Takes DESeq2 results (resinput), padj cut-off (cfpadj), log2fc cutoff (cf_fc) and how many genes to label
  
  #removes NA adjusted pvalue rows, calculates score and makes a rankplot
  resinput <- data.frame(resinput[!is.na(resinput$padj),])
  resinput$mlogpadj = -log(resinput$padj)
  resinput$score = as.numeric(resinput$mlogpadj * resinput$log2FoldChange)
  resinput = resinput[order(resinput$score),]
  
  #remove insignificant data
  sel = resinput[which(resinput$padj < cfpadj),]
  sel = sel[which(sel$log2FoldChange >= cf_l2fc | sel$log2FoldChange <= -cf_l2fc),]
  sel$rank = 1:nrow(sel)
  
  #select top 5 extremes for plot, add layer for coloring
  extr_top = c(head(sel$symbol, top), tail(sel$symbol, top))
  sel$extr = "black"
  sel$extr[(sel$symbol %in% extr_top)] = "royalblue2"
  
  rplot = ggplot(sel, aes(rank, score, label = symbol, col = extr)) + 
    geom_point(alpha = 0.5, size = 1, col = sel$extr) + 
    geom_text_repel(aes(rank, score), label = ifelse(sel$symbol %in% extr_top, sel$symbol,""),
                    colour = sel$extr, size = 2.5, max.overlaps = 100) + theme_bw() +
    ggtitle(paste0("Rankplot score with cut-offs padj < ",cfpadj, " and l2fc > ", cf_l2fc, " or < ", -cf_l2fc))
  
  return(rplot)
  ##add colors and text repel
}

export_glist = function(resinput,cfpadj, cf_l2fc){
  #Takes DESeq2 results an exports after cutoffs are applied like in rankplots
  resinput <- data.frame(resinput[!is.na(resinput$padj),])
  resinput$mlogpadj = -log(resinput$padj)
  resinput$score = as.numeric(resinput$mlogpadj * resinput$log2FoldChange)
  resinput = resinput[order(resinput$score),]
  
  #cutoffs
  sel = resinput[which(resinput$padj < cfpadj),]
  sel = sel[which(sel$log2FoldChange >= cf_l2fc | sel$log2FoldChange <= -cf_l2fc),]
  
  return(sel)
}