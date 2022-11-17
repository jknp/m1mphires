##############
##Make a directory, install all required packages 
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

#########################
##Function definitions
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

##########
##Read quality control plots
#loading bamstats, presorted to match index key
bamstats <- data.frame(read.delim(paste0(path_gcf, "bamstats.txt"),
                               header = T, sep = "\t", strip.white = TRUE),stringsAsFactors = T)
colnames(bamstats) <- c("sample", "readno", "unmapreadno", "fmapped", "nonprim_alignmentsno", "pair", "properpair", "singleread", "fproper", "dups", "fdups", "failed","mapq5no","fmapq5")
bamstats <- bamstats[,-5]

bamstats$key <- index$sample
bamstats$biorep <- index$biorep

#loading gene_exp_stats
expstats <- data.frame(read.delim(paste0(path_gcf, "gene_expstats.txt"),
                               header = T, sep = "\t", strip.white = TRUE),stringsAsFactors = T)
expstats$key <- index$sample
expstats$biorep <- index$biorep

#loading strandstats
strstats <- data.frame(read.delim(paste0(path_gcf, "strandedstats.txt"),
                               header = T, sep = "\t", strip.white = TRUE),stringsAsFactors = T)
strstats$key <- index$sample
strstats$biorep <- index$biorep

#plots
q1 = ggplot(bamstats, aes(singleread, key, fill = biorep)) + geom_bar(stat = "identity") + theme_bw() +
  ylab("") + xlab("Single mapped reads")

q2 = ggplot(bamstats, aes(fmapq5, key, fill = biorep)) + geom_bar(stat = "identity") + theme_bw() +
  ylab("") + xlab("Fraction mapq >= 5")

q3 = ggplot(expstats, aes(X.mtRNA, key, fill = biorep)) + geom_bar(stat = "identity") + theme_bw() +
  ylab("") + xlab("% mitochondrial RNA")

q4 = ggplot(strstats, aes(reads_on_same_strand, key, fill = biorep)) + geom_bar(stat = "identity") + theme_bw() +
  ylab("") + xlab("Fraction on same strand")


#combining plots
p1 = ggarrange(q1,q3, q2, q4, nrow =2, ncol =2, common.legend = T, legend = "bottom")
annotate_figure(p1, top = text_grob("Read number, quality, %mito, reads on same strand", 
               color = "black", size = 14))
ggsave(paste0(path_plot,"221021_QCplots.png"))
ggsave(paste0(path_plot,"221021_QCplots.pdf"))

#################
##Dataframe pre-processing, adding gene annotations etc in one file

#Split the countdata in a table with just the (normalized) readcounts and the gene annotation for ENSG and HGNC:
counts_only <- norm_10m[,2:13]
annot <- data.frame(norm_10m$ensembl_gene_id)
annot$symbol <- mapIds(org.Hs.eg.db, keys=annot$norm_10m.ensembl_gene_id, keytype = "ENSEMBL", column="SYMBOL")
colnames(annot)[1] <- "gene"

annot$name <- paste0(annot$gene, "_", annot$symbol)
colnames(annot) <- c("gene", "symbol", "name")
#manual save in 221021_dendrogram.png


####################
##Clustered heatmap

#get the standard deviation of the expression per gene and order by sd
most_variable_counts <- data.frame(counts_only)
most_variable_counts$sd <- apply(counts_only,1,sd)

#Add annotation, but remove NO_GENE_NAME for these heatmap plots >> 41389 genes from 61806 total remain.
most_variable_counts$name <- annot$name
most_variable_counts <- most_variable_counts[with(most_variable_counts,order(-sd)),]
sel <- most_variable_counts[with(most_variable_counts,order(-sd)),]

#some data class manipulations for aheatmap
mvc = as.matrix(sel[,-14])
row.names(mvc) <- most_variable_counts$name
class(mvc) = "numeric"

#visualize the top50 most variable genes
topno = 50
aheatmap(
  mvc[1:topno, (1:ncol(mvc)-1)],
  labRow = sel$gene[1:topno],
  labCol = index[,"sample"],
  annCol = index[,c("biorep","techrep", "exposure")] #Do these genes correlate with any of the clinical variables?
)
#######################
