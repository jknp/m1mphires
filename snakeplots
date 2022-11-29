#JK 221128

##load packages, install if missing
require(readr)
require(dplyr)
require(DESeq2)
require(ggplot2)
require(ggpubr)
require(ggrepel)
require(PerformanceAnalytics)
require(DESeq2)
require(gridExtra)

#path setup
path = ""
path_plot = paste0(path, "")
path_counts = ""

#reading files
options(readr.num_columns = 0)

resA <- data.frame(read.delim(paste0(path, "221128_wt_vsA_snake.txt"), 
                               header = T, sep = "\t", strip.white = TRUE))
resB <- data.frame(read.delim(paste0(path, "221128_wt_vsB_snake.txt"), 
                              header = T, sep = "\t", strip.white = TRUE))
resC <- data.frame(read.delim(paste0(path, "221128_wt_vsC_snake.txt"), 
                              header = T, sep = "\t", strip.white = TRUE))

#For correlation matrix plot
index <- data.frame(read.delim(paste0(path, "221016_index.txt"), 
                                                 header = T, sep = "\t", strip.white = TRUE))
cts <- data.frame(read.delim(paste0(path_counts, "genecounts.txt"), 
                              header = T, sep = "\t", strip.white = TRUE))



##helper function for rankplot
rankplot = function(resinput, extr){
  require("ggplot2")
  require("ggrepel")
  #Takes DESeq2 results that were ordered and ranked (resinput),
  #and how many genes to label in both extremes
  
  #removes NA adjusted pvalue rows, calculates score and makes a rankplot
  sel <- resinput
  #select top 5 extremes for plot, add layer for coloring
  extr_top = c(head(sel$symbol, extr), tail(sel$symbol, extr))
  sel$extr = "black"
  sel$extr[(sel$symbol %in% extr_top)] = "royalblue2"
  
  rplot = ggplot(sel, aes(rank, log2FoldChange, label = symbol, col = extr)) +
    geom_point(size = 1, col = sel$extr) + 
    geom_text_repel(aes(rank, log2FoldChange), label = ifelse(sel$symbol %in% extr_top, sel$symbol,""),
                    colour = sel$extr, size = 3, max.overlaps = 100) + theme_bw() +
    ggtitle(paste0("Rankplot score with padj < ", 0.01))
  
  return(rplot)
}

#exploratory plot of all three DESeq2 results
extr = 5
p0 <- ggarrange(rankplot(resA,extr), rankplot(resB,extr), rankplot(resC,extr), ncol =3)

##Genelists based on plot
glist1 <- ""

#Extreme 5 of all three
glist1b <- ""

##Genelists based on intersections
### In all three comparisons up and down
glist2 <- ""

### In all three comparisons up and down and in two comparisons up and down
glist3 <- ""


##helper function for snakeplot
snakeplot = function(resinput, glist){
  require("ggplot2")
  require("ggrepel")
  #Takes DESeq2 results that were ordered and ranked (resinput),
  #and how many genes to label in both extremes
  
  sel <- resinput
  #select genelist
  sel$extr = "black"
  sel$extr[(sel$symbol %in% glist)] = "royalblue2"
  
  rplot = ggplot(sel, aes(rank, log2FoldChange, label = symbol, col = extr)) +
    geom_point(size = 1, col = sel$extr) + 
    geom_text_repel(aes(rank, log2FoldChange), label = ifelse(sel$symbol %in% glist, sel$symbol,""),
                    colour = sel$extr, size = 3, max.overlaps = 100) + theme_bw() +
    ggtitle(paste0("Rankplot score with padj < ", 0.01))
  
  return(rplot)
}


p0b <- ggarrange(snakeplot(resA,glist1b), snakeplot(resB,glist1b), snakeplot(resC,glist1b), ncol =3)

#plot with using genelist 2
p1 <- ggarrange(snakeplot(resA,glist2), snakeplot(resB,glist2), snakeplot(resC,glist2), ncol =3)

#plot with using genelist 3
p2 <- ggarrange(snakeplot(resA,glist3), snakeplot(resB,glist3), snakeplot(resC,glist3), ncol =3)

##Correlation plot of these snakeplots


scorplot <- function(index, glist_sel, cts){
  #selection from countdata
  txsel <- cts[which(cts$external_gene_id %in% glist_sel),]
  rownames(txsel) <- txsel$external_gene_id
  
  #Sample index data, intialize DESeq2 object
  index2 <- index
  index2[,1] <- paste0("X", index[,1])
  rownames(coldata) <- index2[,1]
  dds <- DESeqDataSetFromMatrix(countData = txsel[,5:13], 
                                colData = coldata[which(coldata$exposure == 9),],
                                design = ~sample)
  #Normalization
  dds2 <- estimateSizeFactors(dds)
  norm_cts <- counts(dds2, normalized=TRUE) 
  colnames(norm_cts) <- index$sample[4:12]
  
  #Dataframe with untransformed means per replicate for each gene
  avg = data.frame(norm_cts[,1:3])
  colnames(avg) = c("A9","B9","C9")
  avg$A9 <- rowMeans(norm_cts[, 1:3])
  avg$B9 <- rowMeans(norm_cts[, 5:6])
  avg$C9 <- rowMeans(norm_cts[, 7:9])
  avg <- log(avg[,1:3])
  
  #plotting chart
  return(chart.Correlation(avg))
}

##selection on the basis of plot1
glist_sel <- c("CDH12", "GFRA1", "MEF2C","ADGRL2","LINC01446","KLK3")

scorplot(index, glist_sel, cts) #Simple correlation plot based on hits in p1
scorplot(index, glist2, cts) #Using glist2
#scorplot(index, glist3, cts) #Using glist3
