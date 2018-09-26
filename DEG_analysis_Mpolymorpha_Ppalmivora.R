### Marchantia infected with P.palmivora (ARItd) transcriptome analysis

#-------Set up the environment ----
source("http://bioconductor.org/biocLite.R")
biocLite()
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(genefilter)
library(gplots)
library(pheatmap)
library(stringr)
library(UpSetR)

# setwd('~/Marchantia_repo')

#----1.Preliminaries ----
#----download feature counts table ----
fc_table <- read.table("data/FeartureCounts_STAR_Mpolymorpha.tsv")


#----fix column names ----
names(fc_table) <- c("gene_id", 'chr', 'start', 'stop',
                     'strand', 'length',
                     'A1A', 'A1B', 'A1C', 
                     'A2A', 'A2B', 'A2C',
                     'A3A', 'A3B', 'A3C',
                     'A4A', 'A4B', 'A4C',
                     'M1A', 'M1B', 'M1C', 
                     'M2A', 'M2B', 'M2C',
                     'M3A', 'M3B', 'M3C', 
                     'M4A', 'M4B', 'M4C')

#----format count matrix ----
fc_matrix <- fc_table[,c(1,7:30)]
row.names(fc_matrix) <- fc_matrix[,1]
fc_matrix <- fc_matrix[,-c(1)]

#---- prep metadata ----
sample_table <- read.csv("data/sample_table_Mpolymorpha.csv",
                         header = T,
                         row.names = 1)
#to make it shorter:
names(sample_table)[1] <- "Experiment" 
#rename for consistency:
coldata <- sample_table

#put the variable of interest at the end of the formula, the control level is the first level.
coldata$Experiment <- factor(coldata$Experiment, levels = c('mock', 'infected'))

#----2.DEG analysis with DESeq2
# Pair-wise comparissons, all the steps in a function for convenience and modularity

get_DEG <- function(pw_counts = fc_matrix,
                    pw_coldata = coldata,
                    day,
                    pv = 0.001,
                    lfc = 2,
                    res = FALSE){
    # produce DEG analysis for a given time point, output a df with significant results (padj and log2FC)
    
    # quick check: day should be numeric
    if (!(is.numeric(day))) stop('day must be numeric!')
    
    # subset relevant counts and metadata rows
    pw_coldata <- coldata[coldata$Time == paste0(day, 'd'),]
    pw_counts <- select(fc_matrix,
                        contains(as.character(day)))
    
    # DEG analysis:
    pw_dds <- DESeqDataSetFromMatrix(countData = pw_counts,
                                     colData = pw_coldata,
                                     design = ~ Experiment)
    pw_dds <- pw_dds[rowSums(counts(pw_dds)) > 10, ]
    pw_dds <- DESeq(pw_dds)
    pw_res <- results(pw_dds, alpha = 0.001)
    
    if (res == TRUE) {
        return(pw_res) # convenient for volcanoplots
    } else {
        
        pw_res_ordered <- pw_res[order(pw_res$padj),]
        
        # Export significant results: |LFC| >= 2 and adjusted p-value <= 10^-3
        pw_resSig <- as.data.frame(subset(pw_res_ordered, padj < pv & abs(log2FoldChange) >= lfc))
        
        return(pw_resSig)       
    }
}

# helper function to construct volcanoplot from dds object
plot_volcano <- function(res_object,
                         lfc = 2,
                         pval = 0.001){
    tab <- data.frame(logFC = res_object$log2FoldChange,
                      negLogPval = -log10(res_object$padj))
    signGenes <- (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
    
    # plot
    par(mar = c(2, 2, 1, 1))
    plot(tab, pch = 16, cex = 0.6,
         xlim = c(-10,10),
         xlab = '',
         ylab = '',
         cex.axis = 1.0)
    points(tab[signGenes, ], pch = 19, cex = 0.8, col = "#E41A1C") 
    abline(h = -log10(pval), col = "#4DAF4A", lty = 3, lwd = 3) 
    abline(v = c(-lfc, lfc), col = "#377EB8", lty = 3, lwd = 3) 
}

# Now we can run pair-wise DEG tests for each time point independently (comparing infected vs mock):

#---- DEG results for pair-wise comparisons ----
DEG_1d <- get_DEG(day = 1)
DEG_2d <- get_DEG(day = 2)
DEG_3d <- get_DEG(day = 3)
DEG_4d <- get_DEG(day = 4)

#----
par(mfrow = c(4,1))
plot_volcano(get_DEG(day = 1, res = TRUE))
plot_volcano(get_DEG(day = 2, res = TRUE))
plot_volcano(get_DEG(day = 3, res = TRUE))
plot_volcano(get_DEG(day = 4, res = TRUE))

dev.off()

#### upset plots:

all_DEG_ids <- unique(c(rownames(DEG_1d), rownames(DEG_2d),
                        rownames(DEG_3d), rownames(DEG_4d)))

# 1916 up genes

up_1dpi <- as.numeric(all_DEG_ids %in% rownames(DEG_1d[DEG_1d$log2FoldChange >= 2,]))
up_2dpi <- as.numeric(all_DEG_ids %in% rownames(DEG_2d[DEG_2d$log2FoldChange >= 2,]))
up_3dpi <- as.numeric(all_DEG_ids %in% rownames(DEG_3d[DEG_3d$log2FoldChange >= 2,]))
up_4dpi <- as.numeric(all_DEG_ids %in% rownames(DEG_4d[DEG_4d$log2FoldChange >= 2,]))

upset_me <- data.frame(id = all_DEG_ids,
                       up1 = up_1dpi,
                       up2 = up_2dpi,
                       up3 = up_3dpi,
                       up4 = up_4dpi)

upset(upset_me, sets = c('up1', 'up2', 'up3', 'up4'),
      order.by = "freq",
      sets.bar.color = 'black',
      main.bar.color = '#1B9E77',
      sets.x.label = "DEG per time point", 
      text.scale = c(1.5, 1.5, 1.1, 1.2, 1.5))


dev.off()

# down genes

down_1dpi <- as.numeric(all_DEG_ids %in% rownames(DEG_1d[DEG_1d$log2FoldChange <= -2,]))
down_2dpi <- as.numeric(all_DEG_ids %in% rownames(DEG_2d[DEG_2d$log2FoldChange <= -2,]))
down_3dpi <- as.numeric(all_DEG_ids %in% rownames(DEG_3d[DEG_3d$log2FoldChange <= -2,]))
down_4dpi <- as.numeric(all_DEG_ids %in% rownames(DEG_4d[DEG_4d$log2FoldChange <= -2,]))

upset_me_down <- data.frame(id = all_DEG_ids,
                            down1 = down_1dpi,
                            down2 = down_2dpi,
                            down3 = down_3dpi,
                            down4 = down_4dpi)


upset(upset_me_down,
      sets = c('down4', 'down3', 'down2', 'down1'),
      order.by = "freq",
      keep.order = TRUE,
      sets.bar.color = 'black',
      main.bar.color = '#7570B3',
      sets.x.label = "DEG per time point", 
      text.scale = c(1.5, 1.3, 1, 1, 1.2, 1.2))

dev.off()

#### heatmaps of the most significant DEG overview

library(genefilter)
dds <- DESeqDataSetFromMatrix(countData = fc_matrix,
                              colData = coldata,
                              design = ~ Experiment)

vsd <- vst(dds, blind = FALSE) # variance-stabilised counts

topVarGenes <- head(order(rowVars(assay(vsd)), 
                          decreasing = TRUE),
                    1000)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c('Time', 'Experiment')])

pheatmap(mat, cluster_cols = FALSE,
         show_rownames = FALSE,
         color = colorRampPalette(c("navy", "white", "#1B9E77"))(50))

dev.off()














