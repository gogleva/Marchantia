### Nicotiona benthamiana - Ppalmivora (ARItd) Leaf transcriptome analysis

#-------Set up the environment
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

#----1.Preliminaries ----
#----download feature counts table ----

# DELETEME
#setwd('/home/anna/anna/Labjournal/Manuscripts/marchantia_II_plant_response/Marchantia_repo/')

# Load the featurecounts

mycounts <- read_csv("data/FeatureCounts_STAR_Nbenthamiana.csv")

#remove columns that we dont need anymore
cleancounts <- select(mycounts, -c(X2:X6))

#Give the columns appropriate names
names(cleancounts) <- c('rowids', 
                        'A72A', 'A72B', 'A72C',
                        'A14A', 'A14B', 'A14C',
                        'A24A', 'A24B', 'A24C',
                        'A48A', 'A48B', 'A48C', 
                        'M72A', 'M72B', 'M72C',
                        'M14A', 'M14B', 'M14C',
                        'M24A', 'M24B', 'M24C',
                        'M48A', 'M48B', 'M48C')

#re-order columns (grouped from 14-72hpi in ARI vs MOCK treatments)

fc_table <- select(cleancounts,
                   rowids,
                   A14A:A48C,
                   A72A:A72C,
                   M14A:M48C,
                   M72A:M72C)

#These two steps will make sure that the rownames are now the rowids/benthi loci, then delete the redundant rowids columns

row.names(fc_table) <- fc_table$rowids
fc_table <- select(fc_table, -rowids)

#---- prep metadata ----
sample_table <- read.csv("data/sample_table_Nbenthamiana.csv",
                         header = T,
                         row.names = 1)

#to make it shorter:
names(sample_table)[1] <- "Experiment" 
#rename for consistency:
coldata <- sample_table

#put the variable of interest at the end of the formula, the control level is the first level.
coldata$Experiment <- factor(coldata$Experiment, levels = c('mock', 'infected'))

#----2. DEG analysis with DESeq2
# Pair-wise comparissons, all the steps in a function for convenience and modularity

#----4. DEG analysis: pair-wise comparissons, put all the steps in a helper function

get_DEG <- function(pw_counts = fc_table,
                    pw_coldata = coldata,
                    hour,
                    pv = 0.001,
                    lfc = 2,
                    res = FALSE){
    
    # produce DEG analysis for a given time point, output a df with significant results (padj and log2FC)
    
    # quick check: day should be numeric
    if (!(is.numeric(hour))) stop('hour must be numeric!')
    
    # subset relevant counts and metadata rows
    pw_coldata <- coldata[coldata$Time == paste0(hour, 'h'),]
    pw_counts <- select(fc_table,
                        contains(as.character(hour)))
    
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

# volcanoplot from dds object
plot_volcano <- function(res_object,
                         lfc = 2,
                         pval = 0.001){
    tab <- data.frame(logFC = res_object$log2FoldChange,
                      negLogPval = -log10(res_object$padj))
    signGenes <- (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
    
    # plot
    par(mar = c(2, 2, 1, 1))
    plot(tab, pch = 16, cex = 0.6,
         xlim = c(-15,15),
         xlab = '',
         ylab = '',
         cex.axis = 1.0)
    points(tab[signGenes, ], pch = 19, cex = 0.8, col = "#E41A1C") 
    abline(h = -log10(pval), col = "#4DAF4A", lty = 3, lwd = 3) 
    abline(v = c(-lfc, lfc), col = "#377EB8", lty = 3, lwd = 3) 
}

# Now we can run pair-wise DEG test for each time point independently

DEG_14h <- get_DEG(hour = 14)
DEG_24h <- get_DEG(hour = 24)
DEG_48h <- get_DEG(hour = 48)
DEG_72h <- get_DEG(hour = 72)

#----
par(mfrow = c(4,1))
# here we do DEG analysis and do volc. plots simultaneously:
plot_volcano(get_DEG(hour = 14, res = TRUE))
plot_volcano(get_DEG(hour = 24, res = TRUE))
plot_volcano(get_DEG(hour = 48, res = TRUE))
plot_volcano(get_DEG(hour = 72, res = TRUE))

dev.off()

#### upset plots:

all_DEG_ids <- unique(c(rownames(DEG_14h),
                        rownames(DEG_24h),
                        rownames(DEG_48h),
                        rownames(DEG_72h)))

# 1916 up genes
# helper function to extract up/down gene sets in upset-compatible format
make_set <- function(deg_set, 
                     lfc = 2, # default
                     direction){
    # deg_set - deg set we want to convert to upset-compatible set;
    # lfc - lfc thershold
    # direction - 'up' for up-regulated or 'down' for down-regulated genes
    
    if(direction == 'up'){
        res_set <- as.numeric(
            all_DEG_ids %in% rownames(deg_set[deg_set$log2FoldChange >= lfc,]))
    } else {
        res_set <- as.numeric(
            all_DEG_ids %in% rownames(deg_set[deg_set$log2FoldChange <= -lfc,]))                  
    } 
}

upset_me <- data.frame(id = all_DEG_ids,
                       up14 = make_set(DEG_14h, direction = 'up'),
                       up24 = make_set(DEG_24h, direction = 'up'),
                       up48 = make_set(DEG_48h, direction = 'up'),
                       up72 = make_set(DEG_72h, direction = 'up'))

upset(upset_me, sets = c('up14', 'up24', 'up48', 'up72'),
      order.by = "freq",
      sets.bar.color = 'black',
      main.bar.color = '#1B9E77',
      sets.x.label = "DEG per time point", 
      text.scale = c(1.5, 1.5, 1.1, 1.2, 1.5))

dev.off()

# down genes

upset_me_down <- data.frame(id = all_DEG_ids,
                            down14 = make_set(DEG_14h,
                                              direction = 'down'),
                            down24 = make_set(DEG_24h,
                                              direction = 'down'),
                            down48 = make_set(DEG_48h,
                                              direction = 'down'),
                            down72 = make_set(DEG_72h,
                                              direction = 'down'))
upset(upset_me_down,
      sets = c('down72', 'down48', 'down24', 'down14'),
      order.by = "freq",
      keep.order = TRUE,
      sets.bar.color = 'black',
      main.bar.color = '#7570B3',
      sets.x.label = "DEG per time point", 
      text.scale = c(1.5, 1.3, 1, 1, 1.2, 1.2))
