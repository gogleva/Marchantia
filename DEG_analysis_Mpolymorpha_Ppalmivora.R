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
#put the variable of interest at the end of the formula, the control level is the first level.
sample_table$Experiment.type <- factor(sample_table$Experiment.type, levels = c('mock', 'infected'))

names(sample_table)[1] <- "Experiment" 
#rename for consistency:
coldata <- sample_table

#----2.DEG analysis with DESeq2
# Pair-wise comparissons, all the steps in a function for convenience and modularity

get_DEG <- function(pw_counts = fc_matrix,
                    pw_coldata = coldata,
                    day,
                    pv = 0.001,
                    lfc = 1.5,
                    res = FALSE){
    # produce DEG analysis for a given time point, output a df with significant results (padj and log2FC)
    
    # quick check: day should be numeric
    if (!(is.numeric(day))) stop('day must be numeric!')
    
    # subset relevant counts and metadata rows
    pw_coldata <- pw_coldata[pw_coldata$Time == paste0(day, 'd'),]
    pw_counts <- select(pw_counts,
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
        
        pw_resSig$day = day
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
                       up1 = make_set(DEG_1d, direction = 'up'),
                       up2 = make_set(DEG_2d, direction = 'up'),
                       up3 = make_set(DEG_3d, direction = 'up'),
                       up4 = make_set(DEG_4d, direction = 'up'))

upset(upset_me, sets = c('up1', 'up2', 'up3', 'up4'),
      order.by = "freq",
      sets.bar.color = 'black',
      main.bar.color = '#1B9E77',
      sets.x.label = "DEG per time point", 
      text.scale = c(1.5, 1.5, 1.1, 1.2, 1.5))

dev.off()

# down genes

upset_me_down <- data.frame(id = all_DEG_ids,
                            down1 = make_set(DEG_1d,
                                              direction = 'down'),
                            down2 = make_set(DEG_2d,
                                              direction = 'down'),
                            down3 = make_set(DEG_3d,
                                              direction = 'down'),
                            down4 = make_set(DEG_4d,
                                              direction = 'down'))
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

mp_vsd <- as.data.frame(assay(vsd)) %>%
          mutate(id = rownames(vsd))

mp_DEG <- filter(mp_vsd, id %in% all_DEG_ids) %>%
    select(id, everything())

rownames(mp_DEG) <- mp_DEG$id
mp_DEG <- select(mp_DEG, -id)
DEGmat <- mp_DEG - rowMeans(mp_DEG)
anno <- as.data.frame(colData(vsd)[, c('Time', 'Experiment')])

# heatmap of all DEGs!
pheatmap(DEGmat, cluster_cols = FALSE,
         show_rownames = FALSE,
         color = colorRampPalette(c("navy", "white", "#1B9E77"))(50),
         annotation = anno)


# Alluvial plots for functional summaries ----
## alluvial plots for Marchantia

library(alluvial)

# combine all DEG information from pair-wise comparisons in a tidy table:

DEG_summary <- rbind(DEG_1d, DEG_2d, DEG_3d, DEG_4d)
DEG_tidy <- DEG_summary %>%
            mutate(gene_id = rownames(DEG_summary)) %>%
            select(c(gene_id, log2FoldChange, padj, day))

# attach curated annotations
mp_annotation <- read_csv("data/Mpolymorpha_tidy_annotation.csv")
               
# up-regulated genes:
DEG_annotated_UP <- left_join(DEG_tidy, mp_annotation) %>%
                 filter(log2FoldChange > 0) %>%
                 group_by(day, type) %>%
                 summarise(count = n()) %>%
                 spread(key = day, value = count) %>%
                 replace(., is.na(.), 0) %>%
                 filter(type != 0) 

funs_up <- DEG_annotated_UP %>%
    gather(`1`, `2`, `3`, `4`,
           key = 'time_point', value = 'num_genes')
funs_up$time_point <- as.numeric(funs_up$time_point)

# alluvial plot:

set.seed(39) # for nice colours

cols <- hsv(h = sample(1:11/11), s = sample(1:11)/15, v = sample(3:11)/15)

alluvial_ts(funs_up,
            wave = .3,
            ygap = 5,
            col = cols,
            plotdir = 'centred',
            alpha= .9,
            grid = TRUE,
            grid.lwd = 5,
            xmargin = 0.4,
            lab.cex = .8,
            xlab = '',
            ylab = '',
            border = NA,
            axis.cex = .8, 
            leg.cex = .7,
            leg.col='white', 
            title = "Mpolymorpha, DEG functional summary")

# down-regulated genes

# up-regulated genes:
DEG_annotated_DOWN <- left_join(DEG_tidy, mp_annotation) %>%
    filter(log2FoldChange < 0) %>%
    group_by(day, type) %>%
    summarise(count = n()) %>%
    spread(key = day, value = count) %>%
    replace(., is.na(.), 0) %>%
    filter(type != 0) 

funs_down <- DEG_annotated_DOWN %>%
    gather(`1`, `2`, `3`, `4`,
           key = 'time_point', value = 'num_genes')
funs_down$time_point <- as.numeric(funs_down$time_point)

# alluvial plot:

alluvial_ts(funs_down,
            wave = .3,
            ygap = 5,
            col = cols,
            plotdir = 'centred',
            alpha= .9,
            grid = TRUE,
            grid.lwd = 5,
            xmargin = 0.4,
            lab.cex = .8,
            xlab = '',
            ylab = '',
            border = NA,
            axis.cex = .8, 
            leg.cex = .7,
            leg.col='white', 
            title = "Mpolymorpha, DEG functional summary")











