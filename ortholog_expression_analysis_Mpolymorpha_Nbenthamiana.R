### Comapare expression patterns of ortholog genes between M.polymorpha and N.benthamiana during infection with P.palmivora (ARI-tdTomdato)

#---- Set up the environment ----

source("http://bioconductor.org/biocLite.R")
biocLite()
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(genefilter)
library(gplots)
library(pheatmap)
library(stringr)

#---- 1. Load raw counts and sample metadata for M.polymorpha - ARI-tdTomato and N.benthamiana - ARI-tdTomato timecourses. ----

# Marchantia polymorpha + P.palmivora (ARI-tdTomato) ----

#----download feature counts table
mp_fc_table <- read.table("data/FeartureCounts_STAR_Mpolymorpha.tsv")

#----fix column names
names(mp_fc_table) <- c("gene_id", 'chr', 'start', 'stop',
                     'strand', 'length',
                     'A1A', 'A1B', 'A1C', 
                     'A2A', 'A2B', 'A2C',
                     'A3A', 'A3B', 'A3C',
                     'A4A', 'A4B', 'A4C',
                     'M1A', 'M1B', 'M1C', 
                     'M2A', 'M2B', 'M2C',
                     'M3A', 'M3B', 'M3C', 
                     'M4A', 'M4B', 'M4C')

#----format count matrix
mp_fc_matrix <- mp_fc_table[,c(1,7:30)]
rownames(mp_fc_matrix) <- mp_fc_matrix[,1]
mp_fc_matrix <- mp_fc_matrix[,-c(1)]

#---- prep metadata
mp_sample_table <- read.csv("data/sample_table_Mpolymorpha.csv",
                         header = T,
                         row.names = 1)
names(mp_sample_table)[1] <- "Experiment" 
mp_coldata <- mp_sample_table
mp_coldata$Experiment <- factor(mp_coldata$Experiment, levels = c('mock', 'infected'))

# Nicotiana benthamiana + P.palmivora (ARI-tdTomato) ----

# load raw counts
nb_counts <- read_csv("data/FeatureCounts_STAR_Nbenthamiana.csv")

#remove columns that we dont need anymore
cleancounts <- select(nb_counts, -c(X2:X6))

# rename columns
names(cleancounts) <- c('rowids', 
                        'A72A', 'A72B', 'A72C',
                        'A14A', 'A14B', 'A14C',
                        'A24A', 'A24B', 'A24C',
                        'A48A', 'A48B', 'A48C', 
                        'M72A', 'M72B', 'M72C',
                        'M14A', 'M14B', 'M14C',
                        'M24A', 'M24B', 'M24C',
                        'M48A', 'M48B', 'M48C')

# re-order columns (grouped from 14-72hpi in ARI vs MOCK treatments)
nb_fc_table <- select(cleancounts,
                   rowids,
                   A14A:A48C,
                   A72A:A72C,
                   M14A:M48C,
                   M72A:M72C)

#These two steps will make sure that the rownames are now the rowids/benthi loci, then delete the redundant rowids columns

row.names(nb_fc_table) <- nb_fc_table$rowids
nb_fc_matrix <- select(nb_fc_table, -rowids)

#---- prep metadata
nb_sample_table <- read.csv("data/sample_table_Nbenthamiana.csv",
                         header = T,
                         row.names = 1)
names(nb_sample_table)[1] <- "Experiment" 
nb_coldata <- nb_sample_table

#put the variable of interest at the end of the formula, the control level is the first level.
nb_coldata$Experiment <- factor(nb_coldata$Experiment, levels = c('mock', 'infected'))

#---- 2. DEG analysis: pair-wise comparisons between infected-mock for each plant and time point.
# we want to keep LFC values in 4 stages of infection, filter by adjsted p-value. Will use LFC later to compare expression patterns of orthologues genes.

### DEGs for orthologs:
# (wee need LFC in 4 time points, with filtering by p-value, w/o filtering by LFC)

get_DEG_LFC <- function(pw_counts,
                        pw_coldata,
                        species = c('Mpolymorpha', 'Nbenthamina'),
                        day_time, 
                        pv = 0.001
                        ){
    # produce DEG analysis for a given time point, output a df with significant results (padj and log2FC)
    
    # quick check: day_time should be numeric
    if (!(is.numeric(day_time))) stop('day_time must be numeric!')
    
    # subset relevant counts and metadata rows
    if (species == 'Mpolymorpha'){
        pw_coldata <- pw_coldata[pw_coldata$Time == paste0(day_time, 'd'),]
    }
    
    if (species == 'Nbenthamiana'){
        pw_coldata <- pw_coldata[pw_coldata$Time == paste0(day_time, 'h'),]
    }
    
    pw_counts <- select(pw_counts,
                        contains(as.character(day_time)))
    
    # DEG analysis:
    pw_dds <- DESeqDataSetFromMatrix(countData = pw_counts,
                                     colData = pw_coldata,
                                     design = ~ Experiment)
    pw_dds <- pw_dds[rowSums(counts(pw_dds)) > 10, ]
    pw_dds <- DESeq(pw_dds)
    pw_res <- results(pw_dds, alpha = 0.001)
    
    pw_res_ordered <- pw_res[order(pw_res$padj),]
        
    # Filter genes by adjusted p-value <= 10^-3
    pw_resSig <- as.data.frame(subset(pw_res_ordered, padj < pv))
    pw_resSig$gene_id <- rownames(pw_resSig)
    return(pw_resSig)       
}

# test run
mp_1 <- get_DEG_LFC(pw_counts = mp_fc_matrix,
            pw_coldata = mp_sample_table,
            species = 'Mpolymorpha',
            day_time = 1)

nb_14 <- get_DEG_LFC(pw_counts = nb_fc_table,
                    pw_coldata = nb_sample_table,
                    species = 'Nbenthamiana',
                    day_time = 14)






