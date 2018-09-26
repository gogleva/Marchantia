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
cleancounts <- select(nb_mycounts, -c(X2:X6))

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
nb_fc_table <- select(nb_fc_table, -rowids)

#---- prep metadata
nb_sample_table <- read.csv("data/sample_table_Nbenthamiana.csv",
                         header = T,
                         row.names = 1)
names(nb_sample_table)[1] <- "Experiment" 
nb_coldata <- nb_sample_table

#put the variable of interest at the end of the formula, the control level is the first level.
nb_coldata$Experiment <- factor(nb_coldata$Experiment, levels = c('mock', 'infected'))



