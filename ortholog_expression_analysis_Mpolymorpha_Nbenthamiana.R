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
library(ggpmisc) # for nice square numbers

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
mp_sample_table$Experiment.type <- factor(mp_sample_table$Experiment.type, levels = c("mock", "infected"))
names(mp_sample_table)[1] <- "Experiment" 
mp_coldata <- mp_sample_table

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
                   M72A:M72C) %>%
               as.data.frame()

#These two steps will make sure that the rownames are now the rowids/benthi loci, then delete the redundant rowids columns

row.names(nb_fc_table) <- nb_fc_table$rowids
nb_fc_matrix <- select(nb_fc_table, -rowids)

#---- prep metadata
nb_sample_table <- read.csv("data/sample_table_Nbenthamiana.csv",
                         header = T,
                         row.names = 1)
nb_sample_table$Experiment.type <- factor(nb_sample_table$Experiment.type, levels = c("mock", "infected"))
names(nb_sample_table)[1] <- "Experiment" 
nb_coldata <- nb_sample_table

#---- 2. DEG analysis: pair-wise comparisons between infected-mock for each plant and time point.
# we want to keep LFC values in 4 stages of infection, filter by adjsted p-value. Will use LFC later to compare expression patterns of orthologues genes.

### DEGs for orthologs:
# (wee need LFC in 4 time points, with filtering by p-value, w/o filtering by LFC)

get_DEG_LFC <- function(pw_counts,
                        pw_coldata,
                        species = c('Mpolymorpha', 'Nbenthamina'),
                        day_time, 
                        pv = 0.001,
                        set_stage
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
    
    pw_resSig <- pw_resSig %>%
                 mutate(day_time = set_stage) %>%
                 select(gene_id, day_time, everything())
    return(pw_resSig)       
}


# MP LFC ----

mp_deg1 <- get_DEG_LFC(pw_counts = mp_fc_matrix,
                       pw_coldata = mp_sample_table,
                       species = 'Mpolymorpha',
                       day_time = 1,
                       set_stage = 1)
mp_deg2 <- get_DEG_LFC(pw_counts = mp_fc_matrix,
                       pw_coldata = mp_sample_table,
                       species = 'Mpolymorpha',
                       day_time = 2,
                       set_stage = 2)

mp_deg3 <- get_DEG_LFC(pw_counts = mp_fc_matrix,
                       pw_coldata = mp_sample_table,
                       species = 'Mpolymorpha',
                       day_time = 3,
                       set_stage = 3)
mp_deg4 <- get_DEG_LFC(pw_counts = mp_fc_matrix,
                       pw_coldata = mp_sample_table,
                       species = 'Mpolymorpha',
                       day_time = 4,
                       set_stage = 4)
# combine in one object

mp_deg_all <- rbind(mp_deg1, mp_deg2,
                    mp_deg3, mp_deg4) %>%
              mutate(species = 'Mpoly')


# Niben LFC: to do - optimize this ----

nb_deg14 <- get_DEG_LFC(pw_counts = nb_fc_matrix,
                    pw_coldata = nb_sample_table,
                    species = 'Nbenthamiana',
                    day_time = 14,
                    set_stage = 1)

nb_deg24 <- get_DEG_LFC(pw_counts = nb_fc_matrix,
                        pw_coldata = nb_sample_table,
                        species = 'Nbenthamiana',
                        day_time = 24,
                        set_stage = 2)

nb_deg48 <- get_DEG_LFC(pw_counts = nb_fc_matrix,
                        pw_coldata = nb_sample_table,
                        species = 'Nbenthamiana',
                        day_time = 48,
                        set_stage = 3)

nb_deg72 <- get_DEG_LFC(pw_counts = nb_fc_matrix,
                        pw_coldata = nb_sample_table,
                        species = 'Nbenthamiana',
                        day_time = 72,
                        set_stage = 4)

nb_deg_all <- rbind(nb_deg14, nb_deg24,
                    nb_deg48, nb_deg72) %>%
              mutate(species = 'Niben')


# bind lists of potentialy deg genes (nb: here we filter just by adjusted p-value, not by LFC --> gene lists obtained with DEG_analysis_Mpolymorpha_Ppalmivora.R and DEG_analysis_Nbenthamiana_Ppalmivora.R scripts will be slightly different)

deg_nb_mp <- rbind(mp_deg_all,
                   nb_deg_all)


#----Parsing OrthoFinder output----

single_copy_OG <- read_tsv('data/SingleCopyOrthogroups.txt',
                           col_names = FALSE)

names(single_copy_OG) <- 'OG'

orthogroups <- read_delim('data/Orthogroups.csv',
                          delim = '\t')
names(orthogroups) <- c('OG', 'Mpoly', 'Niben')

# extract single copy OG and gene lists:

scp_og <- left_join(single_copy_OG, orthogroups) %>%
    separate(Mpoly, '\\.', into = c('Mpoly', 
                                    'the_rest')) %>%
    separate(Niben, '\\.', into = c('Niben', 'stuff')) %>% select(-c(the_rest, stuff)) %>%
    gather(Mpoly, Niben, key = 'species', value = 'gene_id')

# attach DEG tables with OG id:
ogdeg <- left_join(scp_og, deg_nb_mp)

# prepare for plotting (there might be a better way)

plotdat_mp <- select(ogdeg, c(species, OG, gene_id,
                              day_time, log2FoldChange)) %>%
    filter(species == 'Mpoly') %>%
    select(-species) %>%
    rename(log2FoldChange = 'mp_lfc',
           gene_id = 'MP_gene_id')

plotdat_nb <- select(ogdeg, c(species, OG, gene_id,
                              day_time, log2FoldChange)) %>%
    filter(species == 'Niben') %>%
    select(-species) %>%
    rename(log2FoldChange = 'nb_lfc',
           gene_id = 'NB_gene_id') 

plotdat_mp_nb <- left_join(plotdat_mp, plotdat_nb) %>%
    filter(!is.na(mp_lfc) & !is.na(nb_lfc))

# rough plot without annotations
ggplot(plotdat_mp_nb, aes(x = nb_lfc, y = mp_lfc)) +
    geom_point(aes(alpha = 0.5)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)  +
    facet_wrap(~ day_time)

# attach Marchantia curated annotations
mp_annotation <- read_csv("data/Mpolymorpha_tidy_annotation.csv")
names(mp_annotation)[1] <- 'MP_gene_id'

plotdat_mp_nb_annotated <- left_join(plotdat_mp_nb, mp_annotation)


# full picture
ggplot(plotdat_mp_nb_annotated, aes(x = nb_lfc, y = mp_lfc)) +
    geom_point(aes(color = type)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)  +
    facet_wrap(~day_time)

# remove genes w/o selected annotation
plotdat_mp_nb_annotated %>%
    filter(!is.na(type)) %>%
    ggplot(aes(x = nb_lfc, y = mp_lfc)) +
    geom_point(aes(color = type)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)  +
    facet_wrap(~day_time)

### scatterplots for the paper

plotdat_mp_nb_annotated %>%
    arrange(!is.na(type), type) %>%
    ggplot(aes(x = nb_lfc, y = mp_lfc)) +
    geom_point(aes(color = type)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0)  +
    stat_quadrant_counts(colour = 'black') +
    facet_wrap(~day_time) +
    xlab('LFC N. benthamiana') +
    ylab('LFC M. polymorpha') +
    theme_bw()


## Plot heatmaps for single-copy orthologs, belonging to the falvonoid pathway

# reference table with inforamtion about gene ids, OG and annotation of single-copy genes in Mp and Nb involved in flavonoid pathway (9 genes):

flavonoids_sc <- plotdat_mp_nb_annotated %>%
                 filter(type == 'flavonoid pathway') %>%
                 select(c(OG, MP_gene_id, NB_gene_id, type,
                          description)) %>%
                 distinct()


## helper functions to construct heatmaps for specified gene lists

# helper function to  construct dds object, perform count transformation and extract values for relevant genes to build hetamaps

get_mat <- function(gene_list,
                    fc_matrix,
                    coldata){

    dds <- DESeqDataSetFromMatrix(countData = fc_matrix,
                                     colData = coldata,
                                     design = ~ Experiment)
    # regularised log transformed conts
    rld <- rlog(dds, blind = FALSE)
    rld <- as.data.frame(assay(rld))
    # center values:
    rld_centered <- rld - rowMeans(rld)
    
    # extract only genes of interest:
    
    plot_mat <- rld_centered[gene_list,]
    return(plot_mat)
}

# now, get the gene matrices for Mpolymorpha and Nbenthamiana
mp_mat <- get_mat(gene_list = flavonoids_sc$MP_gene_id,
                  fc_matrix = mp_fc_matrix,
                  coldata = mp_coldata)

nb_mat <- get_mat(gene_list = flavonoids_sc$NB_gene_id,
                  fc_matrix = nb_fc_matrix,
                  coldata = nb_coldata)


# helper function to construct table with functional annotations and OG id

get_annotations <- function(annotation, 
                            plot_mat,
                            id_column
                            ){
    df1 <- as.tibble(plot_mat) %>%
        mutate(id = rownames(plot_mat)) %>%
        select(c(id, everything()))
    names(df1)[1] <- id_column
    df1 <- left_join(df1, flavonoids_sc) %>%
           arrange(description, OG)
    return(df1)
}


mp_an <- get_annotations(annotation = flavonoids_sc,
                plot_mat = mp_mat,
                id_column = 'MP_gene_id')

nb_an <- get_annotations(annotation = flavonoids_sc,
                      plot_mat = nb_mat,
                      id_column = 'NB_gene_id')






