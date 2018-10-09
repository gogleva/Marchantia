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

# add heatmaps for single-copy orthologs

# marchantia:

mp_dds <- DESeqDataSetFromMatrix(countData = mp_fc_matrix,
                                 colData = mp_coldata,
                                 design = ~ Experiment)

mp_vsd <- vst(mp_dds, blind = FALSE) # variance-stabilised counts
mp_vsd <- as.data.frame(assay(mp_vsd))
mp_vsd_centered <- mp_vsd - rowMeans(mp_vsd)

# 1. select genes from single-copy ortholog list
# 2. attach functional annotations
# 3. filter and split values and annotations in 2 separate matrices for pheatmap






# select genes from single-copy ortholog list only
mp_DEG_sc <- filter(mp_vsd, id %in% plotdat_mp_nb_annotated$MP_gene_id) %>%
    select(id, everything())

rownames(mp_DEG_sc) <- mp_DEG_sc$id
mp_DEG_sc <- select(mp_DEG_sc, -id)
mp_DEGmat_sc <- mp_DEG_sc - rowMeans(mp_DEG_sc)

anno <- as.data.frame(colData(mp_dds)[, c('Time', 'Experiment')])

## now, we need to select genes belonging to functional categories we are interested in:

our_cats <- c('flavonoid pathway',
              'trafficking',
              'PR',
              'TF')

fun_ids <- filter(mp_annotation, type %in% our_cats) %>%
    pull('MP_gene_id')

# work on this! we need to group genes by category

mp_DEGmat_sc_plot <- mp_DEGmat_sc[rownames(mp_DEGmat_sc) %in% fun_ids,]

## to do: add descriptions and gaps

an_flav <- data.frame(row.names = rownames(mp_DEGmat_sc_plot), MP_gene_id = rownames(mp_DEGmat_sc_plot)) 

categories <- left_join(an_flav, mp_annotation) %>%
              filter(description != 'PPP1') %>%
              select(c(MP_gene_id, type)) %>%
              arrange(type)

mp_DEGdf_sc_plot <- as.data.frame(mp_DEGmat_sc_plot)
mp_DEGdf_sc_plot$MP_gene_id <- rownames(mp_DEGdf_sc_plot)

mp_DEGmat_sc_plot <- left_join(mp_DEGdf_sc_plot,
                               categories) %>%
                     arrange(type) %>%
                     select(-type)
rownames(mp_DEGmat_sc_plot) <- mp_DEGdf_sc_plot$MP_gene_id

mp_DEGmat_sc_plot <- select(mp_DEGmat_sc_plot, -MP_gene_id) %>% as.matrix()

rownames(categories) <- categories$MP_gene_id
categories <- categories[-1]


pheatmap(mp_DEGmat_sc_plot,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         color = colorRampPalette(c("navy", "white", "#1B9E77"))(50),
         annotation_row = categories,
         cellheight = 7, cellwidth = 7,
         fontsize_row = 8,
         fontsize_col = 8,
         gaps_col=c(12),
         main = 'marchantia single-copy orthologs')



#----------- Clean this later ----




### non 1-1 ortholog relationship

all_og <- read_delim("data/Orthogroups.txt", 
                     delim = ':',
                     col_names = FALSE)

names(all_og) <- c('OG_id', 'genes')

# magic
all_og <- all_og %>%
    mutate(genes = str_split(genes, pattern = ' ')) %>%
    unnest(genes) %>%
    filter(genes != "") %>%
    separate(genes, into = c('gene_id', 'the_rest')) %>%
    select(-the_rest)

## map OG members/size
## OG: n genes in Mp; n genes in Nb

nums_mp <- all_og %>%
    filter(grepl('Mapoly', gene_id)) %>%
    group_by(OG_id) %>% 
    add_count(OG_id) %>%
    ungroup()

nums_nb <- all_og %>%
    group_by(OG_id) %>%
    filter(grepl('Niben', gene_id)) %>% 
    add_count(OG_id) %>%
    ungroup()

nums_both <- full_join(nums_mp, nums_nb, by = 'OG_id') %>%
             rename('n.x' = 'n_Mpoly',
                    'n.y' = 'n_Niben') %>%
             mutate(n_Mpoly = replace_na(n_Mpoly,0),
                    n_Niben = replace_na(n_Niben,0)) 

# overview of gene number changes between Mp and Niben genomes
nums_both %>%
    select(c(OG_id, n_Mpoly, n_Niben)) %>%
    unique() %>%
    ggplot(aes(x = n_Niben, y = n_Mpoly, color = factor(n_Mpoly))) +
    geom_jitter(aes(alpha = 0.1)) +
    xlim(c(0,10)) +
    ylim(c(0,10)) +
    theme(legend.position = 'none')

## Explore flavonoid pathway & trafficking

og_flav <- nums_both %>%
           rename('gene_id.x' = 'MP_gene_id') %>%
           full_join(mp_annotation) %>%
           filter(type == 'flavonoid pathway')

## try heatmaps to visialise members of the same OG, responding (or not during infection)

# these are single-copy orthologs in Niben and Mpoly
View(og_flav %>%
    filter(n_Mpoly == 1 & n_Niben == 1))

# we will try to visualise them





