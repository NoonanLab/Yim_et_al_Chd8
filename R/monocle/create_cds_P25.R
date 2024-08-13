# load packages
require(tidyverse)
require(Seurat)
require(monocle3)

# load data
sdata.P25 <- readRDS('sdata_align_snRNA-seq_singlet.rds')
sdata.P25 <- DietSeurat(sdata.P25, assays = 'RNA')

# set Idents to cluster labels
Idents(sdata.P25) <- sdata.P25$cluster_label

# load gene names table
gene.names <- read_csv('gene_names.csv')


# convert to SCE/CDS -----------------------------------------------------------
temp.sce <- as.SingleCellExperiment(sdata.P25)
logcounts(temp.sce) <- NULL

# get gene metadata
gene.metadata <- data.frame(mouse_id = rownames(temp.sce)) %>%
  left_join(gene.names, by = 'mouse_id') %>% column_to_rownames('mouse_id') %>%
  dplyr::rename(gene_short_name = mouse_name)

# create CDS
cds_P25 <- new_cell_data_set(expression_data = counts(temp.sce),
                             cell_metadata = colData(temp.sce),
                             gene_metadata = gene.metadata)

# fix meta.data
cds_P25$genotype <- relevel(as.factor(cds_P25$genotype), ref = 'wt')
cds_P25$sex <- relevel(as.factor(cds_P25$sex), ref = 'male')

# split by cluster labels passing filtering steps
clusters.exclude <- c('L5-IT_2', 'mCtx', 'Clau', 'Str', 'N?', 
                      'Glia+UL', 'Glia+iN', 'AG+DL', 'AG+N', 'MG+UL', 'MG+DL')
clusters.keep <- setdiff(levels(cds_P25$cluster_label), clusters.exclude)
cds_P25_list <- lapply(clusters.keep, \(x) cds_P25[, cds_P25$cluster_label == x]) %>% set_names(clusters.keep)


# calculate pct.wt and pct.het (non-zero rate in each group) -------------------

# Using Seurat FoldChange function
# NOTE: positive avg_log2FC indicate that a gene is more highly expressed in the first group
CalculateFC <- function(cds, slot = 'data', min.pct = 0.1, ...) {
  m <- counts(cds) %>% CreateAssayObject() %>% NormalizeData() %>%
    FoldChange(cells.1 = colnames(cds)[cds$genotype == 'het'],
               cells.2 = colnames(cds)[cds$genotype == 'wt'],
               slot = slot, ...)
  fData(cds)$avg_log2FC <- m$avg_log2FC
  fData(cds)$pct.wt     <- m$pct.2
  fData(cds)$pct.het    <- m$pct.1
  fData(cds)$keep.gene  <- pmax(m$pct.1, m$pct.2) >= min.pct
  return(cds)
}


# run function and filter genes
cds_P25_filtered <- lapply(cds_P25_list, CalculateFC) %>% 
  lapply(\(x) x[fData(x)$keep.gene, ]) %>%
  lapply(\(x) { fData(x)$keep.gene <- NULL; return(x) })

# save filtered CDS
saveRDS(cds_P25_filtered, 'cds_filtered_P25.rds')
