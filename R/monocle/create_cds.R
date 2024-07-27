# load packages
require(tidyverse)
require(Seurat)
require(monocle3)

# load data
sdata.align <- readRDS('sdata_align_pseudotime_minmin.rds')
downsampled_matrix <- readRDS('downsampled_counts_matrix.rds')

# load gene names table
gene.names <- read_csv('gene_names.csv')

# add downsampled counts, set default assay
sdata.new <- sdata.align
sdata.new[['downsampled']] <- CreateAssayObject(counts = downsampled_matrix)
DefaultAssay(sdata.new) <- 'downsampled'
sdata.new <- DietSeurat(sdata.new, assays = 'downsampled')


# convert to SCE/CDS -----------------------------------------------------------
temp.sce <- as.SingleCellExperiment(sdata.new)
logcounts(temp.sce) <- NULL

# get gene metadata
gene.metadata <- data.frame(mouse_id = rownames(temp.sce)) %>%
  left_join(gene.names, by = 'mouse_id') %>% column_to_rownames('mouse_id') %>%
  dplyr::rename(gene_short_name = mouse_name)

# create CDS
cds <- new_cell_data_set(expression_data = counts(temp.sce),
                         cell_metadata = colData(temp.sce),
                         gene_metadata = gene.metadata)

# label cell types
cds$partition_label <- case_match(cds$pseudotime.bin,
                                  c(0.00, 0.05, 0.10, 0.15, 0.20) ~ 'Radial_Glia',
                                  c(0.25, 0.30, 0.35, 0.40, 0.45) ~ 'IP',
                                  c(0.50, 0.55) ~ 'Early_Neurons',
                                  c(0.60, 0.65, 0.70, 0.75) ~ 'UL_Neurons',
                                  c(0.80, 0.85, 0.90) ~ 'DL_Neurons',
                                  c(0.95) ~ 'Subplate') %>%
  factor(levels = c('Radial_Glia', 'IP', 'Early_Neurons', 
                    'UL_Neurons', 'DL_Neurons', 'Subplate'))

# save CDS
saveRDS(cds, 'cds_downsampled_all.rds')