# load packages
require(tidyverse)
require(Seurat)
require(pbapply)

# load Seurat object and gene names
sdata.hg <- readRDS('sdata_Polioudakis_pseudotime.rds')
gene.names <- read_csv('Polioudakis_ensembl88_genes_edit.csv')

# get metadata, remove cells with NA pseudotime value
meta.data <- sdata.hg@meta.data %>% 
  rownames_to_column(var = 'cell.name') %>% filter(!is.na(pseudotime))

# get data
expr.data <- GetAssayData(sdata.hg)[, meta.data$cell.name]

# apply function over each row (gene)
#   tapply(X, INDEX = pseudotime.bin, FUN = mean)

# compute pct.detected per bin
detected <- expr.data > 0
m0 <- pbapply(detected, 1, tapply, meta.data$pseudotime.bin, mean) 

# filter by pct.detected, max bin per gene
keep.genes <- names(which(apply(m0, 2, max) >= 0.05))

# compute mean.expr per bin
m1 <- pbapply(expr.data[keep.genes, ], 1, tapply, meta.data$pseudotime.bin, mean)

# add gene ids, fix column names
bin.data.out <- t(m1) %>% 
  as_tibble(rownames = 'human_name') %>%
  set_names(str_replace(make.names(names(.)), 'X', 'bin_')) %>%
  left_join(gene.names, by = 'human_name') %>% relocate(human_id)

# saveRDS
saveRDS(bin.data.out, 'Polioudakis_bin_data.rds')