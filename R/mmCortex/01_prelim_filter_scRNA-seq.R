# load packages
require(tidyverse)
require(Seurat)

# load Seurat objects (list)
sdata.list <- readRDS('sdata_mmCortex_list.rds')

# filter
sdata.list.filtered <- lapply(sdata.list, \(obj) {
  # define percent_mito filter
  mt.cutoff <- switch(unique(obj$time_point),
                      E12.5 = 2, E14.5 = 2, E16.0 = 3, E17.5 = 3)
  # filter cells
  obj <- obj %>% 
    subset(nFeature_RNA >= 500) %>%    # remove nFeature_RNA < 500
    subset(nCount_RNA.scale <= 4) %>%  # remove nCount_RNA.scale > 4
    subset(percent_mito <= mt.cutoff)  # remove percent_mito > mt.cutoff
  # remove empty genes
  obj <- obj[rowSums(obj[['RNA']]@counts) > 0, ]
})

# get cell cycle genes
cc.genes.tb <- read_csv('seurat_cc_genes_mouse.csv')
cc.genes.mm <- cc.genes.tb %>% group_by(phase) %>% 
  summarize(genes = list(mouse_id)) %>% pull(genes, phase)

# normalize data and compute cell cycle scores
sdata.list.filtered <- lapply(sdata.list.filtered, NormalizeData) 
sdata.list.filtered <- lapply(sdata.list.filtered, CellCycleScoring,
                              s.features   = cc.genes.mm$S, 
                              g2m.features = cc.genes.mm$G2M)

# saveRDS
saveRDS(sdata.list.filtered, 'sdata_filtered_list.rds')
