# load packages
require(tidyverse)
require(Seurat)

# load Seurat objects (list)
sdata.list <- readRDS('sdata_mmCortex_snRNA-seq_list.rds')

# filter
sdata.list.filtered <- lapply(sdata.list, \(obj) {
  # define percent_mito filter
  mt.cutoff <- 1
  # filter cells
  obj <- obj %>% 
    subset(percent_mito <= mt.cutoff)  # remove percent_mito > mt.cutoff
  # remove empty genes
  obj <- obj[rowSums(obj[['RNA']]@counts) > 0, ]
})

# saveRDS
saveRDS(sdata.list.filtered, 'sdata_snRNA-seq_filtered_list.rds')