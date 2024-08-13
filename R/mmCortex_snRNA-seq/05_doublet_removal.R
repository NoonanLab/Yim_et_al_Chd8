# load packages
require(scDblFinder)
require(Seurat)
require(tidyverse)

# load Seurat object
sdata.align <- readRDS('sdata_align_snRNA-seq_label.rds')

# set Idents to scDblFinder-compatible cluster names
Idents(sdata.align) <- sdata.align$seurat_clusters <- sdata.align$integrated_snn_res.0.5

# convert Seurat object to scDblFinder-compatible SCE object
sdata.sce <- as.SingleCellExperiment(sdata.align, assay = "RNA")

# run scDblFinder with cluster-informed sampling
sdata.sce <- scDblFinder(sdata.sce, clusters = 'seurat_clusters', samples = 'sample_id', verbose = TRUE)

# save SCE scDblFinder output
saveRDS(sdata.sce, 'sdata_align_snRNA-seq_scDblFinder-sce.rds')

# append scDblFinder output to input Seurat object's metadata
sdata.align <- as.Seurat(sdata.sce) %>% .@meta.data %>%
  AddMetaData(sdata.align, ., col.name = 'scDblFinder.class')

# remove doublets from Seurat object
sdata.align <- subset(sdata.align, subset = scDblFinder.class == 'singlet')

# saveRDS
saveRDS(sdata.align, 'sdata_align_snRNA-seq_singlet.rds')
