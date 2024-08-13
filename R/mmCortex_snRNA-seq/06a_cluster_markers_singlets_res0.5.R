# load packages
require(tidyverse)
require(Seurat)

# load Seurat object (singlet-only dataset) and gene names
sdata.align <- readRDS('sdata_align_snRNA-seq_singlet.rds') 
gene.names <- read_csv('gene_names.csv')

# set default assay to 'RNA'
DefaultAssay(sdata.align) <- 'RNA'

# set Idents to cluster labels
Idents(sdata.align) <- sdata.align$cluster_label

# reduce object size
sdata.align[['integrated']] <- NULL

# check, use cluster resolution = 0.5
# identical(Idents(sdata.align), sdata.align$integrated_snn_res.0.5)
# Idents(sdata.align) <- sdata.align$seurat_clusters <- sdata.align$integrated_snn_res.0.5

# compute cluster markers (~2 hr)
cluster.markers <- FindAllMarkers(sdata.align, assay = 'RNA', logfc.threshold = 0.25)

# output table of cluster markers with gene names
cluster.markers.out <- as_tibble(cluster.markers) %>%
  rename(mouse_id = gene) %>%
  group_by(cluster) %>% mutate(rank = rank(-avg_log2FC)) %>%
  relocate(cluster, p_val:p_val_adj, rank) %>%
  left_join(gene.names, by = 'mouse_id')

write_csv(cluster.markers.out, 'integrated_clusterMarkers_singlet_res0.5.csv') 