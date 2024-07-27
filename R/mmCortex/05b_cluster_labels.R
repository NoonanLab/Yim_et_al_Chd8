# load packages
require(tidyverse)
require(Seurat)

# load Seurat object
sdata.align <- readRDS('sdata_align_RefF_prelim.rds')

# set Idents (resolution = 0.4)
Idents(sdata.align) <- sdata.align$seurat_clusters <- sdata.align$integrated_snn_res.0.4

# rename Idents
sdata.align <- RenameIdents(sdata.align,
                            '0'  = 'RG',
                            '3'  = 'RG',
                            '10' = 'RG',
                            '11' = 'RG',
                            '4'  = 'IP',
                            '5'  = 'IP',
                            '2'  = 'neuron_early',
                            '1'  = 'neuron_UL',
                            '6'  = 'neuron_UL',
                            '8'  = 'neuron_UL',
                            '7'  = 'neuron_DL',
                            '9'  = 'interneuron',
                            '14' = 'RBC',
                            '13' = 'vasculature',
                            '15' = 'vasculature',
                            '12' = 'Cajal_Retzius',
                            '17' = 'OPC',       # Other_Olig
                            '16' = 'microglia') # Other_C1q

# save in metadata
sdata.align$cluster_label <- Idents(sdata.align)
sdata.align$cluster_label_0 <- str_c(sdata.align$cluster_label, '_',
                                     sdata.align$integrated_snn_res.0.4)

# set levels
cluster.levels <- sdata.align@meta.data %>% 
  arrange(cluster_label, seurat_clusters) %>% 
  pull(cluster_label_0) %>% unique()
sdata.align$cluster_label_0 <- factor(sdata.align$cluster_label_0, levels = cluster.levels)

# select clusters to compute PHATE (other clusters omitted)
sdata.align$compute_phate <- F
sdata.align$compute_phate[grep('RG|IP|^Neuron', sdata.align$cluster_label)] <- T

# saveRDS
saveRDS(sdata.align, 'sdata_align_RefF_label.rds')

# sdata.align[['SCT']] <- NULL
# sdata.align[['RNA']] <- NULL
# saveRDS(sdata.align, 'sdata_align_RefF_label_minmin.rds')