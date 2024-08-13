# load packages
require(tidyverse)
require(Seurat)

# load Seurat object (pre-doublet identification)
sdata.align <- readRDS('sdata_align_snRNA-seq_prelim.rds')

# set Idents (resolution = 0.5)
Idents(sdata.align) <- sdata.align$seurat_clusters <- sdata.align$integrated_snn_res.0.5

# rename Idents in full dataset
sdata.align <- RenameIdents(sdata.align,
                          '0'  = 'L2/3',
                          '2' = 'L4',
                          '4'  = 'L4/5-IT',
                          '9' = 'L5-IT_1',
                          '24' = 'L5-IT_2',
                          '5' = 'L6-IT',
                          '18'  = 'L5/6-NP',
                          '10'  = 'L5-PT',
                          '1'  = 'L6-CT',
                          '16'  = 'mCtx',
                          '25'  = 'Clau',
                          '20'  = 'iN_1',
                          '19' = 'iN_2', 
                          '13'  = 'iN_3',
                          '8'  = 'iN_4',
                          '17' = 'Str',
                          '12'  = 'Ambig',
                          '3'  = 'AG_1',
                          '15'  = 'AG_2',
                          '11'  = 'OPC',
                          '6'  = 'OL_1',
                          '21'  = 'OL_2',
                          '7' = 'MG',
                          '14' = 'Vasc',
                          '26' = 'Glia+UL',
                          '22' = 'Glia+iN',
                          '27' = 'AG+DL',
                          '29' = 'AG+N',
                          '23' = 'MG+UL',
                          '28' = 'MG+DL')

# save in metadata
sdata.align$cluster_label <- Idents(sdata.align)
sdata.align$cluster_label_0 <- str_c(sdata.align$cluster_label, '_',
                                     sdata.align$integrated_snn_res.0.5)

# set levels
cluster.levels <- sdata.align@meta.data %>% 
  arrange(cluster_label, seurat_clusters) %>% 
  pull(cluster_label_0) %>% unique()
sdata.align$cluster_label_0 <- factor(sdata.align$cluster_label_0, levels = cluster.levels)

# saveRDS
saveRDS(sdata.align, 'sdata_align_snRNA-seq_label.rds')

# sdata.align[['SCT']] <- NULL
# sdata.align[['RNA']] <- NULL
# saveRDS(sdata.align, 'sdata_align_snRNA-seq_label_minmin.rds')