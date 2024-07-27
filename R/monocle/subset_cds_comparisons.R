# load packages
require(tidyverse)
require(monocle3)
require(Seurat)

# load cell_data_set objects
cds <- readRDS('cds_downsampled_all.rds')

# subset CDS by time point
cds_E12 <- cds[, cds$time_point == 'E12.5']
cds_E14 <- cds[, cds$time_point == 'E14.5']
cds_E16 <- cds[, cds$time_point == 'E16.0']
cds_E17 <- cds[, cds$time_point == 'E17.5']

# clear data
rm(cds)


# subset by comparison ---------------------------------------------------------

# primary trajectory (create list)
cds_E12_list <- list(Primary_Trajectory = cds_E12[, cds_E12$phate_trajectory == 'Primary'])
cds_E14_list <- list(Primary_Trajectory = cds_E14[, cds_E14$phate_trajectory == 'Primary'])
cds_E16_list <- list(Primary_Trajectory = cds_E16[, cds_E16$phate_trajectory == 'Primary'])
cds_E17_list <- list(Primary_Trajectory = cds_E17[, cds_E17$phate_trajectory == 'Primary'])

# cell types
cds_E12_celltypes <- lapply(levels(cds_E12$partition_label), \(x) cds_E12[, cds_E12$phate_trajectory == 'Primary' & cds_E12$partition_label == x])
cds_E14_celltypes <- lapply(levels(cds_E14$partition_label), \(x) cds_E14[, cds_E14$phate_trajectory == 'Primary' & cds_E14$partition_label == x])
cds_E16_celltypes <- lapply(levels(cds_E16$partition_label), \(x) cds_E16[, cds_E16$phate_trajectory == 'Primary' & cds_E16$partition_label == x])
cds_E17_celltypes <- lapply(levels(cds_E17$partition_label), \(x) cds_E17[, cds_E17$phate_trajectory == 'Primary' & cds_E17$partition_label == x])

# add names
names(cds_E12_celltypes) <- levels(cds_E12$partition_label)
names(cds_E14_celltypes) <- levels(cds_E14$partition_label)
names(cds_E16_celltypes) <- levels(cds_E16$partition_label)
names(cds_E17_celltypes) <- levels(cds_E17$partition_label)

# combine lists
cds_E12_list <- c(cds_E12_list, cds_E12_celltypes)
cds_E14_list <- c(cds_E14_list, cds_E14_celltypes)
cds_E16_list <- c(cds_E16_list, cds_E16_celltypes)
cds_E17_list <- c(cds_E17_list, cds_E17_celltypes)

# clear data
rm(list = setdiff(ls(), ls(pattern = 'list')))

   
# calculate pct.wt and pct.het (non-zero rate in each group) -------------------

# Using Seurat FoldChange function
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

# run function
cds_E12_list <- lapply(cds_E12_list, CalculateFC)
cds_E14_list <- lapply(cds_E14_list, CalculateFC)
cds_E16_list <- lapply(cds_E16_list, CalculateFC)
cds_E17_list <- lapply(cds_E17_list, CalculateFC)

# filter genes
cds_E12_filtered <- lapply(cds_E12_list, \(x) { x <- x[fData(x)$keep.gene, ]; fData(x)$keep.gene <- NULL; return(x) })
cds_E14_filtered <- lapply(cds_E14_list, \(x) { x <- x[fData(x)$keep.gene, ]; fData(x)$keep.gene <- NULL; return(x) })
cds_E16_filtered <- lapply(cds_E16_list, \(x) { x <- x[fData(x)$keep.gene, ]; fData(x)$keep.gene <- NULL; return(x) })
cds_E17_filtered <- lapply(cds_E17_list, \(x) { x <- x[fData(x)$keep.gene, ]; fData(x)$keep.gene <- NULL; return(x) })

# save filtered CDS
saveRDS(cds_E12_filtered, 'cds_filtered_E12.rds')
saveRDS(cds_E14_filtered, 'cds_filtered_E14.rds')
saveRDS(cds_E16_filtered, 'cds_filtered_E16.rds')
saveRDS(cds_E17_filtered, 'cds_filtered_E17.rds')