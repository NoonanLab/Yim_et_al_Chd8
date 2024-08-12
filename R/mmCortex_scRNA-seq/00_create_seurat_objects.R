# load packages
require(tidyverse)
require(Seurat)

# Load data --------------------------------------------------------------------
data.dir <- '/path/to/filter_multimap_counts/'
setwd(data.dir)

# get directory names
subdir.names <- list.files(data.dir, pattern = 'filtered', 
                           recursive = T, include.dirs = T)

# load count matrices
sdata.counts <- lapply(subdir.names, Read10X, gene.column = 1)
names(sdata.counts) <- str_remove(subdir.names, '/Solo.out/Gene/filtered')

# load sample info
sample.info <- read_csv('mmCortex_sample_info.csv')

# check names in sample info
identical(sort(names(sdata.counts)), sort(sample.info$sample_id))

# reorder counts
sdata.counts <- sdata.counts[sample.info$sample_id]


# Create gene table ------------------------------------------------------------
gene.table <- read_delim(str_c(subdir.names[1], '/features.tsv.gz'), 
                         col_names = c('mouse_id', 'mouse_name'),
                         col_select = 1:2)

# save gene table
# write_csv(gene.table, '10x_mm10_ensembl98_genes.csv')

# get mitochondrial genes
mt.genes <- filter(gene.table, startsWith(mouse_name, 'mt')) %>% pull(mouse_id)


# Create Seurat Objects --------------------------------------------------------

# create objects
sdata0 <- lapply(sdata.counts, CreateSeuratObject)

# set object metadata
sdata.list <- lapply(names(sdata0),
                     function(sample.id, sample.info, mt.genes) {
  # get object, set project name
  obj <- sdata0[[sample.id]]
  obj@project.name <- sample.id
  # get meta.data
  meta.data <- filter(sample.info, sample_id == sample.id) %>% 
    select(time_point:batch) %>% unlist()
  # scale/center nCount_RNA and nFeature_RNA
  obj[['nCount_RNA.scale']]   <- scale(obj$nCount_RNA)[, 1]
  obj[['nFeature_RNA.scale']] <- scale(obj$nFeature_RNA)[, 1]
  # compute percent mitochondrial genes
  obj[['percent_mito']] <- PercentageFeatureSet(obj, features = mt.genes)
  # store barcodes, append sample ID to cell names
  obj[['barcode']] <- Cells(obj)
  obj <- RenameCells(obj, add.cell.id = sample.id)
  # set meta.data
  obj[['orig.ident']] <- obj[['sample_id']] <- sample.id
  obj <- AddMetaData(obj, metadata = as.list(meta.data))
  obj$sex <- factor(obj$sex, levels = c('male', 'female'))
  obj$genotype <- factor(obj$genotype, levels = c('wt', 'het'))
  # remove empty genes
  obj <- obj[which(rowSums(obj[['RNA']]@counts) > 0), ]
  return(obj)
}, sample.info = sample.info, mt.genes = mt.genes) %>% set_names(names(sdata0))

# save as list
saveRDS(sdata.list, 'sdata_mmCortex_list.rds')
