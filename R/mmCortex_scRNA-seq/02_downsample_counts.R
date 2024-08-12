# load packages
require(tidyverse)
require(Seurat)
require(DropletUtils)

# load Seurat objects (list)
sdata.list <- readRDS('sdata_filtered_list.rds')

# get meta.data
meta.data <- lapply(sdata.list, \(obj) obj@meta.data) %>% bind_rows()

# calculate scale factors
x <- meta.data %>% 
  group_by(time_point, sample_id) %>% 
  summarize(median.nCount_RNA = median(nCount_RNA)) %>%
  mutate(scale_factor = min(median.nCount_RNA)/median.nCount_RNA)

# get counts
x.counts <- lapply(sdata.list, \(obj) obj[['RNA']]@counts)[x$sample_id]

# downsample counts
sdata.counts.downsample <- mapply(\(x, prop, seed.use) {
  set.seed(seed.use)
  m <- downsampleMatrix(x = x, prop = prop, bycol = F)
}, x = x.counts, prop = x$scale_factor, seed.use = 0)

# add empty genes
genes <- sort(unique(unlist(sapply(sdata.counts.downsample, rownames))))
sdata.counts.downsample <- lapply(sdata.counts.downsample, \(m){
  genes.empty <- setdiff(genes, rownames(m))
  m.empty <- matrix(0, nrow = length(genes.empty), ncol = ncol(m), 
                    dimnames = list(genes.empty, colnames(m)))
  m.out <- rbind(m, m.empty)[genes, ]
})

# combine count matrices, remove empty genes
m.downsample <- do.call(cbind, sdata.counts.downsample)
m.downsample <- m.downsample[rowSums(m.downsample) > 0, ]

# saveRDS
saveRDS(m.downsample, 'downsampled_counts_matrix.rds')
