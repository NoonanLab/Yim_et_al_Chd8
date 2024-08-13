# load packages
require(Seurat)

# load Seurat objects (list)
sdata.list <- readRDS('sdata_snRNA-seq_filtered_list.rds')

# split by batch
sdata.batch <- merge(sdata.list[[1]], sdata.list[2:length(sdata.list)])
sdata.batch <- SplitObject(sdata.batch, split.by = 'batch')

# SCTransform
# Note: SCTransform replaces NormalizeData, ScaleData, and FindVariableFeatures
sdata.batch <- lapply(sdata.batch, function(obj) {
  obj@project.name <- unique(obj$batch)
  message('[ SCTranform: Batch ', obj@project.name, ' ]')
  vars.to.regress <- c('percent_mito')
  obj <- SCTransform(obj, vars.to.regress = vars.to.regress)
  message(); return(obj)
})

# saveRDS
saveRDS(sdata.batch, 'sdata_snRNA-seq_batch_SCT.rds')
