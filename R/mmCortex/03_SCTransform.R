# load packages
require(Seurat)

# load Seurat objects (list)
sdata.list <- readRDS('sdata_filtered_list.rds')

# split by batch
sdata.batch <- merge(sdata.list[[1]], sdata.list[2:length(sdata.list)])
sdata.batch <- SplitObject(sdata.batch, split.by = 'batch')

# SCTransform
# Note: SCTransform replaces NormalizeData, ScaleData, and FindVariableFeatures
sdata.batch <- lapply(sdata.batch, function(obj) {
  obj@project.name <- unique(obj$batch)
  message('[ SCTranform: Batch ', obj@project.name, ' ]')
  # for scRNA-seq dataset
  obj$CC.Difference <- obj$S.Score - obj$G2M.Score  
  vars.to.regress <- c('CC.Difference', 'percent_mito') 
  # for snRNA-seq dataset
  vars.to.regress <- c('percent_mito') 
  # for all datasets
  obj <- SCTransform(obj, vars.to.regress = vars.to.regress)
  message(); return(obj)
})

# saveRDS
saveRDS(sdata.batch, 'sdata_batch_SCT.rds')
