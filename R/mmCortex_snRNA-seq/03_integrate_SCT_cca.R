# load packages
require(Seurat)

# load Seurat objects (list)
sdata.batch <- readRDS('sdata_snRNA-seq_batch_SCT.rds') # output from 'SCTransform.R'

# select integration features and prepare object
features    <- SelectIntegrationFeatures(sdata.batch) # default, nfeatures = 2000
sdata.batch <- PrepSCTIntegration(sdata.batch, anchor.features = features)

# compute anchors, CCA, without reference
batch.anchors <- FindIntegrationAnchors(object.list = sdata.batch, 
                                        anchor.features = features,
                                        normalization.method = 'SCT',
                                        reduction = 'cca', dims = 1:30) # default

# integrate data
sdata.align <- IntegrateData(anchorset = batch.anchors,
                             normalization.method = 'SCT',
                             dims = 1:30)

# fix meta.data
sdata.align$genotype <- relevel(factor(sdata.align$genotype), ref = 'wt')
sdata.align$sex <- relevel(factor(sdata.align$sex), ref = 'male')


# Prelim Analysis --------------------------------------------------------------
# DO NOT scale after integrating with SCT
dims.use <- 1:17

sdata.align <- RunPCA(sdata.align, verbose = F) # default npcs = 50
sdata.align <- RunUMAP(sdata.align, dims = dims.use, min.dist = 0.3, n.neighbors = 40)

# compute prelim clusters
sdata.align <- FindNeighbors(sdata.align, dims = dims.use) # default: k.param = 20, prune.SNN = 1/15
sdata.align <- FindClusters(sdata.align, resolution = seq(0.2, 0.8, by = 0.1))

# set Ident, resolution = 0.5
Idents(sdata.align) <- sdata.align$seurat_clusters <- sdata.align$integrated_snn_res.0.5

# saveRDS
saveRDS(sdata.align, 'sdata_align_snRNA-seq_prelim.rds')