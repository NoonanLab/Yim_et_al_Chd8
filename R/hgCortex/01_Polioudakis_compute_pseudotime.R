# load packages
require(tidyverse)
require(Seurat)
require(phateR)
require(princurve)
require(scales)

# raw counts and metadata downloaded from solo.bmap.ucla.edu/shiny/webapp/
# load data
load('sc_dev_cortex_geschwind/raw_counts_mat.rdata')
meta.data <- read.csv('sc_dev_cortex_geschwind/cell_metadata.csv', row.names = 1)

# create Seurat object
sdata.hg <- CreateSeuratObject(counts = raw_counts_mat[, rownames(meta.data)], 
                               meta.data = meta.data)

# set 'Cluster' as Idents
Idents(sdata.hg) <- sdata.hg$Cluster
sdata.hg$Cluster <- Idents(sdata.hg)

# compute cell cycle difference
sdata.hg$CC.Difference <- sdata.hg$S_phase_score - sdata.hg$G2M_phase_score

# select clusters for computing PHATE
#     Include: vRG oRG PgS PgG2M IP ExN ExM ExM-U ExDp1 ExDp2
#     Exlude:  InMGE InCGE OPC End Per Mic
sdata.hg$compute_phate <- F
sdata.hg$compute_phate[grep('RG$|^Pg|IP|^Ex', sdata.hg$Cluster)] <- T

# normalize data, find variable features
sdata.hg <- NormalizeData(sdata.hg)
sdata.hg <- FindVariableFeatures(sdata.hg)

# regress out variables
vars.to.regress <- c('Number_UMI', 'Donor', 'Library', 'Percentage_mitochondrial', 'CC.Difference')
sdata.hg <- ScaleData(sdata.hg, vars.to.regress = vars.to.regress)

# compute PCA and UMAP
sdata.hg <- RunPCA(sdata.hg, verbose = F)
sdata.hg <- RunUMAP(sdata.hg, dims = 1:40)

# saveRDS
# saveRDS(sdata.hg, 'sdata_Polioudakis_prelim.rds')

# get PCA
cells.use <- filter(sdata.hg@meta.data, compute_phate) %>% rownames()
dims.use  <- 1:40
pca.embedding <- Embeddings(sdata.hg, 'pca')[cells.use, dims.use]


# Compute PHATE ----------------------------------------------------------------
phate.out <- phate(pca.embedding, knn = 40, t = 40, gamma = 0, seed = 1)

# save PHATE embedding
# saveRDS(phate.out, 'Polioudakis_phate_dims40_k40_t40_gamma0.rds')

# create empty matrix for PHATE
colnames(phate.out$embedding) <- str_c('PHATE_', 1:ncol(phate.out$embedding))
phate.embedding <- matrix(NA, nrow = ncol(sdata.hg), ncol = ncol(phate.out$embedding))
colnames(phate.embedding) <- colnames(phate.out$embedding)
rownames(phate.embedding) <- colnames(sdata.hg)

# store PHATE embedding
phate.embedding[rownames(phate.out$embedding), ] <- phate.out$embedding
sdata.hg[['phate']] <- CreateDimReducObject(embeddings = phate.embedding,
                                            assay = 'RNA', key = 'PHATE_')


# Compute Pseudotime -----------------------------------------------------------

# fit principal curve
curve.out <- phate.out$embedding %>% principal_curve(approx_points = 100)

# rescale principal curve to original PHATE coordinates
pseudotime.curve <- as_tibble(curve.out$s, rownames = 'cell') %>%
  mutate(PHATE_1 = rescale(PHATE_1, to = range(phate.out$embedding[, 1]), from = c(0, 1))) %>%
  mutate(PHATE_2 = rescale(PHATE_2, to = range(phate.out$embedding[, 2]), from = c(0, 1))) %>%
  column_to_rownames('cell') %>% as.matrix()

# create empty matrix for principal curve
curve.embedding <- matrix(NA, nrow = nrow(phate.embedding), ncol = ncol(phate.embedding))
rownames(curve.embedding) <- rownames(phate.embedding)
colnames(curve.embedding) <- paste0('princurve_', 1:ncol(curve.embedding))

# store principal curve in seurat object
curve.embedding[rownames(pseudotime.curve), ] <- pseudotime.curve
sdata.hg[['principal.curve']] <- CreateDimReducObject(embeddings = curve.embedding,
                                                      assay = 'RNA', key = 'princurve_')

# test correlation with SOX2
sox2.data <- FetchData(sdata.hg, 'SOX2', cells = names(curve.out$lambda))
sox2.cor <- cor(unlist(sox2.data), curve.out$lambda)

# create empty pseudotime variable
pseudotime <- rep(NA, ncol(sdata.hg))
names(pseudotime) <- Cells(sdata.hg)

# if SOX2 is correlated with pseudotime, reverse pseudotime (else do not reverse)
if (sox2.cor > 0) {
  pseudotime[(names(curve.out$lambda))] <- rescale(-curve.out$lambda)
} else {
  pseudotime[(names(curve.out$lambda))] <- rescale(curve.out$lambda)
}

# store pseudotime
sdata.hg$pseudotime <- pseudotime

# compute pseudotime bin (n = 20)
nbin <- 20
breaks <- seq(0, 1, length.out = nbin + 1)
pseudotime.cut <- cut(sdata.hg$pseudotime, breaks, include.lowest = T)
pseudotime.bin <- sapply(strsplit(as.character(pseudotime.cut), ','), '[', 1)
pseudotime.bin <- as.numeric(gsub('\\(|\\[', '', pseudotime.bin))

# store pseudotime bin
sdata.hg$pseudotime.cut <- pseudotime.cut
sdata.hg$pseudotime.bin <- pseudotime.bin

# saveRDS
saveRDS(sdata.hg, 'sdata_Polioudakis_pseudotime.rds')