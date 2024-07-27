# load packages
require(tidyverse)
require(Seurat)
require(phateR) # version 1.0.7
require(princurve)
require(scales)

# load Seurat object
sdata.align <- readRDS('sdata_align_RefF_label.rds')

# get PCA
cells.use <- filter(sdata.align@meta.data, compute_phate) %>% rownames()
dims.use  <- 1:50
pca.embedding <- Embeddings(sdata.align, 'pca')[cells.use, dims.use]


# Compute PHATE  ---------------------------------------------------------------
phate.out <- phate(pca.embedding, knn = 80, t = 80, gamma = 0, n.jobs = -1, seed = 1)

# save PHATE embedding
saveRDS(phate.out, 'phate_dims50_k80_t80_gamma0.rds')

# create empty matrix for PHATE
phate.embedding <- matrix(NA, nrow = ncol(sdata.align), ncol = ncol(phate.out$embedding))
rownames(phate.embedding) <- colnames(sdata.align)
colnames(phate.embedding) <- str_c('PHATE_', 1:ncol(phate.embedding))

# store PHATE embedding
phate.embedding[rownames(phate.out$embedding), ] <- phate.out$embedding
sdata.align[['phate']] <- CreateDimReducObject(embeddings = phate.embedding,
                                               assay = 'integrated', key = 'PHATE_')

# label primary and secondary PHATE trajectories
sdata.align@meta.data <- sdata.align@meta.data %>%
  mutate(phate_trajectory = case_when(cluster_label_0 == 'RG_10' ~ 'Secondary',
                                      compute_phate ~ 'Primary'))


# Compute Pseudotime -----------------------------------------------------------

# subset primary trajectory from PHATE embedding
phate.subset <- FetchData(sdata.align, 
                          vars = c('PHATE_1', 'PHATE_2', 'phate_trajectory')) %>%
  filter(phate_trajectory == 'Primary') %>% select(-phate_trajectory)

# rescale and fit principal curve
curve.out <- phate.subset %>%
  mutate(across(everything(), rescale)) %>% 
  as.matrix() %>% principal_curve(approx_points = 100)

# rescale principal curve to original PHATE coordinates
pseudotime.curve <- as_tibble(curve.out$s, rownames = 'cell') %>%
  mutate(PHATE_1 = rescale(PHATE_1, to = range(phate.subset[, 1]), from = c(0, 1)),
         PHATE_2 = rescale(PHATE_2, to = range(phate.subset[, 2]), from = c(0, 1))) %>%
  column_to_rownames('cell') %>% as.matrix()

# create empty matrix for principal curve
curve.embedding <- matrix(NA, nrow = nrow(phate.embedding), ncol = ncol(phate.embedding))
rownames(curve.embedding) <- rownames(phate.embedding)
colnames(curve.embedding) <- paste0('princurve_', 1:ncol(curve.embedding))

# store principal curve in seurat object
curve.embedding[rownames(pseudotime.curve), ] <- pseudotime.curve
sdata.align[['principal.curve']] <- CreateDimReducObject(embeddings = curve.embedding,
                                                         assay = 'integrated', key = 'princurve_')

# test correlation with Sox2 = ENSMUSG00000074637
sox2.data <- FetchData(sdata.align, vars = 'ENSMUSG00000074637', 
                       cells = names(curve.out$lambda))
sox2.cor <- cor(sox2.data[, 1], curve.out$lambda)

# create empty pseudotime variable
pseudotime <- rep(NA, ncol(sdata.align))
names(pseudotime) <- Cells(sdata.align)

# if Sox2 is correlated with pseudotime, reverse pseudotime
if (sox2.corr > 0) {
  pseudotime[(names(curve.out$lambda))] <- rescale(-curve.out$lambda)
} else {
  pseudotime[(names(curve.out$lambda))] <- rescale(curve.out$lambda)
}

# store pseudotime
sdata.align$pseudotime <- pseudotime

# compute pseudotime bin (n = 20)
nbin <- 20
breaks <- seq(0, 1, length.out = nbin + 1)
pseudotime.cut <- cut(sdata.align$pseudotime, breaks, include.lowest = T)
pseudotime.bin <- sapply(strsplit(as.character(pseudotime.cut), ','), '[', 1)
pseudotime.bin <- as.numeric(gsub('\\(|\\[', '', pseudotime.bin))

# store pseudotime bin
sdata.align$pseudotime.cut <- pseudotime.cut
sdata.align$pseudotime.bin <- pseudotime.bin


# saveRDS ----------------------------------------------------------------------

saveRDS(sdata.align, 'sdata_align_pseudotime.rds')

sdata.align[['SCT']] <- NULL
saveRDS(sdata.align, 'sdata_align_pseudotime_min.rds')

sdata.align[['RNA']] <- NULL
saveRDS(sdata.align, 'sdata_align_pseudotime_minmin.rds')
