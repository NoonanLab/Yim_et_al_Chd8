# load packages
require(tidyverse)
require(Seurat)
require(pbapply)
require(tictoc)

# load data and gene names
sdata.align  <- readRDS('sdata_align_pseudotime_minmin.rds')
m.downsample <- readRDS('downsampled_counts_matrix.rds')
gene.names <- read_csv('gene_names.csv')

# get metadata, remove cells with NA pseudotime value
meta.data <- sdata.align@meta.data %>% 
  rownames_to_column(var = 'cell.name') %>% filter(!is.na(pseudotime))

# subset and normalize downsampled data
norm.data <- m.downsample[, meta.data$cell.name] %>%
  CreateAssayObject() %>% NormalizeData() %>% GetAssayData()

# split data by time_point and genotype
group.data <- meta.data %>% group_by(time_point, genotype) %>%
  select(cell.name, pseudotime.bin) %>%
  nest(meta.data = c(cell.name, pseudotime.bin)) %>%
  mutate(expr.data = lapply(meta.data, \(x) norm.data[, x$cell.name]),
         detected  = lapply(expr.data, \(x) x > 0))

# for each time_point and genotype ... (~1 hr 15 min)
tic(); bin.data <- mapply(function(meta.data, expr.data, detected) { 
  # apply function over each row (gene)
  #   tapply(X, INDEX = pseudotime.bin, FUN = mean)
  m0 <- pbapply(expr.data[, meta.data$cell.name], 1,
                tapply, meta.data$pseudotime.bin, mean) # mean.expr per bin
  m1 <- pbapply(detected[, meta.data$cell.name], 1,
                tapply, meta.data$pseudotime.bin, mean) # pct.detected per bin
  m2 <- tibble(mean.expr = list(t(m0)), pct.detected = list(t(m1)))
}, meta.data = group.data$meta.data, expr.data = group.data$expr.data, 
   detected = group.data$detected, SIMPLIFY = F) %>% bind_rows(); toc()

# combine tables
bin.data.out <- select(group.data, time_point, genotype) %>% bind_cols(bin.data)

# filter tables
bin.data.filtered <- bin.data.out %>%
  # get max(pct.detected) for each gene across all bins, select max.pct >= 0.05
  transmute(mean.expr, max.pct = lapply(pct.detected, apply, 1, max),
            genes.keep = lapply(max.pct, \(x) names(which(x >= 0.05)))) %>%
  # get union of wild type and het genes (gene is present in either group)
  group_by(time_point) %>% 
  mutate(genes.keep = list(sort(unique(unlist(genes.keep))))) %>% 
  rowwise() %>% transmute(time_point, genotype, 
                          x = list(as.data.frame(mean.expr[genes.keep, ]))) %>%
  # add gene names, fix column names
  ungroup() %>% mutate(x = lapply(x, rownames_to_column, 'mouse_id')) %>% 
  unnest(x, names_repair = make.names) %>% 
  set_names(str_replace(names(.), 'X', 'bin_')) %>%
  # add gene names
  left_join(gene.names, by = 'mouse_id') %>% 
  relocate(time_point:mouse_id, mouse_name)

# saveRDS
saveRDS(bin.data.out, 'bin_data.rds')
saveRDS(bin.data.filtered, 'bin_data_filtered.rds')
