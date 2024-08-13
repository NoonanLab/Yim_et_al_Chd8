# load packages
require(tidyverse)
require(Seurat)

# load Seurat object
sdata.align <- readRDS('sdata_align_snRNA-seq_singlet.rds')

# for each sample, count number of singlets per cell type
cluster.data <- sdata.align %>%
  FetchData(vars = c('sample_id', 'genotype', 'sex', 'batch', 'cluster_label')) %>%
  mutate(cluster_label = droplevels(cluster_label)) %>%
  group_by(batch, genotype, sex) %>%
  dplyr::count(sample_id, cluster_label, .drop = F) %>%
  arrange(cluster_label, batch, genotype, sex) %>%
  dplyr::rename(cell_type = cluster_label, n_nuclei = n) 
  
# output cluster data
cluster.data.out <- cluster.data %>%
  group_by(sample_id) %>% mutate(pct = n_nuclei / sum(n_nuclei)) %>%
  pivot_wider(names_from = cell_type, values_from = c(n_nuclei, pct)) %>%
  set_names(str_remove(names(.), 'n_nuclei_'))

# write CSV
write_csv(cluster.data.out, 'cluster_data_singlet_all.csv')

# prepare to subset clusters to those passing cluster-based filtering steps, excluding clusters that:
	# (1) contain < 200 singlet nuclei; 
	# (2) disproportionately consist of nuclei from one batch (≥75%); or
	# (3) cannot be assigned a cell type identity, due to ambiguous marker genes
clusters.exclude <- c('L5-IT_2', 'mCtx', 'Clau', 'Str', 'N?',
                      'Glia+UL', 'Glia+iN', 'AG+DL', 'AG+N', 'MG+UL', 'MG+DL')
clusters.keep <- setdiff(levels(sdata.align$cluster_label), clusters.exclude)

# for each sample, count number of singlets per cell type that passed cluster filtering
cluster.data.filter <- sdata.align %>%
  FetchData(vars = c('sample_id', 'genotype', 'sex', 'cluster_label')) %>%
  filter(cluster_label %in% clusters.keep) %>%
  mutate(cluster_label = droplevels(cluster_label)) %>%
  group_by(genotype, sex) %>%
  dplyr::count(sample_id, cluster_label, .drop = F) %>%
  arrange(cluster_label, genotype, sex) %>%
  dplyr::rename(cell_type = cluster_label, n_nuclei = n) 

# output filtered cluster data
cluster.data.filter.out <- cluster.data.filter %>%
  group_by(sample_id) %>% mutate(pct = n_nuclei / sum(n_nuclei)) %>%
  pivot_wider(names_from = cell_type, values_from = c(n_nuclei, pct)) %>%
  set_names(str_remove(names(.), 'n_nuclei_'))
  
# write CSV for filtered dataset
write_csv(cluster.data.filter.out, 'cluster_data_singlet_filtered.csv')

# run t-test on clusters passing all filtering criteria
cluster.test <- cluster.data.filter %>%
  group_by(genotype) %>% nest(m = c(genotype, sex, sample_id, n_nuclei)) %>%
  mutate(res = lapply(m, \(x) t.test(n_nuclei ~ genotype, data = x))) %>%
  mutate(pval = sapply(res, \(x) x$p.value)) %>%
  # compute mean and se per group
  mutate(se = lapply(m, summarize, mean = mean(n_nuclei), 
                     se = sqrt(var(n_nuclei)/length(n_nuclei)))) %>% 
  unnest(se)

# output results
cluster.test.out <- cluster.test %>%
  pivot_wider(names_from = genotype, values_from = c(mean, se)) %>%
  select(cell_type, starts_with('mean'), starts_with('se'), pval)

# write CSV
write_csv(cluster.test.out, 'cluster_ttest_singlet_filtered.csv')
