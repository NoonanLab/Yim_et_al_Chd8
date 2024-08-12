# load packages
require(tidyverse)
require(Seurat)

# load Seurat object
sdata.align <- readRDS('sdata_align_RefF_label_minmin.rds')

# for each sample, count number of cells per cell type
cluster.data <- sdata.align %>%
  FetchData(vars = c('sample_id', 'time_point', 'genotype', 'sex', 'cluster_label')) %>%
  group_by(time_point, genotype, sex) %>%
  dplyr::count(sample_id, cluster_label, .drop = F) %>%
  arrange(time_point, cluster_label, genotype, sex) %>%
  dplyr::rename(cell_type = cluster_label, n_cells = n) %>%
  # rename clusters
  mutate(cell_type = case_match(cell_type,
                                'RG' ~ 'Radial Glia',
                                'neuron_early' ~ 'Early Neuron',
                                'neuron_UL' ~ 'UL Neuron',
                                'neuron_DL' ~ 'DL Neuron',
                                'interneuron' ~ 'Interneuron',
                                'vasculature' ~ 'Vasculature',
                                'Cajal_Retzius' ~ 'Cajal-Retzius',
                                'microglia' ~ 'Microglia',
                                .default = cell_type))

# output cluster data
cluster.data.out <- cluster.data %>%
  group_by(sample_id) %>% mutate(pct = n_cells / sum(n_cells)) %>%
  pivot_wider(names_from = cell_type, values_from = c(n_cells, pct)) %>%
  set_names(str_remove(names(.), 'n_cells_'))
  
# write CSV
write_csv(cluster.data.out, 'cluster_data.csv')


# run t-test
cluster.test <- cluster.data %>%
  group_by(genotype) %>% nest(m = c(genotype, sex, sample_id, n_cells)) %>%
  mutate(res = lapply(m, \(x) t.test(n_cells ~ genotype, data = x))) %>%
  mutate(pval = sapply(res, \(x) x$p.value)) %>%
  # compute mean and se per group
  mutate(se = lapply(m, summarize, mean = mean(n_cells), 
                     se = sqrt(var(n_cells)/length(n_cells)))) %>% 
  unnest(se)

# output results
cluster.test.out <- cluster.test %>%
  pivot_wider(names_from = genotype, values_from = c(mean, se)) %>%
  select(time_point:cell_type, starts_with('mean'), starts_with('se'), pval)

# write CSV
write_csv(cluster.test.out, 'cluster_ttest.csv')
