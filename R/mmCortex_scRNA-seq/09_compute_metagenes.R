# load packages
require(tidyverse)
require(seriation)

# load data
bin.paa.sax <- readRDS('bin_paa8_alpha9.rds')

# get SAX
sax.tb <- bin.paa.sax %>% 
  filter(!na.expression) %>% # remove empty genes
  transmute(time_point, genotype, mouse_id, mouse_name, sax9) %>% 
  unnest_wider(col = sax9, names_sep = '_')

# subset wildtype, no E12.5
sax.tb.wt <- filter(sax.tb, time_point != 'E12.5' & genotype == 'wt')


# cluster trajectories ---------------------------------------------------------
k <- 25
set.seed(1); sax.cl <- kmeans(select(sax.tb.wt, starts_with('sax')), 
                              centers = k, iter.max = 20)

# compute correlation distance and linkage tree
r <- cor(t(sax.cl$centers))
d <- as.dist(1 - r)
hc <- hclust(d, method = 'complete') %>% seriation:::reorder.hclust(d)

# cut linkage tree
hc.cut <- cutree(hc, h = 0.6)

# re-number clusters
hc.tb <- tibble(h.clust = unique(hc.cut[hc$order]), h.new = 1:max(hc.cut)) %>%
  left_join(tibble(k.clust = as.integer(names(hc.cut)), h.clust = hc.cut), by = 'h.clust') %>%
  left_join(tibble(k.clust = hc$order, k.order = 1:k), by = 'k.clust') %>%
  # label metagenes, A through J
  left_join(tibble(h.new = 1:10, metagene = str_split_1('ABCDEFGHIJ', '')), by = 'h.new')

# make cluster table
clusters.tb <- sax.tb.wt %>% add_column(k.clust = sax.cl$cluster) %>% 
  left_join(hc.tb, by = 'k.clust') %>% arrange(k.order)

# compute metagene centers
centers.new <- clusters.tb %>% group_by(metagene) %>% 
  summarize(across(starts_with('sax'), mean))

# compute correlation from all genes (SAX) to hclust centers -- exclude E12.5
# NOTE: genes with flat line trajectories will return NA correlation
sax.cor <- filter(sax.tb, time_point != 'E12.5') %>% 
  unite('rowname', time_point, genotype, mouse_id, mouse_name) %>%
  column_to_rownames() %>% t() %>%
  cor(t(column_to_rownames(centers.new, 'metagene')))

# find maximum correlation between trajectory and metagene center
metagenes.cor <- as_tibble(sax.cor, rownames = 'rowname') %>%
  separate(rowname, into = c('time_point', 'genotype', 'mouse_id', 'mouse_name'), sep = '_') %>%
  add_column(r = apply(sax.cor, 1, max),
             metagene = apply(sax.cor, 1, \(x) names(which(x == max(x))))) %>%
  left_join(select(clusters.tb, time_point:mouse_id, k.order),
            by = c('time_point', 'genotype', 'mouse_id')) %>%
  relocate(time_point:mouse_name, metagene, r)

# output metagene assignments
metagenes.out <- metagenes.cor %>% select(time_point:r) %>%
  pivot_wider(names_from = genotype, values_from = c(metagene, r)) %>%
  arrange(time_point, metagene_wt, -r_wt) %>%
  rename(cor_wt = r_wt, cor_het = r_het)

# write CSV
write_csv(metagenes.out, 'metagenes_excludeE12.csv')

# save variables
# save(sax.cl, hc, clusters.tb, centers.new, file = 'metagene_cluster_data.Rdata')
write_csv(centers.new, 'metagene_centers.csv')
