# load packages
require(tidyverse)
require(seriation)
require(jmotif)

# load data
bin.data <- readRDS('Polioudakis_bin_data.rds')

# load mouse metagene centers
centers.mm <- read_csv('mmCortex_metagene_centers.csv')

# scale data
bin.data.scale <- bin.data %>% select(-human_id) %>% 
  column_to_rownames('human_name') %>% t() %>% scale() %>% t() %>%
  as_tibble(rownames = 'human_name') %>%
  left_join(select(bin.data, human_name, human_id), by = 'human_name') %>%
  relocate(human_id)

# SAX variables
paa.size <- 8
alphabet.size <- 9

# compute SAX transform
bin.paa8.sax9 <- bin.data.scale %>% rowwise() %>%
  transmute(human_id, human_name, 
            bin.mean.scale = list(c_across(starts_with('bin')))) %>%
  mutate(paa = list(paa(bin.mean.scale, paa_num = paa.size)),
         sax.word = list(series_to_chars(paa, a_size = alphabet.size)),
         sax9 = list(letters_to_idx(sax.word))) %>%
  mutate(sax.word = str_c(sax.word, collapse = '')) %>%
  ungroup()

# saveRDS
#saveRDS(bin.paa8.sax9, 'Polioudakis_paa8_alpha9.rds')

# get SAX
sax.tb <- bin.paa8.sax9 %>% 
  select(human_name, sax9) %>% unnest_wider(sax9, names_sep = '_')


# cluster trajectories ---------------------------------------------------------
k <- 25
set.seed(1); sax.cl <- kmeans(select(sax.tb, starts_with('sax')), 
                              centers = k, iter.max = 20)

# compute correlation distance and linkage tree
r <- cor(t(sax.cl$centers))
d <- as.dist(1 - r)
hc <- hclust(d, method = 'complete') %>% seriation:::reorder.hclust(d)

# cut linkage tree
hc.cut <- cutree(hc, h = 0.6)

# re-number clusters
hc.tb <- tibble(h.clust = unique(hc.cut[hc$order]), 
                h.label = str_c('h', 1:max(hc.cut))) %>%
  left_join(tibble(k.clust = as.integer(names(hc.cut)), h.clust = hc.cut), by = 'h.clust') %>%
  left_join(tibble(k.clust = hc$order, k.order = 1:k), by = 'k.clust')

# make cluster table
clusters.tb <- sax.tb %>% add_column(k.clust = sax.cl$cluster) %>% 
  left_join(hc.tb, by = 'k.clust') %>% arrange(k.order) %>%
  mutate(h.label = factor(h.label, levels = unique(h.label)))

# compute metagene centers
centers.new <- clusters.tb %>% group_by(h.label) %>% 
  summarize(across(starts_with('sax'), mean))

# compute correlation from all genes (SAX) to hclust centers
# NOTE: genes with flat line trajectories will return NA correlation
sax.cor <- sax.tb %>% column_to_rownames('human_name') %>% t() %>%
  cor(t(column_to_rownames(centers.new, 'h.label')))

# find maximum correlation between trajectory and metagene center
metagenes.cor <- as_tibble(sax.cor, rownames = 'human_name') %>%
  add_column(r = apply(sax.cor, 1, max),
             h.label = apply(sax.cor, 1, \(x) names(which(x == max(x))))) %>%
  left_join(select(bin.data, human_name, human_id), by = 'human_name') %>%
  relocate(human_id,human_name, h.label, r)


# compare human metagenes to mouse metagenes -----------------------------------

# compute correlation between metagene centers
centers.cor <- cor(t(column_to_rownames(centers.new, 'h.label')),
                   t(column_to_rownames(centers.mm, 'metagene')))

# label human metagenes
centers.label <- centers.cor %>% as_tibble(rownames = 'h.label') %>% 
  pivot_longer(-h.label, names_to = 'mouse_metagene', values_to = 'r') %>%
  # top human metagene per mouse metagene
  group_by(mouse_metagene) %>% mutate(top.m = rank(-r) == 1) %>%
  # top mouse metagene per human metagene
  group_by(h.label) %>% mutate(n = sum(top.m), rank.r = rank(-r)) %>%
  filter(top.m | (n == 0 & rank.r == 1)) %>% slice_max(r) %>%
  # manual label
  mutate(human_metagene = case_match(h.label, 'h8' ~ 'K', 'h6' ~ 'L', 
                                     .default = mouse_metagene),
         human_metagene = str_c(human_metagene, "'")) %>%
  arrange(human_metagene)

# output metagene assignments
metagenes.out <- metagenes.cor %>%
  left_join(select(centers.label, h.label, human_metagene), by = 'h.label') %>%
  transmute(human_id, human_name, h.label, human_metagene, cor = r) %>%
  arrange(human_metagene, -cor) %>%
  mutate(ASD = human_id %in% genes.asd$human_id,
         DDD = human_id %in% genes.ddd$human_id,
         CHD8_Target = human_id %in% genes.chd8$gene_id)

# output labels, correlation, number of genes
labels.out <- centers.label %>% select(-n) %>% ungroup() %>%
  mutate(mouse_metagene = case_when(top.m ~ mouse_metagene),
         r = case_when(top.m ~ r)) %>%
  left_join(count(metagenes.out, h.label, human_metagene),
            by = c('h.label', 'human_metagene')) %>%
  transmute(human_metagene, original_label = h.label,
            n_genes = n, mouse_metagene, cor = r)

# write CSV
write_csv(metagenes.out, 'Polioudakis_metagenes.csv')
write_csv(labels.out, 'Polioudakis_metagene_labels.csv')

# save variables
save(sax.cl, hc, clusters.tb, centers.new, file = 'Polioudakis_metagene_cluster_data.Rdata')