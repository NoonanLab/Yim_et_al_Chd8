# load packages
require(tidyverse)
require(jmotif)

# load data
bin.data <- readRDS('bin_data_filtered.rds')

# scale data
bin.data.scale <- bin.data %>% 
  unite('rowname', time_point, genotype, mouse_id, mouse_name) %>%
  column_to_rownames() %>% t() %>% scale() %>% t() %>% 
  as_tibble(rownames = 'rowname') %>%
  separate(rowname, into = c('time_point', 'genotype', 'mouse_id', 'mouse_name'), 
           sep = '_') %>%
  mutate(genotype = factor(genotype, levels = levels(bin.data$genotype)))

# SAX variables
paa.size <- 8
alphabet.size <- 9

# compute SAX transform
bin.paa8.sax9 <- bin.data.scale %>% rowwise() %>%
  transmute(time_point, genotype, mouse_id, mouse_name,
            bin.mean.scale = list(c_across(starts_with('bin')))) %>%
  mutate(paa = list(paa(bin.mean.scale, paa_num = paa.size)),
         sax.word = list(series_to_chars(paa, a_size = alphabet.size)),
         sax9 = list(letters_to_idx(sax.word))) %>%
  mutate(sax.word = str_c(sax.word, collapse = '')) %>%
  mutate(na.expression = any(is.na(bin.mean.scale))) %>%
  ungroup()

saveRDS(bin.paa8.sax9, 'bin_paa8_alpha9.rds')