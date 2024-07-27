# load packages
require(tidyverse)
require(Seurat)
require(biomaRt)

# create gene table
cc.genes.tb <- tibble(phase = names(cc.genes.updated.2019),
                      human_name = cc.genes.updated.2019) %>% 
  unnest(human_name) %>%
  mutate(phase = str_to_upper(str_remove(phase, '.genes')))

# retrieve Ensembl IDs and mouse orthologs (dec2021)
mart.out.1 <- getBM(attributes = c('external_gene_name',
                                   'ensembl_gene_id',
                                   'mmusculus_homolog_ensembl_gene',
                                   'mmusculus_homolog_associated_gene_name',
                                   'mmusculus_homolog_orthology_type'), 
                    filters = 'hgnc_symbol', values = cc.genes.tb$human_name,
                    mart = useMart(biomart = 'ensembl', 
                                   dataset = 'hsapiens_gene_ensembl', 
                                   host = 'https://dec2021.archive.ensembl.org'))

# get missing orthologs (oct2014)
get.genes <- mart.out.1 %>% add_count(external_gene_name) %>%
  filter(mmusculus_homolog_ensembl_gene == '', n == 1)
mart.out.2 <- getLDS(attributes = c('external_gene_name',
                                    'ensembl_gene_id',
                                    'mmusculus_homolog_ensembl_gene',
                                    'mmusculus_homolog_orthology_type'),
                     filters = 'ensembl_gene_id', 
                     values = get.genes$ensembl_gene_id,
                     mart = useMart(biomart = 'ensembl', 
                                    dataset = 'hsapiens_gene_ensembl', 
                                    host = 'https://oct2014.archive.ensembl.org'),
                     attributesL = 'external_gene_name',
                     martL = useMart(biomart = 'ensembl', 
                                     dataset = 'mmusculus_gene_ensembl', 
                                     host = 'https://oct2014.archive.ensembl.org'))

# combine tables
cc.genes.tb.out <- bind_rows(
  set_names(mart.out.2, c('human_name', 'human_id', 'mouse_id', 'homology_type', 'mouse_name')),
  set_names(mart.out.1, c('human_id', 'human_name', 'mouse_id', 'mouse_name', 'homology_type'))
) %>%
  right_join(cc.genes.tb, by = 'human_name') %>% 
  relocate(phase, human_name, human_id, mouse_name) %>% 
  arrange(phase, human_name) %>%
  filter(mouse_name != '')

# write CSV
write_csv(cc.genes.tb.out, 'seurat_cc_genes_mouse.csv')