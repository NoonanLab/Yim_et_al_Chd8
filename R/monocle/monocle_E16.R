# load packages
require(monocle3)
require(tidyverse)
require(tictoc)

# load CDS files
message('loading CDS file')
cds_E16_filtered <- readRDS('cds_filtered/cds_filtered_E16.rds')

# run monocle::fit_models
tic('E16.0'); message('run monocle::fit_models')
cds_E16.fits <- cds_E16_filtered %>%
  lapply(fit_models, model_formula_str = '~genotype + sex', verbose = T)
cds_E16.fit_coefs <- cds_E16.fits %>%
  lapply(coefficient_table) %>% 
  lapply(select, -model, -model_summary)
toc()

# save model fits
message('saving')
saveRDS(cds_E16.fit_coefs, 'monocle_fitcoefs_E16.rds')

message('done.')