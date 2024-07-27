# load packages
require(monocle3)
require(tidyverse)
require(tictoc)

# load CDS files
message('loading CDS file')
cds_P25_filtered <- readRDS('cds_filtered/cds_filtered_P25.rds')

# run monocle::fit_models
tic('P25'); message('run monocle::fit_models')
cds_P25.fits <- cds_P25_filtered %>%
  lapply(fit_models, model_formula_str = '~genotype + sex + batch', verbose = T)
cds_P25.fit_coefs <- cds_P25.fits %>%
  lapply(coefficient_table) %>% 
  lapply(select, -model, -model_summary)
toc()

# save model fits
message('saving')
saveRDS(cds_P25.fit_coefs, 'monocle_fitcoefs_P25.rds')

message('done.')
