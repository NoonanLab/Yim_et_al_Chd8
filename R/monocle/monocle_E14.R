# load packages
require(monocle3)
require(tidyverse)
require(tictoc)

# load CDS files
message('loading CDS file')
cds_E14_filtered <- readRDS('cds_filtered/cds_filtered_E14.rds')

# run monocle::fit_models
tic('E14.5'); message('run monocle::fit_models')
cds_E14.fits <- cds_E14_filtered %>%
  lapply(fit_models, model_formula_str = '~genotype + sex + batch', verbose = T)
cds_E14.fit_coefs <- cds_E14.fits %>%
  lapply(coefficient_table) %>% 
  lapply(select, -model, -model_summary)
toc()

# save model fits
message('saving')
saveRDS(cds_E14.fit_coefs, 'monocle_fitcoefs_E14.rds')

message('done.')