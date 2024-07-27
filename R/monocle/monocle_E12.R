# load packages
require(monocle3)
require(tidyverse)
require(tictoc)

# load CDS files
message('loading CDS file')
cds_E12_filtered <- readRDS('cds_filtered/cds_filtered_E12.rds')

# run monocle::fit_models
tic('E12.5'); message('run monocle::fit_models')
cds_E12.fits <- cds_E12_filtered %>%
  lapply(fit_models, model_formula_str = '~genotype + sex + batch', verbose = T)
cds_E12.fit_coefs <- cds_E12.fits %>%
  lapply(coefficient_table) %>% 
  lapply(select, -model, -model_summary)
toc()

# save model fits
message('saving')
saveRDS(cds_E12.fit_coefs, 'monocle_fitcoefs_E12.rds')

message('done.')