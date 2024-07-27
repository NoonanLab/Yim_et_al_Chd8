# load packages
require(monocle3)
require(tidyverse)
require(tictoc)

# load CDS files
message('loading CDS file')
cds_E17_filtered <- readRDS('cds_filtered/cds_filtered_E17.rds')

# run monocle::fit_models 
tic('E17.5'); message('run monocle::fit_models')
cds_E17.fits <- cds_E17_filtered %>%
  lapply(fit_models, model_formula_str = '~genotype + sex + batch', verbose = T)
cds_E17.fit_coefs <- cds_E17.fits %>%
  lapply(coefficient_table) %>% 
  lapply(select, -model, -model_summary)
toc()

# save model fits
message('saving')
saveRDS(cds_E17.fit_coefs, 'monocle_fitcoefs_E17.rds')

message('done.')