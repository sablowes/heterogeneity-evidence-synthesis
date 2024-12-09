# leave-one-group-out (logo) for location-scale models
# fit to peng meta-analysis

## NB: this script was run on a scientific computing cluster

library(tidyverse)
library(brms)


# load model fits
load('/data/idiv_chase/sablowes/evid-synth/results/peng-meta-models.Rdata')

# leave one group out cv
cv10g_m1 <- kfold(m1, group = 'Study')
cv10g_sd_linear <- kfold(peng_sd_linear, group = 'Study')
cv10g_sigma_extent <- kfold(peng_sigma_extent, group = 'Study')
cv10g_sigma_linear <- kfold(peng_sigma_linear, group = 'Study')
cv10g_sigma_study <- kfold(peng_sigma_study, group = 'Study')


save(cv10g_m1, 
     cv10g_sd_linear, cv10g_sigma_extent,
     cv10g_sigma_linear, cv10g_sigma_study,
     file = Sys.getenv('OFILE'))

save(cv10g_sigma_extent, file = 'native-exotic-richness-relationships/model-fits-CV-results/peng-logo-sigma-extent.Rdata')
