# leave-one-group-out (logo) for location-scale models
# need to execute init-dir.R first
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# data from local directory
dat <- read_csv(paste0(wkdir, 'native-exotic-richness-relationships/data/DataS2_MetaAnalysisDatabase.csv'))

# start with the model used by Peng et al.

# load model fit, do logo and save
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m1.Rdata'))
# leave one group out cv
plan(multisession, workers = 8)
cv10g_m1 <- kfold(peng_m1, group = 'Study')
save(cv10g_m1, 
     file = paste0(wkdir, '/native-exotic-richness-relationships/model-fits-CV-results/peng-m1-logo.Rdata'))

# model two
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m2.Rdata'))
plan(multisession, workers = 8)
cv10g_m2 <- kfold(peng_sigma_linear, group = 'Study')
save(cv10g_m2, 
     file = paste0(wkdir, '/native-exotic-richness-relationships/model-fits-CV-results/peng-m2-logo.Rdata'))

# model three
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m3.Rdata'))
plan(multisession, workers = 8)
cv10g_m3 <- kfold(peng_sigma_extent, group = 'Study')
save(cv10g_m3, 
     file = paste0(wkdir, '/native-exotic-richness-relationships/model-fits-CV-results/peng-m3-logo.Rdata'))

# model four
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m4.Rdata'))
plan(multisession, workers = 8)
cv10g_m4 <- kfold(peng_sigma_study, group = 'Study')
save(cv10g_m4, 
     file = paste0(wkdir, '/native-exotic-richness-relationships/model-fits-CV-results/peng-m4-logo.Rdata'))

# model five
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m5.Rdata'))
plan(multisession, workers = 8)
cv10g_m5 <- kfold(peng_sigma_study_corr, group = 'Study')
save(cv10g_m5, 
     file = paste0(wkdir, '/native-exotic-richness-relationships/model-fits-CV-results/peng-m5-logo.Rdata'))


# alternate parameterisation of model 2
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m2-alt-sd-linear.Rdata'))
# leave one group out cv
plan(multisession, workers = 8)
cv10g_m2_alt <- kfold(peng_sd_linear, group = 'Study')
save(cv10g_m2_alt, 
     file = paste0(wkdir, '/native-exotic-richness-relationships/model-fits-CV-results/peng-m2-alt-sd-linear-logo.Rdata'))