# location-scale meta-analysis for Peng et al 2019
# need to execute init-dir.R first
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# data from local directory
dat <- read_csv(paste0(wkdir, 'native-exotic-richness-relationships/data/DataS2_MetaAnalysisDatabase.csv'))

# start by fitting the model used by Peng et al.
peng_m1 <- brm(z | se(sqrt(var_z)) ~ ln_Grain + (1 | Study / Case),
          data = dat,
          prior = c(prior(normal(0.3,1), class = Intercept),
                    prior(normal(0,1), class = b),
                    prior(normal(0,1), class = sd)),
          cores = 4, chains = 4,
          backend = 'cmdstanr')

# kfold cross validation (in parallel on my local machine)
plan(multisession, workers = 8)
peng_m1_kfold <- kfold(peng_m1,
                       folds = 'stratified',
                       group = 'Study', k = 10)

save(peng_m1, peng_m1_kfold,
     file = paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m1.Rdata'))


# model 1.2: sigma (i.e., scale parameter) as a function of grain, 
peng_sigma_linear <- brm(bf(z | se(sqrt(var_z), sigma = TRUE) ~ ln_Grain + (1 | Study),
                            sigma ~ ln_Grain),
                         prior = c(prior(normal(0.3,1), class = Intercept),
                                   prior(normal(0,1), class = b),
                                   prior(normal(0,1), class = Intercept, dpar = sigma),
                                   prior(normal(0,1), class = b, dpar = sigma),
                                   prior(normal(0,1), class = sd)),
                         data = dat,
                         cores = 4, chains = 4,
                         control = list(adapt_delta = 0.99),
                         backend = 'cmdstanr',
                         save_pars = save_pars(all = TRUE))

# kfold cross validation (in parallel on my local machine)
plan(multisession, workers = 8)
peng_sigma_linear_kfold <- kfold(peng_sigma_linear,
                                 folds = 'stratified',
                                 group = 'Study', k = 10)

save(peng_sigma_linear, peng_sigma_linear_kfold,
     file = paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m2.Rdata'))

# model 1.3: sigma ~ f(extent)
peng_sigma_extent <- brm(bf(z | se(sqrt(var_z), sigma = TRUE) ~ ln_Grain + 
                              (1 | Study),
                            sigma ~ 0 + Extent),
                         prior = c(prior(normal(0.3,1), class = Intercept),
                                   prior(normal(0,1), class = b),
                                   prior(normal(0,1), class = b, dpar = sigma),
                                   prior(normal(0,1), class = sd)),
                         data = dat,
                         cores = 4, chains = 4,
                         backend = 'cmdstanr',
                         control = list(adapt_delta = 0.99))

# kfold cross validation (in parallel on my local machine)
plan(multisession, workers = 8)
peng_sigma_extent_kfold <- kfold(peng_sigma_extent,
                                 folds = 'stratified',
                                 group = 'Study',
                                 k = 10)

save(peng_sigma_extent, peng_sigma_extent_kfold,
     file = paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m3.Rdata'))


# model 1.4: varying study level residuals
peng_sigma_study <- brm(bf(z | se(sqrt(var_z), sigma = TRUE) ~ ln_Grain + (1 | Study),
             sigma ~ 1 + (1 | Study)),
          prior = c(prior(normal(0.3,1), class = Intercept),
                    prior(normal(0,1), class = b),
                    prior(normal(0,1), class = Intercept, dpar = sigma),
                    prior(normal(0,1), class = sd),
                    prior(normal(0,1), class = sd, dpar = sigma)),
          data = dat,
          cores = 4, chains = 4,
          control = list(adapt_delta = 0.99),
          backend = 'cmdstanr',
          save_pars = save_pars(all = TRUE))

# kfold cross validation (in parallel on my local machine)
plan(multisession, workers = 8)
peng_sigma_study_kfold <- kfold(peng_sigma_study,
                                 folds = 'stratified',
                                 group = 'Study', k = 10)

save(peng_sigma_study, peng_sigma_study_kfold,
     file = paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m4.Rdata'))
