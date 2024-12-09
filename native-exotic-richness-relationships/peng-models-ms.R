# location-scale meta-analysis for Peng et al 2019
# need to execute init-dir.R first
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# data from local directory
dat <- read_csv(paste0(wkdir, 'native-exotic-richness-relationships/data/DataS2_MetaAnalysisDatabase.csv'))

# cluster directory
# dat <- read_csv('/data/idiv_chase/sablowes/evid-synth/data/DataS2_MetaAnalysisDatabase.csv')

# start by fitting the model used by Peng et al.
peng_m1 <- brm(z | se(sqrt(var_z)) ~ ln_Grain + (1 | Study / Case),
          data = dat,
          prior = c(prior(normal(0.3,1), class = Intercept),
                    prior(normal(0,1), class = b),
                    prior(normal(0,1), class = sd)),
          cores = 4, chains = 4,
          save_pars = save_pars(all = TRUE))

peng_m1_kfold <- kfold(peng_m1,
                       folds = 'stratified',
                       group = 'Study', k = 10)

save(peng_m1, peng_m1_kfold,
     file = paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m1.Rdata'))

# model 1.2: sigma (i.e., parameter associated with scale) as a function of grain, 
peng_sigma_linear <- brm(bf(z | se(sqrt(var_z), sigma = TRUE) ~ ln_Grain + (1 | Study),
                            sigma ~ ln_Grain),
                         prior = c(prior(normal(0.3,1), class = Intercept),
                                   prior(normal(0,1), class = b),
                                   prior(normal(0,1), class = Intercept, dpar = sigma),
                                   prior(normal(0,1), class = b, dpar = sigma),
                                   prior(normal(0,1), class = sd)),
                         data = dat,
                         cores = 4, chains = 4,
                         control = list(adapt_delta = 0.95),
                         save_pars = save_pars(all = TRUE))

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

plan(multisession, workers = 8)
peng_sigma_extent_kfold <- kfold(peng_sigma_extent,
                                 folds = 'stratified',
                                 group = 'Study',
                                 k = 10)

save(peng_sigma_extent, peng_sigma_extent_kfold,
     file = paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m3.Rdata'))


# model 1.4: varying study level residuals (via the scale parameter, sigma)
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
          save_pars = save_pars(all = TRUE))

plan(multisession, workers = 8)
peng_sigma_study_kfold <- kfold(peng_sigma_study,
                                 folds = 'stratified',
                                 group = 'Study', k = 10)

save(peng_sigma_study, peng_sigma_study_kfold,
     file = paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m4.Rdata'))


# grain size as a predictor of case-level varying 
# intercepts (omega)
# adapted from: 
# https://discourse.mc-stan.org/t/brms-hacking-linear-predictors-for-random-effect-standard-deviations/34162?u=martinmodrak
# this is akin to existing location-scale models for meta-analysis
# (e.g., Williams et al. 2021, Viechtbauer & Lopez-Lopez 2022, Rodriguez et al. 2023)
peng_sd_linear <- brm(bf(z | se(sqrt(var_z)) ~ muy + vint1 + exp(logphi) * vint2,
                         muy ~ ln_Grain,
                         vint1 ~ 0 + (1 | Study),
                         logphi ~ ln_Grain,
                         vint2 ~ 0 + (1 | Study:Case),
                         nl = TRUE),
                      data = dat,
                      prior = c(prior(normal(0,1), class = b, nlpar = logphi),
                                prior(normal(0.3,1), class = b, coef = Intercept, nlpar = muy),
                                prior(normal(0,1), class = b, coef = ln_Grain, nlpar = muy),
                                prior(normal(0,1), class = sd, nlpar = vint1),
                                prior(constant(1), class = sd, nlpar = vint2)),
                      cores = 4, chains = 4)

peng_sd_linear_kfold <- kfold(peng_sd_extent,
                              folds = 'stratified',
                              group = 'Study',
                              k = 10)

save(peng_sd_linear, peng_sd_linear_kfold,
     file = paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m5-sd-linear.Rdata'))
