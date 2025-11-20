# code to fit models to patch-scale effort-standardised species richness as 
# function of habitat fragment size

# need to execute init-dir.R first
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# load the data: see Chase et al 2020 Nature (and associated code repo) for full details
frag <- read_csv(paste0(wkdir, 'fragment-size-richness-relationships/data/2_biodiv_frag_fcont_10_mabund_as_is.csv'))

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

#----- simplest model: diversity as a function of fragment size; 
# allow fragment size (slope) to vary by study (varying intercept)----

# sample effort standardised species richness 
Sstd_lognorm_fragSize <- brm(S_std_mean ~ c.lfs + (c.lfs | dataset_label), 
                             # two observations have zero for the response, remove for lognormal distribution
                             data = frag %>% filter(S_std_mean>0),
                             family = lognormal(), # our standardised richness are not integer values
                             prior = c(prior(normal(2.5,1), class = Intercept),
                                       prior(normal(0,1), class = b),
                                       prior(normal(0,1), class = sd),
                                       prior(normal(0,1), class = sigma)),
                             cores = 4, chains = 4,
                             iter = 4000, thin = 2,
                             backend = 'cmdstanr')

# refit with varying (study-level) residual variation
Sstd_lognorm_fragSize_sigma <- brm(bf(S_std_mean ~ c.lfs + (c.lfs | dataset_label),
                                      sigma ~ 1 + (1 | dataset_label)),
                                   # two observations have zero for the response, remove for lognormal distribution
                                   data = frag %>% filter(S_std_mean>0),
                                   family = lognormal(), # our standardised richness are not integer values
                                   prior = c(prior(normal(2.5,1), class = Intercept),
                                             prior(normal(0,1), class = b),
                                             prior(normal(0,1), class = sd),
                                             prior(normal(0,1), class = sd, dpar = sigma),
                                             prior(normal(0,1), class = Intercept, dpar = sigma)),
                                   cores = 4, chains = 4,
                                   iter = 4000, thin = 2,
                                   backend = 'cmdstanr')

save(Sstd_lognorm_fragSize, Sstd_lognorm_fragSize_sigma,
     file = paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/Sstd-m1-m2.Rdata'))

# refit with varying residual variation allowed to change as (log) linear
# function of fragment size
Sstd_lognorm_fragSize_sigma_fs <- brm(bf(S_std_mean ~ c.lfs + (c.lfs | dataset_label),
                                         sigma ~ c.lfs + (c.lfs | dataset_label)),
                                      # two observations have zero for the response, remove for lognormal distribution
                                      data = frag %>% filter(S_std_mean>0),
                                      family = lognormal(), # our standardised richness are not integer values
                                      prior = c(prior(normal(2.5,1), class = Intercept),
                                                prior(normal(0,1), class = b),
                                                prior(normal(0,1), class = sd),
                                                prior(normal(0,1), class = sd, dpar = sigma),
                                                prior(normal(0,1), class = Intercept, dpar = sigma),
                                                prior(normal(0,1), class = b, dpar = sigma)),
                                      cores = 4, chains = 4,
                                      iter = 4000, thin = 2,
                                      control = list(adapt_delta = 0.952),
                                      backend = 'cmdstanr')

save(Sstd_lognorm_fragSize_sigma_fs,
     file = paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/Sstd-m3.Rdata'))

# allow study-level parameters for location and scale to covary
Sstd_lognorm_fragSize_sigma_cor <- brm(bf(S_std_mean ~ c.lfs + (c.lfs | p | dataset_label),
                                          sigma ~ 1 + (1 | p | dataset_label)),
                                       # two observations have zero for the response, remove for lognormal distribution
                                       data = frag %>% filter(S_std_mean>0),
                                       family = lognormal(), # our standardised richness are not integer values
                                       prior = c(prior(normal(2.5,1), class = Intercept),
                                                 prior(normal(0,1), class = b),
                                                 prior(normal(0,1), class = sd),
                                                 prior(normal(0,1), class = sd, dpar = sigma),
                                                 prior(normal(0,1), class = Intercept, dpar = sigma)),
                                       cores = 4, chains = 4,
                                       iter = 4000, thin = 2,
                                       backend = 'cmdstanr')

save(Sstd_lognorm_fragSize_sigma_cor,
     file = paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/Sstd-m4.Rdata'))

# varying residuals as a function of fragment size,
# with correlation among varying effects for location and scale parameters
Sstd_lognorm_fragSize_sigma_fs_cor <- brm(bf(S_std_mean ~ c.lfs + (c.lfs | p | dataset_label),
                                             sigma ~ c.lfs + (c.lfs | p | dataset_label)),
                                          # two observations have zero for the response, remove for lognormal distribution
                                          data = frag %>% filter(S_std_mean>0),
                                          family = lognormal(), # our standardised richness are not integer values
                                          prior = c(prior(normal(2.5,1), class = Intercept),
                                                    prior(normal(0,1), class = b),
                                                    prior(normal(0,1), class = sd),
                                                    prior(normal(0,1), class = sd, dpar = sigma),
                                                    prior(normal(-1,1), class = Intercept, dpar = sigma),
                                                    prior(normal(0,1), class = b, dpar = sigma)),
                                          cores = 8, chains = 8,
                                          iter = 2000, thin = 2,
                                          control = list(adapt_delta = 0.9),
                                          backend = 'cmdstanr')

save(Sstd_lognorm_fragSize_sigma_fs_cor,
     file = paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/Sstd-m5.Rdata'))
