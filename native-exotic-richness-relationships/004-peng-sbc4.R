# sbc to test location-scale models

# set up for sbc
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# want to base simulations on data similar to those in Peng et al. meta-analysis
dat <- read_csv(paste0(wkdir, '/native-exotic-richness-relationships/data/DataS2_MetaAnalysisDatabase.csv'))

# use continuous predictor, want it to look like ln_Grain in the Peng et al data
# and start with the sample sizes in Peng et al.
peng_template_data4 <- tibble(z = dat$z,
                              var_z = dat$var_z,
                              x = dat$ln_Grain,
                              Study = dat$Study)

# priors are needed to fully define the data generating process
peng_priors4 <- prior(normal(0.17,0.05), class = Intercept) + # mu
  prior(normal(0.04,0.02), class = b) +
  prior(normal(-1.6, 0.25), class = Intercept, dpar = sigma) + 
  prior(normal(0,0.2), class = sd, group = Study) + 
  prior(normal(0,0.7), class = sd, dpar = sigma)

# use brms to simulate data (with same model as Peng et al. meta-regression)
peng_generator4 <- SBC::SBC_generator_brms(bf(z | se(sqrt(var_z), sigma = TRUE) ~ 
                                                x + (1 | Study),
                                              sigma ~ 1 + (1 | Study)),
                                           data = peng_template_data4,
                                           prior = peng_priors4, 
                                           control = list(adapt_delta = 0.95),
                                           thin = 200, warmup = 15000, 
                                           refresh = 5000,
                                           backend = 'cmdstanr')

# generate datasets 
peng_datasets4 <- generate_datasets(peng_generator4, 500)

# backend required for sbc: these are the models fit to the simulated data 
peng_backend4 <- SBC_backend_brms_from_generator(peng_generator4, 
                                                 chains = 4, thin = 1,
                                                 warmup = 1000, iter = 2000,
                                                 prior = c(prior(normal(0.3,1), class = Intercept),
                                                           prior(normal(0,1), class = b),
                                                           prior(normal(0,1), class = Intercept, dpar = sigma),
                                                           prior(normal(0,1), class = sd),
                                                           prior(normal(0,1), class = sd, dpar = sigma)),
                                                 control = list(adapt_delta = 0.99),
                                                 backend = 'cmdstanr')
# do sbc
peng_results4 <- compute_SBC(peng_datasets4, peng_backend4)

# remove fits with divergent transitions
peng_ok4 <- peng_results4$backend_diagnostics %>% 
  filter(n_divergent == 0) %>% 
  pull(sim_id)

# all fits have NAs in the Rhats (due to how the model is parameterised)
# get the rhats that are not NA for inspection
max_rhat4 <- c()
for(i in 1:length(peng_results4$fits[peng_ok4])){
  print(paste(i, ' in ', length(peng_ok4)))
  max_rhat4[i] <- max(rhat(peng_results4$fits[peng_ok4][[i]]), na.rm = TRUE)
}

hist(max_rhat4)

peng_ok4 <- peng_ok4[max_rhat4 < 1.05]
length(peng_ok4)

# some visual diagnostics
# tidy up the labels for the plots (facets)
labels <- as_labeller(c("Intercept" = "beta[0]",
                        "b_x" = "beta[1]",
                        "sd_Study__Intercept" = "tau",
                        "b_sigma_Intercept" = "beta[0]^sigma",
                        "sd_Study__sigma_Intercept" = "zeta"), 
                      default = label_parsed)

pdf(paste0(wkdir, 'native-exotic-richness-relationships/figures/peng_sbc4.pdf'),
    width = 7, height = 5)

plot_rank_hist(peng_results4[peng_ok4],
               variables = c("Intercept", "b_x", 
                             "b_sigma_Intercept",
                             "sd_Study__Intercept",
                             "sd_Study__sigma_Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "sd_Study__sigma_Intercept")),
             labeller = labels)

plot_coverage(peng_results4[peng_ok4],
              variables = c("Intercept", "b_x", 
                            "b_sigma_Intercept",
                            "sd_Study__Intercept",
                            "sd_Study__sigma_Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "sd_Study__sigma_Intercept")),
             labeller = labels)

plot_ecdf(peng_results4[peng_ok4],
               variables = c("Intercept", "b_x", 
                             "b_sigma_Intercept",
                             "sd_Study__Intercept",
                             "sd_Study__sigma_Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "sd_Study__sigma_Intercept")),
             labeller = labels)

plot_sim_estimated(peng_results4[peng_ok4],
                   variables = c("Intercept", "b_x", 
                                 "b_sigma_Intercept",
                                 "sd_Study__Intercept",
                                 "sd_Study__sigma_Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "sd_Study__sigma_Intercept")),
             labeller = labels, scales = 'free')

dev.off()

# check size before saving
print(object.size(peng_results4[peng_ok4]), 
      units = "Gb")

peng_sbc_results_4 <- peng_results4[peng_ok4]

save(peng_sbc_results_4,
     file = paste0(wkdir, '../model_fits/sbc-results/peng_sbc_results_goodfits_1.4.Rdata'))
