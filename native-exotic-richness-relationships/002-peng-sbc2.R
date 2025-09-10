# sbc to test location-scale models

# set up for sbc
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# want to base simulations on data similar to those in Peng et al. meta-analysis
dat <- read_csv(paste0(wkdir, '/native-exotic-richness-relationships/data/DataS2_MetaAnalysisDatabase.csv'))


# use continuous predictor, want it to look like ln_Grain in the Peng et al data
# these are the sample sizes in Peng et al. (spoiler: things don't look good)
peng_template_data2 <- tibble(z = dat$z,
                              var_z = dat$var_z,
                              x = dat$ln_Grain,
                              Study = dat$Study,
                              Case = dat$Case)

# priors are needed to fully define the data generating process
peng_priors2 <- prior(normal(0.13,0.06), class = "Intercept") + # mu
  prior(normal(0.04,0.02), class = 'b') +
  prior(normal(-1.1, 0.2), class = Intercept, dpar = sigma) + 
  prior(normal(-0.03, 0.03), class = b, dpar = sigma) + 
  prior(normal(0,0.22), class = sd, group = Study, lb = 0)

# use brms to simulate data (with same model as Peng et al. meta-regression)
peng_generator2 <- SBC::SBC_generator_brms(bf(z | se(sqrt(var_z), 
                                                     sigma = TRUE) ~ x + 
                                                (1 | Study),
                                              sigma ~ x),
                                           data = peng_template_data2,
                                           prior = peng_priors2,
                                           thin = 200, warmup = 15000, 
                                           refresh = 20000)

# generate enough datasets to discover problems if they exist
peng_datasets2 <- generate_datasets(peng_generator2, 200)

# backend required for sbc: these are the models fit to the simulated data 
peng_backend2 <- SBC_backend_brms_from_generator(peng_generator2, 
                                                 chains = 8, 
                                                 prior = c(prior(normal(0.3,1), class = Intercept),
                                                           prior(normal(0,1), class = b),
                                                           prior(normal(0,1), class = Intercept, dpar = sigma),
                                                           prior(normal(0,1), class = b, dpar = sigma),
                                                           prior(normal(0,1), class = sd)),
                                                 warmup = 1000, iter = 4000,
                                                 control = list(adapt_delta = 0.99))
# do sbc
peng_results2 <- compute_SBC(peng_datasets2, 
                             peng_backend2)

# some bad diagnostics (well lots actually)
peng_ok2 <- peng_results2$backend_diagnostics %>% 
  filter(n_divergent==0) %>% 
  pull(sim_id)

# all fits have NAs in the Rhats (due to the fixed sigma parameter)
# want to find fits with rhats > 1.05
max_rhat2 <- c()
for(i in 1:length(peng_results2$fits[peng_ok2])){
  print(i)
  max_rhat2[i] <- max(rhat(peng_results2$fits[peng_ok2][[i]]), na.rm = TRUE)
}

hist(max_rhat2)

peng_ok2 <- peng_ok2[max_rhat2 < 1.15]
length(peng_ok2)

# some visual diagnostics
# tidy up the facet labels
labels <- as_labeller(c("Intercept" = "beta[0]",
                        "b_x" = "beta[1]",
                        "sd_Study__Intercept" = "tau",
                        "b_sigma_Intercept" = "beta[0]^sigma",
                        "b_sigma_x" = "beta[1]^sigma"), 
                      default = label_parsed)

pdf(paste0(wkdir, 'native-exotic-richness-relationships/figures/peng_sbc2.pdf'),
    width = 7, height = 5)

plot_rank_hist(peng_results2[peng_ok2],
               variables = c("Intercept", "b_x", 
                             "b_sigma_Intercept",
                             "b_sigma_x",
                             "sd_Study__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "b_sigma_x")),
             labeller = labels)

plot_coverage(peng_results2[peng_ok2],
              variables = c("Intercept", "b_x", 
                            "b_sigma_Intercept",
                            "b_sigma_x",
                            "sd_Study__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "b_sigma_x")),
             labeller = labels)

plot_coverage_diff(peng_results2[peng_ok2],
              variables = c("Intercept", "b_x", 
                            "b_sigma_Intercept",
                            "b_sigma_x",
                            "sd_Study__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "b_sigma_x")),
             labeller = labels)

plot_ecdf(peng_results2[peng_ok2],
               variables = c("Intercept", "b_x", 
                             "b_sigma_Intercept",
                             "b_sigma_x",
                             "sd_Study__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "b_sigma_x")),
             labeller = labels)

plot_ecdf_diff(peng_results2[peng_ok2],
          variables = c("Intercept", "b_x", 
                        "b_sigma_Intercept",
                        "b_sigma_x",
                        "sd_Study__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "b_sigma_x")),
             labeller = labels)
plot_sim_estimated(peng_results2[peng_ok2],
                   variables = c("Intercept", "b_x", 
                                 "b_sigma_Intercept",
                                 "b_sigma_x",
                                 "sd_Study__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "b_sigma_Intercept",
                       "b_sigma_x")),
             labeller = labels, scales = 'free')

dev.off()

print(object.size(peng_results2[peng_ok2]), units = "Mb")

peng_sbc_results_2 <- peng_results2[peng_ok2]

save(peng_sbc_results_2,
     file = paste0(wkdir, '../model_fits/sbc-results/peng_sbc_results_goodfits_2.Rdata'))
