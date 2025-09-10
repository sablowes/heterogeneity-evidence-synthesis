# sbc to test location-scale models

# set up for sbc
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# want to base simulations on data similar to those in Peng et al. meta-analysis
dat <- read_csv(paste0(wkdir, '/native-exotic-richness-relationships/data/DataS2_MetaAnalysisDatabase.csv'))

# use continuous predictor, want it to look like ln_Grain in the Peng et al data
# these are the sample sizes in Peng et al. (spoiler: things don't look good)
peng_template_data1 <- tibble(z = dat$z,
                              var_z = dat$var_z,
                              x = dat$ln_Grain,
                              Study = dat$Study,
                              Case = dat$Case)

# priors are needed to fully define the data generating process
peng_priors1 <- prior(normal(0.13,0.06), class = "Intercept") + # prior for mu (intercept)
  prior(normal(0.04,0.01), class = 'b') + # prior for mu (slope)
  prior(normal(0,0.25), class = sd, group = Study, coef = Intercept) +
  prior(normal(0,0.2), class = sd, group = Study:Case, coef = Intercept)

# use brms to simulate data (with same model as Peng et al. meta-regression)
peng_generator1 <- SBC::SBC_generator_brms(z | se(var_z) ~ x + (1 | Study/Case), 
                                           data = peng_template_data1, 
                                           prior = peng_priors1, 
                                           thin = 100, warmup = 15000, 
                                           refresh = 2000,
                                           control = list(adapt_delta = 0.9),
                                           backend = 'cmdstanr')

# generate enough to pick up problems if they exist
# NB: fits to simulated data can have problems (e.g., due to the combination of 
# simulated effect size and known se, so some rhats > 1.05, and potentially many
# divergent transitions. So, we'll end up discarding some of these
# Vignettes show how rejection sampling can be used to prevent unrealistic data being included
# in these datasets (e.g., many abs(z) > 3), though we'd need to consider the 
# combination of simulated effect size and the known se it gets assigned here. 
# Instead, I take a brute force approach (i.e., many simulated data sets), and 
# will discard poor fits to the simulated data
peng_datasets1 <- generate_datasets(peng_generator1, 500)

# backend required for sbc: these are the models fit to the simulated data 
peng_backend1 <- SBC_backend_brms_from_generator(peng_generator1, 
                                                 chains = 2, thin = 1,
                                                 prior = c(prior(normal(0.3,1), class = Intercept),
                                                           prior(normal(0,1), class = b),
                                                           prior(normal(0,1), class = sd)),
                                                 warmup = 1000, iter = 2000,
                                                 control = list(adapt_delta = 0.95),
                                                 backend = 'cmdstanr')
# do sbc
peng_results1 <- compute_SBC(peng_datasets1, peng_backend1)

# some bad diagnostics
peng_ok1 <- peng_results1$backend_diagnostics %>% 
  filter(n_divergent == 0) %>% 
  pull(sim_id)

# all fits have NAs in the Rhats (due to the fixed sigma parameter)
# want to find fits with rhats > 1.05
max_rhat1 <- c()
for(i in 1:length(peng_results1$fits[peng_ok1])){
  print(i, ' of ', length(peng_ok1))
  max_rhat1[i] <- max(rhat(peng_results1$fits[peng_ok1][[i]]), na.rm = TRUE)
}

hist(max_rhat1)

peng_ok1 <- peng_ok1[max_rhat1 < 1.1]
length(peng_ok1)


# tidy up the labels for the plots (facets)
labels <- as_labeller(c("Intercept" = "beta[0]",
            "b_x" = "beta[1]",
            "sd_Study__Intercept" = "tau",
            "sd_Study:Case__Intercept" = "omega"), default = label_parsed)

# some visual diagnostics
pdf(paste0(wkdir, 'native-exotic-richness-relationships/figures/peng_sbc1.pdf'),
    width = 7, height = 5)

plot_rank_hist(peng_results1[peng_ok1],
               variables = c("Intercept", "b_x", 
                             "sd_Study__Intercept",
                             "sd_Study:Case__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "sd_Study:Case__Intercept")),
             labeller = labels)

plot_coverage(peng_results1[peng_ok1],
              variables = c("Intercept", "b_x", 
                            "sd_Study__Intercept",
                            "sd_Study:Case__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "sd_Study:Case__Intercept")),
             labeller = labels)

plot_ecdf(peng_results1[peng_ok1],
               variables = c("Intercept", "b_x", 
                             "sd_Study__Intercept",
                             "sd_Study:Case__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "sd_Study:Case__Intercept")),
             labeller = labels)

plot_sim_estimated(peng_results1[peng_ok1],
                   variables = c("Intercept", "b_x", 
                                 "sd_Study__Intercept",
                                 "sd_Study:Case__Intercept")) +
  facet_wrap(~factor(variable, 
                     c("Intercept", "b_x", 
                       "sd_Study__Intercept",
                       "sd_Study:Case__Intercept")),
             labeller = labels, scales = 'free')

dev.off()

# with the bad fits this is ~ 20Gb, but < 0.5Gb for these fits (i.e.,
# no divergent transitions, rhats < 1.1)
print(object.size(peng_results1[peng_ok1]), units = "Mb")

peng_sbc_results_1 <- peng_results1[peng_ok1]

save(peng_sbc_results_1,
     file = paste0(wkdir, '../model_fits/sbc-results/peng_sbc_results_goodfits_1.1-extra.Rdata'))
      
     