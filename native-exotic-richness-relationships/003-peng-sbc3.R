# sbc to test location-scale models

# set up for sbc
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# want to base simulations on data similar to those in Peng et al. meta-analysis
dat <- read_csv(paste0(wkdir, '/native-exotic-richness-relationships/data/DataS2_MetaAnalysisDatabase.csv'))


# use continuous predictor, want it to look like ln_Grain in the Peng et al data
# these are the sample sizes in Peng et al. (spoiler: things don't look good)
peng_template_data3a <- tibble(z = dat$z,
                              var_z = dat$var_z,
                              x = dat$ln_Grain,
                              Extent = dat$Extent,
                              Study = dat$Study)

# priors are needed to fully define the data generating process
peng_priors3a <- prior(normal(0.13,0.05), class = "Intercept") + # mu
  prior(normal(0.04,0.02), class = 'b') +
  prior(normal(-1,1), class = b, dpar = sigma) + 
  prior(normal(0,0.2), class = sd, group = Study, lb = 0)

# use brms to simulate data (with same model as Peng et al. meta-regression)
peng_generator3a <- SBC::SBC_generator_brms(bf(z | se(sqrt(var_z), sigma = TRUE) ~ x + 
                                                 (1 | Study),
                                               sigma ~ 0 + Extent),
                                           data = peng_template_data3a,
                                           prior = peng_priors3a, 
                                           thin = 100, warmup = 10000, 
                                           refresh = 2000)

# generate enough datasets to discover problems if they exist
peng_datasets3a <- generate_datasets(peng_generator3a, 200)

# backend required for sbc: these are the models fit to the simulated data 
peng_backend3a <- SBC_backend_brms_from_generator(peng_generator3a, 
                                                 chains = 4, thin = 1,
                                                 prior = c(prior(normal(0.3,1), class = Intercept),
                                                           prior(normal(0,1), class = b),
                                                           prior(normal(0,1), class = b, dpar = sigma),
                                                           prior(normal(0,1), class = sd)),
                                                 control = list(adapt_delta = 0.99),
                                                 warmup = 1000, iter = 2000)
# do sbc
peng_results3a <- compute_SBC(peng_datasets3a, peng_backend3a)

# keep fits with no divergent transitions
peng_ok3a <- peng_results3a$backend_diagnostics %>% 
  filter(n_divergent<10) %>% 
  pull(sim_id)

# all fits have NAs in the Rhats (due to the fixed sigma parameter)
# want to find fits with rhats > 1.05
max_rhat3 <- c()
for(i in 1:length(peng_results3a$fits[peng_ok3a])){
  print(i)
  max_rhat3[i] <- max(rhat(peng_results3a$fits[peng_ok3a][[i]]), na.rm = TRUE)
}

hist(max_rhat3)

peng_ok3a <- peng_ok3a[max_rhat3 < 1.05]
length(peng_ok3a)

# some visual diagnostics
# peng_results3a$fits[[1]] %>% variables
labels <- as_labeller(c("Intercept" = "beta[0]",
                        "b_x" = "beta[1]",
                        "sd_Study__Intercept" = "tau",
                        "b_sigma_Extent0M10" = "beta[0-10]^sigma",
                        "b_sigma_Extent10M100" = "beta[10-100]^sigma",
                        "b_sigma_Extent102M103" = "beta[10^2-10^3]^sigma",
                        "b_sigma_Extent103M104" = "beta[10^3-10^4]^sigma",
                        "b_sigma_Extent104M105" = "beta[10^4-10^5]^sigma",
                        "b_sigma_Extent105M106" = "beta[10^5-10^6]^sigma",
                        "b_sigma_Extent106M" = "beta[10^6]^sigma"), 
                      default = label_parsed)

pdf(paste0(wkdir, 'native-exotic-richness-relationships/figures/peng_sbc3.pdf'),
    width = 7, height = 5)

plot_rank_hist(peng_results3a[peng_ok3a],
               variables = c("Intercept", "b_x", 
                             "sd_Study__Intercept",
                             "b_sigma_Extent0M10",
                             "b_sigma_Extent10M100",
                             "b_sigma_Extent102M103",
                             "b_sigma_Extent103M104",
                             "b_sigma_Extent104M105",
                             "b_sigma_Extent105M106",
                             "b_sigma_Extent106M")) +
  facet_wrap(~factor(variable,
                     levels = c("Intercept", "b_x", 
                                "sd_Study__Intercept",
                                "b_sigma_Extent0M10",
                                "b_sigma_Extent10M100",
                                "b_sigma_Extent102M103",
                                "b_sigma_Extent103M104",
                                "b_sigma_Extent104M105",
                                "b_sigma_Extent105M106",
                                "b_sigma_Extent106M")),
                     labeller = labels)

plot_coverage_diff(peng_results3a[peng_ok3a],
              variables = c("Intercept", "b_x", 
                            "b_sigma_Extent0M10",
                            "b_sigma_Extent10M100",
                            "b_sigma_Extent102M103",
                            "b_sigma_Extent103M104",
                            "b_sigma_Extent104M105",
                            "b_sigma_Extent105M106",
                            "b_sigma_Extent106M",
                            "sd_Study__Intercept"))+
  facet_wrap(~factor(variable,
                     levels = c("Intercept", "b_x", 
                                "sd_Study__Intercept",
                                "b_sigma_Extent0M10",
                                "b_sigma_Extent10M100",
                                "b_sigma_Extent102M103",
                                "b_sigma_Extent103M104",
                                "b_sigma_Extent104M105",
                                "b_sigma_Extent105M106",
                                "b_sigma_Extent106M")),
             labeller = labels)

plot_ecdf(peng_results3a[peng_ok3a],
          variables = c("Intercept", "b_x", 
                        "b_sigma_Extent0M10",
                        "b_sigma_Extent10M100",
                        "b_sigma_Extent102M103",
                        "b_sigma_Extent103M104",
                        "b_sigma_Extent104M105",
                        "b_sigma_Extent105M106",
                        "b_sigma_Extent106M",
                        "sd_Study__Intercept"))+
  facet_wrap(~factor(variable,
                     levels = c("Intercept", "b_x", 
                                "sd_Study__Intercept",
                                "b_sigma_Extent0M10",
                                "b_sigma_Extent10M100",
                                "b_sigma_Extent102M103",
                                "b_sigma_Extent103M104",
                                "b_sigma_Extent104M105",
                                "b_sigma_Extent105M106",
                                "b_sigma_Extent106M")),
             labeller = labels)

plot_sim_estimated(peng_results3a[peng_ok3a],
                   variables = c("Intercept", "b_x", 
                                 "b_sigma_Extent0M10",
                                 "b_sigma_Extent10M100",
                                 "b_sigma_Extent102M103",
                                 "b_sigma_Extent103M104",
                                 "b_sigma_Extent104M105",
                                 "b_sigma_Extent105M106",
                                 "b_sigma_Extent106M",
                                 "sd_Study__Intercept")) +
  facet_wrap(~factor(variable,
                     levels = c("Intercept", "b_x", 
                                "sd_Study__Intercept",
                                "b_sigma_Extent0M10",
                                "b_sigma_Extent10M100",
                                "b_sigma_Extent102M103",
                                "b_sigma_Extent103M104",
                                "b_sigma_Extent104M105",
                                "b_sigma_Extent105M106",
                                "b_sigma_Extent106M")),
             labeller = labels, scales = 'free')

dev.off()

print(object.size(peng_results3a[peng_ok3a]), units = "Mb")

peng_sbc_results_3a <- peng_results3a[peng_ok3a]

save(peng_sbc_results_3a,
     file = paste0(wkdir, '../model_fits/sbc-results/peng_sbc_results_goodfits_3a.Rdata'))
