# sbc to test location-scale models for richness ~ f(frag_size)

# set up for sbc
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# sbc needs a data generating process. 

# To use brms, we need a "template dataset", and then we can use brms 
# to define the data generating process.

# The predictor(s) will be used for data generation; 
# the response (standardised richness) values will be ignored, 
# but need to be present and of the correct data type
set.seed(123123)

# want to use empirical predictor data; 
frag <- read_csv('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/fragment-size-richness-relationships/data/2_biodiv_frag_fcont_10_mabund_as_is.csv')

# add mean centred (log) fragsize
frag$c.lfs <- log(frag$frag_size_num) - mean(log(frag$frag_size_num))

frag_template_data <- tibble(S_std_mean = frag$S_std_mean, 
                             c.lfs = frag$c.lfs,
                             dataset_label = frag$dataset_label) %>% 
  # lognormal response distribution, remove zeroes
  filter(S_std_mean > 0)

# model 2.4: same priors as model 2.2, with default prior on correlations
frag_generator4 <- SBC::SBC_generator_brms(bf(S_std_mean ~ c.lfs + 
                                                (c.lfs | p | dataset_label), 
                                              sigma ~ 1 + (1 | p | dataset_label)),
                                           family = lognormal(),
                                           data = frag_template_data, 
                                           prior = frag_prior2, 
                                           thin = 100, warmup = 15000, 
                                           refresh = 2000)

# generate the "fake" data sets: 
frag_datasets4 <- generate_datasets(frag_generator4, 200)

frag_backend4 <- SBC_backend_brms_from_generator(frag_generator4, chains = 2, 
                                                 thin = 1,
                                                 warmup = 1000, iter = 2000,
                                                 prior = c(prior(normal(2.5,1), class = Intercept),
                                                           prior(normal(0,1), class = b),
                                                           prior(normal(0,1), class = sd),
                                                           prior(normal(0,1), class = sd, dpar = sigma),
                                                           prior(normal(-1,1), class = Intercept, dpar = sigma)),
                                                 control = list(adapt_delta = 0.9),
                                                 backend = 'cmdstanr')

frag_results4 <- compute_SBC(frag_datasets4, frag_backend4)

ok_frag4 <- frag_results4$backend_diagnostics %>% 
  filter(n_divergent == 0) %>% 
  pull(sim_id)

# rhat diagnostics in sbc output have NAs (not sure why)
# want to find fits with rhats > 1.05
max_rhat4 <- c()
for(i in 1:length(frag_results4$fits[ok_frag4])){
  print(i)
  max_rhat4[i] <- max(rhat(frag_results4$fits[ok_frag4][[i]]), na.rm = TRUE)
}

# remove the some more problematic fits
ok_frag4 <- ok_frag4[max_rhat4 < 1.05]

frag_sbc_results_4 <- frag_results4[ok_frag4]


# plot diagnostics
label4 <- as_labeller(c("b_Intercept" = "beta[0]",
                        "b_c.lfs" = "beta[1]",
                        "sd_dataset_label__Intercept" = "sigma[0]",
                        "sd_dataset_label__c.lfs" = "sigma[1]",
                        "b_sigma_Intercept" = "beta[0]^sigma",
                        "sd_dataset_label__sigma_Intercept" = "sigma[0]^sigma",
                        "cor_dataset_label__Intercept__c.lfs" = "rho[sigma[0]~sigma[1]]"), 
                      default = label_parsed)



pdf(paste0(wkdir, 'fragment-size-richness-relationships/figures/frag_sbc4.pdf'),
    width = 7, height = 5)

plot_rank_hist(frag_sbc_results_4,
               variables = c("b_c.lfs",#"b_Intercept", 
                             "b_sigma_Intercept",
                             # "sd_dataset_label__Intercept",
                             "sd_dataset_label__c.lfs", 
                             "sd_dataset_label__sigma_Intercept")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "b_sigma_Intercept",
                       "sd_dataset_label__sigma_Intercept")),
             labeller = label4)

plot_coverage_diff(frag_sbc_results_4,
              variables = c("b_c.lfs",#"b_Intercept", 
                            "b_sigma_Intercept",
                            # "sd_dataset_label__Intercept",
                            "sd_dataset_label__c.lfs", 
                            "sd_dataset_label__sigma_Intercept",
                            "cor_dataset_label__Intercept__c.lfs")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "b_sigma_Intercept",
                       "sd_dataset_label__sigma_Intercept",
                       "cor_dataset_label__Intercept__c.lfs")),
             labeller = label4)

plot_ecdf_diff(frag_sbc_results_4,
          variables = c("b_c.lfs",#"b_Intercept", 
                        "b_sigma_Intercept",
                        # "sd_dataset_label__Intercept",
                        "sd_dataset_label__c.lfs", 
                        "sd_dataset_label__sigma_Intercept",
                        "cor_dataset_label__Intercept__c.lfs")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "b_sigma_Intercept",
                       "sd_dataset_label__sigma_Intercept",
                       "cor_dataset_label__Intercept__c.lfs")),
             labeller = label4)

plot_sim_estimated(frag_sbc_results_4,
                   variables = c("b_c.lfs",#"b_Intercept", 
                                 "b_sigma_Intercept",
                                 # "sd_dataset_label__Intercept",
                                 "sd_dataset_label__c.lfs", 
                                 "sd_dataset_label__sigma_Intercept",
                                 "cor_dataset_label__Intercept__c.lfs")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "b_sigma_Intercept",
                       "sd_dataset_label__sigma_Intercept",
                       "cor_dataset_label__Intercept__c.lfs")),
             labeller = label4, scales = 'free')
dev.off()

# save 
print(object.size(frag_sbc_results_4), units = "Gb")

length(frag_sbc_results_4) # 181

save(frag_sbc_results_4,
     file = paste0(wkdir, '../model_fits/sbc-results/frag_sbc_results_goodfits_2.4.Rdata'))

