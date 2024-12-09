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

frag_prior2 <- prior(normal(2.5,0.1), class = Intercept) + 
  prior(normal(0.05, 0.01), class = 'b') + 
  prior(normal(1,0.1), class = sd, coef = Intercept, group = dataset_label) +
  prior(normal(0.05,0.05), class = sd, coef = c.lfs, group = dataset_label) +
  prior(normal(-1.2, 0.1), class = Intercept, dpar = sigma) + 
  prior(normal(0.4,0.1), class = sd, dpar = sigma, lb = 0)

frag_generator2 <- SBC::SBC_generator_brms(bf(S_std_mean ~ c.lfs + 
                                                (c.lfs | dataset_label), 
                                              sigma ~ 1 + (1 | dataset_label)),
                                           family = lognormal(),
                                           data = frag_template_data, 
                                           prior = frag_prior2, 
                                           thin = 100, warmup = 15000, 
                                           refresh = 2000)

# generate the "fake" data sets: 
frag_datasets2 <- generate_datasets(frag_generator2, 200)

frag_backend2 <- SBC_backend_brms_from_generator(frag_generator2, chains = 2, 
                                                 thin = 1,
                                                 warmup = 1000, iter = 2000,
                                                 prior = c(prior(normal(2.5,1), class = Intercept),
                                                           prior(normal(0,1), class = b),
                                                           prior(normal(0,1), class = sd),
                                                           prior(normal(0,1), class = sd, dpar = sigma),
                                                           prior(normal(0,1), class = Intercept, dpar = sigma)),
                                                 control = list(adapt_delta = 0.9),
                                                 backend = 'cmdstanr')

frag_results2 <- compute_SBC(frag_datasets2, frag_backend2)

ok_frag2 <- frag_results2$backend_diagnostics %>% 
  filter(n_divergent==0) %>% 
  pull(sim_id)

# rhat diagnostics in sbc output have NAs (not sure why)
# want to find fits with rhats > 1.05
max_rhat2 <- c()
for(i in 1:length(frag_results2$fits[ok_frag2])){
  print(i)
  max_rhat2[i] <- max(rhat(frag_results2$fits[ok_frag2][[i]]), na.rm = TRUE)
}

ok_frag2 <- ok_frag2[max_rhat2 < 1.05]

# visual diagnostics
label2 <- as_labeller(c("b_Intercept" = "beta[0]",
                        "b_c.lfs" = "beta[1]",
                        "sd_dataset_label__Intercept" = "sigma[0]",
                        "sd_dataset_label__c.lfs" = "sigma[1]",
                        "cor_dataset_label__Intercept__c.lfs" = "rho",
                        "b_sigma_Intercept" = "beta[0]^sigma",
                        "sd_dataset_label__sigma_Intercept" = "sigma[0]^sigma"), 
                      default = label_parsed)

pdf(paste0(wkdir, 'fragment-size-richness-relationships/figures/frag_sbc2.pdf'),
    width = 7, height = 5)

plot_rank_hist(frag_sbc_results_2,
               variables = c("b_c.lfs", #"b_Intercept"
                             "b_sigma_Intercept",
                             # "sd_dataset_label__Intercept",
                             "sd_dataset_label__c.lfs", 
                             # "cor_dataset_label__Intercept__c.lfs",
                             "sd_dataset_label__sigma_Intercept")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "cor_dataset_label__Intercept__c.lfs",
                       "b_sigma_Intercept",
                       "sd_dataset_label__sigma_Intercept")),
             labeller = label2)

plot_coverage(frag_sbc_results_2,
              variables = c("b_c.lfs", #"b_Intercept", 
                            "b_sigma_Intercept",
                            # "sd_dataset_label__Intercept",
                            "sd_dataset_label__c.lfs", 
                            # "cor_dataset_label__Intercept__c.lfs",
                            "sd_dataset_label__sigma_Intercept")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "cor_dataset_label__Intercept__c.lfs",
                       "b_sigma_Intercept",
                       "sd_dataset_label__sigma_Intercept")),
             labeller = label2)

plot_ecdf(frag_sbc_results_2,
          variables = c("b_c.lfs", #"b_Intercept", 
                        "b_sigma_Intercept",
                        # "sd_dataset_label__Intercept",
                        "sd_dataset_label__c.lfs", 
                        # "cor_dataset_label__Intercept__c.lfs",
                        "sd_dataset_label__sigma_Intercept")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "cor_dataset_label__Intercept__c.lfs",
                       "b_sigma_Intercept",
                       "sd_dataset_label__sigma_Intercept")),
             labeller = label2)

plot_sim_estimated(frag_sbc_results_2,
                   variables = c("b_c.lfs", #
                                 "b_sigma_Intercept",
                                 # "sd_dataset_label__Intercept",
                                 "sd_dataset_label__c.lfs", 
                                 # "cor_dataset_label__Intercept__c.lfs",
                                 "sd_dataset_label__sigma_Intercept")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs", 
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "cor_dataset_label__Intercept__c.lfs",
                       "b_sigma_Intercept",
                       "sd_dataset_label__sigma_Intercept")),
             labeller = label2, scales = 'free')

dev.off()

# save 
print(object.size(frag_results2[ok_frag2]), units = "Gb")

frag_sbc_results_2 <- frag_results2[ok_frag2]

save(frag_sbc_results_2,
     file = paste0(wkdir, '../model_fits/sbc-results/frag_sbc_results_goodfits_2.2.Rdata'))

