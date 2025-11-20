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

# model 2.3
frag_priors3 <- prior(normal(2.5,0.1), class = Intercept) +
  prior(normal(0.05, 0.02), class = b) + 
  prior(normal(1,0.1), class = sd, coef = Intercept, group = dataset_label) +
  prior(normal(0.05,0.05), class = sd, coef = c.lfs, group = dataset_label) +
  prior(normal(-1.25, 0.1), class = Intercept, dpar = sigma) +
  prior(normal(0,0.05), class = b, dpar = sigma, coef = c.lfs) +
  prior(normal(0.5,0.05), class = sd, dpar = sigma, coef = Intercept, group = dataset_label) +
  prior(normal(0.05,0.05), class = sd, dpar = sigma, coef = c.lfs, group = dataset_label) 


frag_generator3 <- SBC::SBC_generator_brms(bf(S_std_mean ~ c.lfs + (c.lfs | dataset_label),
                                              sigma ~ c.lfs + (c.lfs | dataset_label)),
                                           family = lognormal(),
                                           data = frag_template_data, 
                                           prior = frag_priors3, 
                                           thin = 100, warmup = 15000, 
                                           refresh = 2000)

frag_datasets3 <- generate_datasets(frag_generator3, 200)

frag_backend3 <- SBC_backend_brms_from_generator(frag_generator3, chains = 2,
                                                 thin = 1,
                                                 warmup = 1000, iter = 2000,
                                                 prior = c(prior(normal(2.5,1), class = Intercept),
                                                           prior(normal(0,1), class = b),
                                                           prior(normal(0,1), class = sd),
                                                           prior(normal(0,1), class = sd, dpar = sigma),
                                                           prior(normal(-1,1), class = Intercept, dpar = sigma),
                                                           prior(normal(0,1), class = b, dpar = sigma)),
                                                 backend = 'cmdstanr',
                                                 control = list(adapt_delta = 0.9))

frag_results3 <- compute_SBC(frag_datasets3, frag_backend3)

# identify fits with no divergent transitions
ok_frag3 <- frag_results3$backend_diagnostics %>% 
  filter(n_divergent==0) %>% 
  pull(sim_id)

# rhat diagnostics in sbc output have NAs
max_rhat3 <- c()
for(i in 1:length(frag_results3$fits[ok_frag3])){
  print(i)
  max_rhat3[i] <- max(rhat(frag_results3$fits[ok_frag3][[i]]), na.rm = TRUE)
}

# remove the rhat > 1.05
ok_frag3 <- ok_frag3[max_rhat3 < 1.05]

# visual diagnostics
label3 <- as_labeller(c("b_Intercept" = "beta[0]",
                        "b_c.lfs" = "beta[1]",
                        "sd_dataset_label__Intercept" = "sigma[0]",
                        "sd_dataset_label__c.lfs" = "sigma[1]",
                        "b_sigma_Intercept" = "beta[0]^sigma",
                        "b_sigma_c.lfs" = "beta[1]^sigma",
                        "sd_dataset_label__sigma_Intercept" = "sigma[0]^sigma",
                        "sd_dataset_label__sigma_c.lfs" = "sigma[1]^sigma"), 
                      default = label_parsed)

frag_sbc_results_3 <- frag_results3[ok_frag3]

pdf(paste0(wkdir, 'fragment-size-richness-relationships/figures/frag_sbc3.pdf'),
    width = 7, height = 5)

plot_rank_hist(frag_sbc_results_3,
               variables = c("b_c.lfs", #"b_Intercept", 
                             "b_sigma_Intercept",
                             "b_sigma_c.lfs",
                             # "sd_dataset_label__Intercept",
                             "sd_dataset_label__c.lfs", 
                             "sd_dataset_label__sigma_Intercept",
                             "sd_dataset_label__sigma_c.lfs")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "b_sigma_Intercept",
                       "b_sigma_c.lfs",
                       "sd_dataset_label__sigma_Intercept",
                       "sd_dataset_label__sigma_c.lfs")),
             labeller = label3)


plot_coverage(frag_sbc_results_3,
              variables = c("b_c.lfs", #"b_Intercept", 
                            "b_sigma_Intercept",
                            "b_sigma_c.lfs",
                            # "sd_dataset_label__Intercept",
                            "sd_dataset_label__c.lfs", 
                            "sd_dataset_label__sigma_Intercept",
                            "sd_dataset_label__sigma_c.lfs")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "b_sigma_Intercept",
                       "b_sigma_c.lfs",
                       "sd_dataset_label__sigma_Intercept",
                       "sd_dataset_label__sigma_c.lfs")),
             labeller = label3)

plot_ecdf(frag_sbc_results_3,
          variables = c("b_c.lfs", #"b_Intercept", 
                        "b_sigma_Intercept",
                        "b_sigma_c.lfs",
                        # "sd_dataset_label__Intercept",
                        "sd_dataset_label__c.lfs", 
                        "sd_dataset_label__sigma_Intercept",
                        "sd_dataset_label__sigma_c.lfs")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "b_sigma_Intercept",
                       "b_sigma_c.lfs",
                       "sd_dataset_label__sigma_Intercept",
                       "sd_dataset_label__sigma_c.lfs")),
             labeller = label3)

plot_sim_estimated(frag_sbc_results_3,
                   variables = c("b_c.lfs", #"b_Intercept", 
                                 "b_sigma_Intercept",
                                 "b_sigma_c.lfs",
                                 # "sd_dataset_label__Intercept",
                                 "sd_dataset_label__c.lfs", 
                                 "sd_dataset_label__sigma_Intercept",
                                 "sd_dataset_label__sigma_c.lfs")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "b_sigma_Intercept",
                       "b_sigma_c.lfs",
                       "sd_dataset_label__sigma_Intercept",
                       "sd_dataset_label__sigma_c.lfs")),
             labeller = label3, scales = 'free')

dev.off()

# save 
print(object.size(frag_results3[ok_frag3]), units = "Gb")

save(frag_sbc_results_3,
     file = paste0(wkdir, '../model_fits/sbc-results/frag_sbc_results_goodfits_2.3.Rdata'))
