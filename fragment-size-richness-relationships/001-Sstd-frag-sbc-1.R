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

# priors also required to fully define the data generating process,
# model 2.1 in ms
frag_prior1 <- prior(normal(2.4,0.1), class = Intercept) + 
  prior(normal(0.06, 0.01), class = b) + 
  prior(normal(1,0.1), class = sd, coef = Intercept, group = dataset_label) +
  prior(normal(0,0.05), class = sd, coef = c.lfs, group = dataset_label) +
  prior(normal(0.3,0.05), class = sigma, lb = 0)

frag_generator1 <- SBC::SBC_generator_brms(S_std_mean ~ c.lfs + 
                                             (c.lfs | dataset_label), 
                                           family = lognormal(),
                                           data = frag_template_data, 
                                           prior = frag_prior1, 
                                           thin = 100, warmup = 15000, 
                                           refresh = 2000)

# generate the "fake" data sets: 
frag_datasets1 <- generate_datasets(frag_generator1, 500)

frag_backend1 <- SBC_backend_brms_from_generator(frag_generator1, chains = 2,
                                                 thin = 1,
                                                 warmup = 1000, iter = 2000,
                                                 prior = c(prior(normal(2.5,1), class = Intercept),
                                                           prior(normal(0,1), class = b),
                                                           prior(normal(0,1), class = sd),
                                                           prior(normal(0,1), class = sigma)),
                                                 control = list(adapt_delta = 0.9),
                                                 backend = 'cmdstanr')

frag_results1 <- compute_SBC(frag_datasets1, frag_backend1)

# inspect divergent transitions
hist(frag_results1$backend_diagnostics$n_divergent)

# retain fits with no divergent transitions
ok_frag1 <- frag_results1$backend_diagnostics %>% 
  filter(n_divergent == 0) %>% 
  pull(sim_id)

# rhat diagnostics in sbc output have NAs
# get non-NA rhats for inspection
max_rhat1 <- c()
for(i in 1:length(frag_results1$fits[ok_frag1])){
  print(i)
  max_rhat1[i] <- max(rhat(frag_results1$fits[ok_frag1][[i]]), na.rm = TRUE)
}

# reduce to fits with max rhat < 1.05
ok_frag1 <- ok_frag1[max_rhat1 < 1.05]

# neat labels for visual diagnostics
label1 <- as_labeller(c("b_Intercept" = "beta[0]",
                        "b_c.lfs" = "beta[1]",
                        "sd_dataset_label__Intercept" = "sigma[0]",
                        "sd_dataset_label__c.lfs" = "sigma[1]",
                        "cor_dataset_label__Intercept__c.lfs" = "rho",
                        "sigma" = "sigma"), 
                      default = label_parsed)

pdf(paste0(wkdir, 'fragment-size-richness-relationships/figures/frag_sbc1.pdf'),
    width = 7, height = 3)

plot_rank_hist(frag_sbc_results_1,
               variables = c("b_c.lfs", #"b_Intercept", 
                             #"sd_dataset_label__Intercept",
                             "sd_dataset_label__c.lfs", 
                             #"cor_dataset_label__Intercept__c.lfs",
                             "sigma")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "cor_dataset_label__Intercept__c.lfs",
                       "sigma")),
             labeller = label1)

plot_coverage(frag_sbc_results_1,
              variables = c("b_c.lfs",#"b_Intercept", 
                            # "sd_dataset_label__Intercept",
                            "sd_dataset_label__c.lfs", 
                            # "cor_dataset_label__Intercept__c.lfs",
                            "sigma")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "cor_dataset_label__Intercept__c.lfs",
                       "sigma")),
             labeller = label1)

plot_ecdf(frag_sbc_results_1,
          variables = c(#"b_Intercept", 
            "b_c.lfs",
            # "sd_dataset_label__Intercept",
            "sd_dataset_label__c.lfs", 
            # "cor_dataset_label__Intercept__c.lfs",
            "sigma")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "cor_dataset_label__Intercept__c.lfs",
                       "sigma")),
             labeller = label1)

plot_sim_estimated(frag_sbc_results_1,
                   variables = c("b_c.lfs",#"b_Intercept", 
                                 # "sd_dataset_label__Intercept",
                                 "sd_dataset_label__c.lfs", 
                                 # "cor_dataset_label__Intercept__c.lfs",
                                 "sigma")) +
  facet_wrap(~factor(variable, 
                     c("b_Intercept", "b_c.lfs",
                       "sd_dataset_label__Intercept",
                       "sd_dataset_label__c.lfs", 
                       "cor_dataset_label__Intercept__c.lfs",
                       "sigma")),
             labeller = label1, scales = 'free')
dev.off()

# check size and save 
print(object.size(frag_results1[ok_frag1]), units = "Gb")

frag_sbc_results_1 <- frag_results1[ok_frag1]

save(frag_sbc_results_1,
     file = paste0(wkdir, '../model_fits/sbc-results/frag_sbc_results_goodfits_2.1.Rdata'))