## code to set working directory, 
# load packages, source a custom function from Yates et al 2021 Ecology
#  https://doi.org/10.1002/ecy.3475
# for calculating (and preparing for plotting) 
# the modified one-standard-error rule, and
# set up sbc (following sbc vignette)

wkdir <- '~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/'
setwd(wkdir)

library(tidyverse)
library(brms)
library(cmdstanr)
library(tidybayes)
library(ggridges)
library(SBC)
library(future)


make_plot_data <- function(metric_data, levels = names(metric_data)){
  best_model <- metric_data %>% 
    map_dbl(mean) %>% 
    which.max() %>% 
    names()
  tibble(model = factor(names(metric_data), levels = levels), 
         metric = metric_data %>% 
           map_dbl(mean),
         metric_diff = metric - max(metric),
         se = metric_data %>% 
           map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_diff = metric_data %>% 
           map(~ .x - metric_data[[best_model]]) %>% 
           map_dbl(~ .x %>% {sd(.)/sqrt(length(.))}),
         se_mod = sqrt(1 -cor(metric_data)[best_model,])*se[best_model])
}

# To use cmdstan (set to false to use rstan instead)
use_cmdstanr <- getOption("SBC.vignettes_cmdstanr", TRUE) 

if(use_cmdstanr) {
  library(cmdstanr)
} else {
  library(rstan)
  rstan_options(auto_write = TRUE)
}

# for work on my laptop 
options(mc.cores = parallel::detectCores() - 2)

# The fits are generally relatively fast,
# so forcing a minimum chunk size to reduce overhead of
# paralellization can decrease computation time.
options(SBC.min_chunk_size = 5)
