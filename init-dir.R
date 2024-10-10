## code to set working directory,
# load packages and source a custom function from Yates et al 2021
# for calculating mod modified one-standard-error rule 

wkdir <- '~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/'
setwd(wkdir)

library(tidyverse)
library(brms)
library(tidybayes)
library(ggridges)


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
