# standardized_richness ~ f(fragment size) results of location-scale models
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')

# model fits
load(paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/Sstd-m1-m2.Rdata'))
load(paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/Sstd-m3.Rdata'))
load(paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/Sstd-m4.Rdata'))
load(paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/Sstd-m5.Rdata'))
     
# stratified kfold (k = 10)
load(paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/sstd_frag_kf10-4725045.Rdata'))

# leave-one-group-out cv
load(paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/sstd-logo-cv10g-1-4418331.Rdata'))
load(paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/sstd-logo-cv10g-2-4418336.Rdata'))
load(paste0(wkdir, 'fragment-size-richness-relationships/model-fits-CV-results/sstd-logo-cv10g-3-4418334.Rdata'))


# kfold
sstd_frag_kfold <- tibble(m1 = Sstd_fragSize_kf10$pointwise[,'elpd_kfold'],
                        m2 = Sstd_fragSize_sigma_kf10$pointwise[,'elpd_kfold'],
                        m3 = Sstd_fragSize_sigma_fs_kf10$pointwise[,'elpd_kfold'],
                        m4 = Sstd_fragSize_sigma_cor_kf10$pointwise[,'elpd_kfold'],
                        m5 = Sstd_fragSize_sigma_fs_cor_kf10$pointwise[,'elpd_kfold'],
) 

# logo 
sstd_frag_logo <- tibble(m1 = cv10g_sstd$pointwise[,'elpd_kfold'], 
                        m2 = cv10g_sstd_sigma$pointwise[,'elpd_kfold'], 
                        m3 = cv10g_sstd_sigma_fs$pointwise[,'elpd_kfold'], 
                        m4 = cv10g_sstd_sigma_cor$pointwise[,'elpd_kfold'], 
                        m5 = cv10g_sstd_sigma_fs_cor$pointwise[,'elpd_kfold'])


ab_panels <-
bind_rows(
  make_plot_data(sstd_frag_kfold) %>% 
    mutate(cv = 'kfold',
           metric = 'Richness'),
  make_plot_data(sstd_frag_logo) %>% 
    mutate(cv = 'Leave-one-group-out',
           metric = 'Richness')
) %>% 
  ggplot() +
  facet_grid(~factor(cv,
                     levels = c('kfold',
                                'Leave-one-group-out'),
                     labels = c('(a) K-fold (k = 10) cross validation',
                                '(b) Leave-one-group-out cross validation')), 
             scales = 'free_x', labeller = label_wrap_gen()) + 
  geom_point(aes(x = model, y = metric_diff, 
                 colour = metric,
                 group = metric),
             position = position_dodge(width = 0.5)) +
  geom_linerange(aes(x = model, ymin = metric_diff - se_mod, 
                     ymax = metric_diff + se_mod,
                     colour = metric,
                     group = metric),
                 position = position_dodge(width = 0.5)) +
  scale_x_discrete(breaks = c('m1', 
                              'm2',
                              'm3',
                              'm4',
                              'm5'),
                   labels = c(expression(paste(bold('Model 2.1: '), sigma, ' ~ 1')),
                              expression(paste(bold('Model 2.2: '), log(sigma),
                                               ' = 1 + (1 | Study)')),
                              expression(paste(bold('Model 2.3: '), log(sigma),
                                               ' = log(fragment size) + (log(fragment size) | Study)')),
                              expression(paste(bold('Model 2.4: '),
                                               log(sigma), ' = 1 + (1 |', rho,'| Study)')),
                              expression(paste(bold('Model 2.5: '), log(sigma),
                                               ' = log(fragment size) + (log(fragment size) |', rho,'| Study)')))) +
  scale_color_manual(values = c('richness' = 'black'),
                     guide = 'none') +
  labs(x = 'Model',
       y = expression(paste(Delta, bar(elpd)))) +
  coord_flip() + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = 'black', fill = NA),
        strip.text = element_text(hjust = 0, size = 9),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8))

# visualise variation in the relationship between patch size and residuals
fs_dat <- read_csv(paste0(wkdir, 'fragment-size-richness-relationships/data/2_biodiv_frag_fcont_10_mabund_as_is.csv'))
fs_dat <- fs_dat %>% 
  mutate(c.lfs = log(frag_size_num) - mean(log(frag_size_num)))
# metadata for fragmentation studies
fs_meta <- read_delim(paste0(wkdir, 'fragment-size-richness-relationships/data/new_meta_2_merge.csv'), 
                   delim = ';') %>% 
  rename(dataset_label = dataset_id)

seed <- 123

fs_plot <- fs_dat %>% 
  expand(c.lfs) %>% 
  nest(c.lfs = c.lfs)

sigma_fs_pop <- Sstd_lognorm_fragSize_sigma_fs_cor %>% 
  spread_draws(b_sigma_Intercept,
               b_sigma_c.lfs,
               seed = seed,
               ndraws = 100)

c_panel <- sigma_fs_pop %>% 
  select(-starts_with('.')) %>% 
  bind_cols(fs_plot) %>% 
  unnest(c.lfs) %>% 
  left_join(fs_dat %>% 
              distinct(c.lfs, frag_size_num)) %>% 
  mutate(yhat = exp(b_sigma_Intercept + b_sigma_c.lfs * c.lfs)) %>% 
  group_by(c.lfs) %>% 
  mutate(rep = 1:n()) %>% 
  ggplot() + 
  geom_line(aes(x = frag_size_num, y = yhat, group = rep),
            linewidth = 0.25, alpha = 0.25) +
  geom_line(aes(x = frag_size_num, y = yhat),
            stat = StatSummary,
            fun = median,  
            linewidth = 1) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", 
                                              scales::math_format(10^.x))) +
  labs(x = 'Fragment size [hectares]',
       y = expression(sigma),
       tag = '(c)') +
  theme_minimal() +
  theme(panel.border = element_rect(colour = 'black', fill = NA),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.tag = element_text(size = 9))

study_sigma <- Sstd_lognorm_fragSize_sigma %>% 
  spread_draws(b_c.lfs,
               r_dataset_label[dataset_label, effect],
               b_sigma_Intercept,
               r_dataset_label__sigma[dataset_label, term],
               seed = seed,
               ndraws = 1000) %>% 
  ungroup() %>% 
  left_join(fs_meta) %>% 
  mutate(model = 'varying sigma')

# taxa
d_panel <- ggplot() +
  geom_density_ridges_gradient(data = study_sigma,
                               aes(x = exp(b_sigma_Intercept + 
                                             r_dataset_label__sigma),
                                   y = taxa,
                                   # fill = stat(quantile)
                               ),
                               # quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9,
                               linetype = 0) +
  geom_rect(data = study_sigma %>% 
              distinct(b_sigma_Intercept) %>% 
              summarise(q05 = quantile(b_sigma_Intercept, prob = 0.05),
                        q95 = quantile(b_sigma_Intercept, prob = 0.95)),
            aes(xmin = exp(q05), xmax = exp(q95),
                ymin = -Inf, ymax = Inf),
            alpha = 0.2) +
  geom_point(data = study_sigma,
             aes(x = exp(b_sigma_Intercept + r_dataset_label__sigma), 
                 y = taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.25, shape = 18, colour = 'black') +
  geom_vline(data = study_sigma %>% 
               distinct(b_sigma_Intercept) %>% 
               summarise(med = median(b_sigma_Intercept)),
             aes(xintercept = exp(med)),
             linetype = 2) +
  geom_text(data = study_sigma %>%
              group_by(taxa) %>% 
              summarise(n_study = n_distinct(dataset_label)) %>% 
              ungroup(),
            aes(x=0.75, y=taxa, 
                label=paste('n[study] == ', n_study)),
            size=2.5,
            nudge_y = 0.25, parse = T) +
  labs(y = 'Taxon group',
       x = expression(sigma),
       tag = '(d)') +
  coord_cartesian(xlim = c(0, 1)) +
  theme_minimal() +
  theme(panel.border = element_rect(colour = 'black', fill = NA),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.tag = element_text(size = 9))
  
cowplot::plot_grid(ab_panels,
                   cowplot::plot_grid(c_panel,
                                      d_panel, 
                                      nrow = 1),
                   nrow = 2,
                   labels = 'Figure 2')

# ggsave(paste0(wkdir, 'fragment-size-richness-relationships/figures/Fig2.pdf'),
#        width = 180, height = 180, units = 'mm')

# want to check that all models support ecosystem decay for 
# standardized species richness
bind_rows(
  as_draws_df(Sstd_lognorm_fragSize, variable = 'b_c.lfs') %>% 
  mutate(model = 'm1'),
  as_draws_df(Sstd_lognorm_fragSize_sigma, variable = 'b_c.lfs') %>% 
    mutate(model = 'm2'),
  as_draws_df(Sstd_lognorm_fragSize_sigma_fs, variable = 'b_c.lfs') %>% 
    mutate(model = 'm3'),
  as_draws_df(Sstd_lognorm_fragSize_sigma_cor, variable = 'b_c.lfs') %>% 
    mutate(model = 'm4'),
  as_draws_df(Sstd_lognorm_fragSize_sigma_fs_cor, variable = 'b_c.lfs') %>% 
    mutate(model = 'm5')) %>% 
  ggplot() + 
  geom_density_ridges_gradient(aes(x = b_c.lfs,
                                   y = model,
                                   fill = stat(quantile)
                                   ),
                               quantiles = c(0.05, 0.95),
                               calc_ecdf = T,
                               scale = 0.9,
                               linetype = 0) +
  geom_point(aes(x = b_c.lfs, 
                 y = model),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.25, shape = 18, colour = 'black') +
  geom_vline(xintercept = 0,
             linetype = 1) +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#cccccc'),
                    guide = 'none') +
  scale_y_discrete(breaks = c('m1', 
                              'm2',
                              'm3',
                              'm4',
                              'm5'),
                   labels = c(expression(paste(bold('Model 2.1: '), sigma, ' ~ 1')),
                              expression(paste(bold('Model 2.2: '), log(sigma),
                                               ' = 1 + (1 | Study)')),
                              expression(paste(bold('Model 2.3: '), log(sigma),
                                               ' = log(fragment size) + (log(fragment size) | Study)')),
                              expression(paste(bold('Model 2.4: '),
                                               log(sigma), ' = 1 + (1 |', rho,'| Study)')),
                              expression(paste(bold('Model 2.5: '), log(sigma),
                                               ' = log(fragment size) + (log(fragment size) |', rho,'| Study)')))) +
  labs(x = expression(paste('Fragment size slope (', beta[1], ')')),
       y = 'Model') +
  theme_minimal() +
  theme(axis.title.x = element_text())

# ggsave(paste0(wkdir, 'fragment-size-richness-relationships/figures/FigSx-decay.pdf'),
#        width = 200, height = 100, units = 'mm')

# plot correlations from model 2.5
corr_post <- Sstd_lognorm_fragSize_sigma_fs_cor %>% 
  gather_draws(`cor_.*`,
               regex = TRUE) 

corr_post %>% 
  rename(value = .value,
         variable = .variable) %>% 
  select(-starts_with('.')) %>% 
  mutate(variable = str_remove(variable, 'cor_dataset_label__')) %>% 
  ggplot() +
  facet_wrap(~factor(variable,
                     levels = c('c.lfs__sigma_c.lfs',
                                'c.lfs__sigma_Intercept',
                                'Intercept__c.lfs',
                                'Intercept__sigma_c.lfs',
                                'Intercept__sigma_Intercept',
                                'sigma_Intercept__sigma_c.lfs'),
                     labels = c(expression(paste('(a)  ', 
                                                 rho[paste(sigma[1],sigma[1]^sigma)])),
                                expression(paste('(b)  ', 
                                                 rho[paste(sigma[1],sigma[0]^sigma)])), 
                                expression(paste('(c)  ', 
                                                 rho[paste(sigma[0],sigma[1])])),
                                expression(paste('(d)  ', 
                                                 rho[paste(sigma[0],sigma[1]^sigma)])),
                                expression(paste('(e)  ', 
                                                 rho[paste(sigma[0],sigma[0]^sigma)])),
                                expression(paste('(f)  ', 
                                                 rho[paste(sigma[0]^sigma,sigma[1]^sigma)])))),
             labeller = label_parsed,
             scales = 'free') +  
  geom_density_ridges_gradient(aes(x = value,
                                   y = 1,
                                   fill = after_stat(quantile)),
                               quantiles = c(0.05, 0.95),
                               calc_ecdf = T,
                               scale = 0.9,
                               linetype = 0) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#cccccc', '#969696', '#cccccc'),
                    labels = c('< 5 %', '5-95%', '> 95%')) +  
  labs(x = 'Correlation',
       y = 'Posterior density') +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, hjust = 0),
        legend.position = 'top')


# ggsave(paste0(wkdir, 'fragment-size-richness-relationships/figures/FigSx-correlations.pdf'),
#        width = 290, height = 200, units = 'mm')

# check interpretation of correlation, calculate study-level slopes
study_slopes <- Sstd_lognorm_fragSize_sigma_fs_cor %>% 
  spread_draws(b_c.lfs,
               r_dataset_label[dataset_label, term],
               b_sigma_c.lfs,
               r_dataset_label__sigma[dataset_label, term],
               seed = seed, 
               ndraws = 1000) %>% 
  ungroup()
  
study_slopes %>% 
  filter(term == 'c.lfs') %>% 
  mutate(mu_slope = b_c.lfs + r_dataset_label,
         sigma_slope = b_sigma_c.lfs + r_dataset_label__sigma) %>% 
  group_by(dataset_label) %>% 
  summarise(median_mu_slope = median(mu_slope),
            median_sigma_slope = median(sigma_slope)) %>% 
  ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_point(aes(x = median_mu_slope, y = median_sigma_slope)) +
  labs(x = expression(paste('Study-level ecosystem decay (', 
                            beta[1] + beta[`1i`], ')')),
       y = expression(paste('Study-level residual slope (', 
                            beta[1]^sigma + beta[`1i`]^sigma, ')'))) +
  theme_minimal()

# ggsave(paste0(wkdir,'fragment-size-richness-relationships/figures/FigSx-slope-slope.pdf'),
#        width = 130, height = 130, units = 'mm')
