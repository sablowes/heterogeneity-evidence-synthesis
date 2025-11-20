# results figures for native-exotic species richness relationship case study
source('~/Dropbox/1current/evidence-synthesis-heterogeneity/heterogeneity-evidence-synthesis/init-dir.R')
seed = 123

# data from local directory
dat <- read_csv(paste0(wkdir, 'native-exotic-richness-relationships/data/DataS2_MetaAnalysisDatabase.csv'))

load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m1.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m2.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m3.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m4.Rdata'))


# stratified kfold
# need to create df of pointwise estimates for Yates fn
peng.kfold <- tibble(m1 = peng_m1_kfold$pointwise[,'elpd_kfold'], 
                     m2 = peng_sigma_linear_kfold$pointwise[,'elpd_kfold'],
                     m3 = peng_sigma_extent_kfold$pointwise[,'elpd_kfold'],
                     m4 = peng_sigma_study_kfold$pointwise[,'elpd_kfold'])

# leave-one-group-out cv
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m1-logo.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m2-logo.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m3-logo.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m4-logo.Rdata'))

peng.logo <- tibble(m1 = cv10g_m1$pointwise[,'elpd_kfold'],
                    m2 = cv10g_m2$pointwise[,'elpd_kfold'],
                    m3 = cv10g_m3$pointwise[,'elpd_kfold'],
                    m4 = cv10g_m4$pointwise[,'elpd_kfold'])

peng_ms_results <-
bind_rows(make_plot_data(peng.kfold) %>% 
            mutate(cv = 'kfold'),
          make_plot_data(peng.logo) %>%
            mutate(cv = 'Leave-one-group-out')
          ) %>% 
  ggplot() +
  facet_grid(~factor(cv,
                     levels = c('kfold',
                                'Leave-one-group-out'),
                     labels = c('(a) Stratified k-fold (k = 10) cross validation',
                                '(b) Leave-one-group-out cross validation')), 
             scales = 'free_x', labeller = label_wrap_gen(width = 120)) + 
  geom_point(aes(x = model, y = metric_diff)) +
  geom_linerange(aes(x = model, ymin = metric_diff - se_mod, 
                     ymax = metric_diff + se_mod)) +
  scale_x_discrete(breaks = c('m1', 'm2', 'm3', 'm4', 'm5'),
                   labels = c(
                     expression(paste(bold('Model 1.1: '), sigma, ' ~ 1')),
                     expression(paste(bold('Model 1.2: '), log(sigma), ' ~ log(grain)')),
                     expression(paste(bold('Model 1.3: '), log(sigma), ' ~ extent')),
                     expression(paste(bold('Model 1.4: '), log(sigma), ' ~ 1 + (1|Study)')))
  ) +
  labs(x = 'Model',
       y = expression(paste(Delta, bar(elpd)))) +
  coord_flip() +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = 'black'),
        strip.text = element_text(hjust = 0, size = 9),
        axis.ticks = element_line(colour = 'black'),
        axis.text.y = element_text(hjust = 0, size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 8))

# want to plot estimate of decline in residual variation with grain size, 
# and residual variation as a function of extent
sigma_linear_decay <- peng_sigma_linear %>% 
  spread_draws(b_sigma_Intercept,
               b_sigma_ln_Grain,
               seed = seed,
               ndraws = 1000)

sigma_linear_decay %>% 
  summarise(q50 = quantile(b_sigma_ln_Grain, prob = 0.5),
            q05 = quantile(b_sigma_ln_Grain, prob = 0.05),
            q95 = quantile(b_sigma_ln_Grain, prob = 0.95))

sigma_extent <- peng_sigma_extent %>% 
  gather_draws(b_sigma_Extent0M10,
               b_sigma_Extent10M100,
               b_sigma_Extent102M103,
               b_sigma_Extent103M104,
               b_sigma_Extent104M105,
               b_sigma_Extent105M106,
               b_sigma_Extent106M,
               seed = seed,
               ndraws = 1000)

sigma_extent$.variable <- factor(sigma_extent$.variable,
                                 levels = c('b_sigma_Extent0M10',
                                            'b_sigma_Extent10M100',
                                            'b_sigma_Extent102M103',
                                            'b_sigma_Extent103M104',
                                            'b_sigma_Extent104M105',
                                            'b_sigma_Extent105M106',
                                            'b_sigma_Extent106M'))

ln_Grain = tibble(ln_Grain = seq(min(dat$ln_Grain), 
                                 max(dat$ln_Grain), 
                                 length.out = 1000)) %>% 
  nest(ln_Grain = ln_Grain)

# todo: better legend with parameter estimates, and model terms 
# https://stackoverflow.com/questions/21442629/superscript-and-subscript-the-same-character-in-an-expression
peng_sigma_grain_plot <-
  sigma_linear_decay %>% 
  select(-starts_with('.')) %>% 
  rename(intercept = b_sigma_Intercept,
         slope = b_sigma_ln_Grain) %>% 
  mutate(model = 'sigma') %>% 
  mutate(int = median(intercept),
         sl = median(slope),
         q95_int = quantile(intercept, probs = 0.95),
         q5_int = quantile(intercept, probs = 0.05),
         q95_sl = quantile(slope, probs = 0.95),
         q5_sl = quantile(slope, probs = 0.05)) %>% 
  ungroup() %>% 
  bind_cols(ln_Grain) %>% 
ggplot() +
  geom_line(data = . %>% 
              mutate(rep = 1:n()) %>% 
              slice(1:100) %>% 
              ungroup() %>% 
              distinct(intercept, slope, ln_Grain, rep, 
                       .keep_all = TRUE) %>% 
              unnest(ln_Grain),
            aes(x = ln_Grain,
                y = exp(intercept + slope * ln_Grain), 
                group = interaction(rep)),
            linewidth = 0.25,
            alpha = 0.25) +
  geom_line(data = . %>% 
              distinct(int, sl, ln_Grain, model, .keep_all = TRUE) %>% 
              unnest(ln_Grain),
            aes(x = ln_Grain,
                y = exp(int + sl * ln_Grain),
                colour = model), 
            linewidth = 1) +
    scale_colour_manual(name = '',
                        labels = c(expression(paste('log(',sigma,') ~ log(grain)'))),
                        values = c('sigma' = 'black')) +
  labs(y = expression(paste('Residual variation (', sigma, ')')),
       x = expression(paste('Grain size [log(', m^2, ')]')),
       tag = '(c)') +
  theme_minimal() + 
  theme(legend.position = 'none',
        legend.justification = c(1,1),
        panel.border = element_rect(fill = NA, colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.tag = element_text(size = 9)) 

peng_sigma_extent_plot <-
  ggplot() + 
  geom_violin(data = sigma_extent,
              aes(x = .variable, y = exp(.value))) + 
  geom_boxplot(data = sigma_extent,
              aes(x = .variable, y = exp(.value)),
              width = 0.1) + 
  geom_text(data = dat %>%
              group_by(Extent) %>%
              summarise(n_obs = n_distinct(Study,Case)) %>%
              ungroup() %>%
              mutate(.variable = unique(sigma_extent$.variable)),
            aes(x=.variable,
                y=-0.01,
                label=paste('n == ', n_obs)),
            size=2.5,
            parse = T) +
  scale_x_discrete(breaks = c('b_sigma_Extent0M10',
                              'b_sigma_Extent10M100',
                              'b_sigma_Extent102M103',
                              'b_sigma_Extent103M104',
                              'b_sigma_Extent104M105',
                              'b_sigma_Extent105M106',
                              'b_sigma_Extent106M'),
                   labels = c('(0-10)', '[10-100)', 
                              expression(paste('[',10^2, '-', 10^3, ')')), 
                              expression(paste('[',10^3, '-', 10^4, ')')),
                              expression(paste('[',10^4, '-', 10^5, ')')),
                              expression(paste('[',10^5, '-', 10^6, ')')),
                              expression(paste('[',10^6, '- )')))) +
  # scale_y_continuous(trans = 'exp') +
  labs(x = expression(paste('Extent [', km^2, ']')),
       y = expression(paste('Residual variation (', sigma, ')')),
       tag = '(d)') +
  theme_minimal() + 
  theme(panel.border = element_rect(fill = NA, colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text.x = element_text(hjust = 1, angle = 30, size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8),
        plot.tag = element_text(size = 9))

cowplot::plot_grid(peng_ms_results,
                   cowplot::plot_grid(peng_sigma_grain_plot,
                                      peng_sigma_extent_plot,
                                      nrow = 1),
                   nrow = 2,
                   labels = 'Figure 1')

# ggsave('native-exotic-richness-relationships/figures/Fig1.pdf',
#        width = 180, height = 180, units = 'mm')


# compare parameter estimates
m1_pars <- peng_m1 %>% 
  gather_draws(`b_.*`, `sd_.*`,
               regex = TRUE,
               seed = seed, 
               ndraws = 1000) %>% 
  ungroup() %>% 
  mutate(model = '1.1')

peng_sigma_grain_pars <- peng_sigma_linear %>% 
  gather_draws(`b_.*`, `sd_.*`,
               regex = TRUE,
               seed = seed, 
               ndraws = 1000) %>% 
  ungroup() %>% 
  mutate(model = '1.2')

peng_sigma_extent_pars <- peng_sigma_extent %>% 
  gather_draws(`b_.*`, `sd_.*`,
               regex = TRUE,
               seed = seed, 
               ndraws = 1000) %>% 
  ungroup() %>% 
  mutate(model = '1.3')

peng_sigma_study_pars <- peng_sigma_study %>% 
  gather_draws(`b_.*`, `sd_.*`,
               regex = TRUE,
               seed = seed, 
               ndraws = 1000) %>% 
  ungroup() %>% 
  mutate(model = '1.4')

bind_rows(m1_pars,
          peng_sigma_grain_pars,
          peng_sigma_extent_pars,
          peng_sigma_study_pars) %>% 
  # distinct(.variable)
  filter(.variable %in% c('b_Intercept', 'b_ln_Grain',
                          'b_muy_Intercept', 'b_muy_ln_Grain',
                          'sd_Study__Intercept', 'sd_Study__vint1_Intercept'
                          # 'b_logphi_Intercept', 'b_logphi_Intercept',
                          # 'b_sigma_Intercept', 'b_sigma_ln_Grain'
                          )) %>% 
  mutate(.variable = str_replace(.variable, 'b_muy', 'b'),
         .variable = str_replace(.variable, 'sd_Study__vint1_Intercept', 'sd_Study__Intercept')) %>% 
  # distinct(.variable, model)
  ggplot() +
  stat_summary(aes(x = .variable, y = .value, 
                 group = model,
                 colour = model),
               fun = 'median',
               fun.min = function(x) quantile(x, probs = 0.05),
               fun.max = function(x) quantile(x, probs = 0.95),
               position = position_dodge(width = 0.4)) +
  scale_colour_manual(values = c('1.1' = '#003f5c',
                                 '1.2' = '#444e86',
                                 '1.3' = '#955196',
                                 '1.4' = '#dd5182')) +
  scale_x_discrete(labels = c(expression(beta[0]),
                              expression(beta[1]),
                              expression(tau)),
                   name = 'Parameter') +
  labs(y = 'Estimate') +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.direction = 'horizontal',
        panel.border = element_rect(fill = NA, colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text = element_text(size = 12))

# ggsave('native-exotic-richness-relationships/figures/FigSx-shared-pars.pdf',
#        width = 125, height = 125, units = 'mm')  

