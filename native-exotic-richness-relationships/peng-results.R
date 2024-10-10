# results figures for native-exotic richness relationships

seed = 123

dat <- read_csv(paste0(wkdir, 'native-exotic-richness-relationships/data/DataS2_MetaAnalysisDatabase.csv'))


load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m1.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m2.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m3.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m4.Rdata'))
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-m5.Rdata'))

# stratified kfold
# need to create df of pointwise estimates for Yates fn
peng.kfold <- tibble(m1 = peng_m1_kfold$pointwise[,'elpd_kfold'], 
                     m2 = peng_sd_linear_kfold$pointwise[,'elpd_kfold'],
                     m3 = peng_sd_extent_kfold$pointwise[,'elpd_kfold'],
                     m4 = peng_sigma_linear_kfold$pointwise[,'elpd_kfold'],
                     m5 = peng_sigma_study_kfold$pointwise[,'elpd_kfold']) 

# leave-one-group-out cv
load(paste0(wkdir, 'native-exotic-richness-relationships/model-fits-CV-results/peng-logo.Rdata'))

peng.logo <- tibble(m1 = cv10g_m1$pointwise[,'elpd_kfold'],
                    m2 = cv10g_sd_linear$pointwise[,'elpd_kfold'],
                    m3 = cv10g_sd_extent$pointwise[,'elpd_kfold'],
                    m4 = cv10g_sigma_linear$pointwise[,'elpd_kfold'],
                    m5 = cv10g_sigma_study$pointwise[,'elpd_kfold'])

peng_ms_results <-
bind_rows(make_plot_data(peng.kfold) %>% 
            mutate(cv = 'kfold'),
          make_plot_data(peng.logo) %>% 
            mutate(cv = 'Leave-one-group-out')) %>% 
  filter(model!='m0') %>% 
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
                   labels = c(expression(paste('Fixed heterogeneity')),
                              expression(paste(log(omega), ' ~ log(grain)')),
                              expression(paste(log(omega), ' ~ extent')),
                              expression(paste(log(sigma), ' ~ log(grain)')),
                              expression(paste(log(sigma), ' ~ 1 + (1|Study)')))
  ) +
  labs(x = 'Model',
       y = expression(paste(Delta, bar(elpd[kfold])))) +
  coord_flip() +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, colour = 'black'),
        strip.text = element_text(hjust = 0, size = 10),
        axis.ticks = element_line(colour = 'black'),
        axis.text.y = element_text(hjust = 0))

# want to plot the two estimates of decline with grain size, and residual
# variation as a function of extent
sd_linear_decay <- peng_sd_linear %>% 
  spread_draws(b_logphi_Intercept,
               b_logphi_ln_Grain,
               seed = seed, 
               ndraws = 1000)

sigma_linear_decay <- peng_sigma_linear %>% 
  spread_draws(b_sigma_Intercept,
               b_sigma_ln_Grain,
               seed = seed,
               ndraws = 1000)

sd_extent <- peng_sd_extent %>% 
  gather_draws(`sd_Study:Case__Intercept:Extent[0-10]`,
               `sd_Study:Case__Intercept:Extent[10-100]`,
               `sd_Study:Case__Intercept:Extent[102-103]`,
               `sd_Study:Case__Intercept:Extent[103-104]`,
               `sd_Study:Case__Intercept:Extent[104-105]`,
               `sd_Study:Case__Intercept:Extent[105-106]`,
               `sd_Study:Case__Intercept:Extent[106-]`,
               seed = seed,
               ndraws = 1000)


ln_Grain = tibble(ln_Grain = seq(min(dat$ln_Grain), 
                                 max(dat$ln_Grain), 
                                 length.out = 1000)) %>% 
  nest(ln_Grain = ln_Grain)

# todo: better legend with parameter estimates, and model terms 
# https://stackoverflow.com/questions/21442629/superscript-and-subscript-the-same-character-in-an-expression
peng_sd_sigma_plot <-
bind_rows(sd_linear_decay %>% 
            select(-starts_with('.')) %>% 
            rename(intercept = b_logphi_Intercept,
                   slope = b_logphi_ln_Grain) %>% 
            mutate(model = 'sd'),
          sigma_linear_decay %>% 
            select(-starts_with('.')) %>% 
            rename(intercept = b_sigma_Intercept,
                   slope = b_sigma_ln_Grain) %>% 
            mutate(model = 'sigma')) %>% 
  group_by(model) %>% 
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
              group_by(model) %>% 
              mutate(rep = 1:n()) %>% 
              slice(1:100) %>% 
              ungroup() %>% 
              distinct(intercept, slope, ln_Grain, model, rep, 
                       .keep_all = TRUE) %>% 
              unnest(ln_Grain),
            aes(x = ln_Grain,
                y = exp(intercept + slope * ln_Grain), 
                colour = model, 
                # linetype = model,
                group = interaction(model, rep)),
            linewidth = 0.25,
            alpha = 0.25) +
  geom_line(data = . %>% 
              distinct(int, sl, ln_Grain, model, .keep_all = TRUE) %>% 
              unnest(ln_Grain),
            aes(x = ln_Grain,
                y = exp(int + sl * ln_Grain), 
                colour = model, 
                # linetype = model,
                group = model), 
            linewidth = 1) +
  scale_colour_manual(labels = c(expression(paste('log(',sigma,') ~ log(grain)')),
                               expression(paste('log(',omega,') ~ log(grain)'))),
                      values = c('sd' = '#003f5c',
                                 'sigma' = '#bc5090')) +
  labs(y = 'Residual variation [standard deviations]',
       x = expression(paste('Grain size [log(', m^2, ')]')),
       tag = '(c)') +
  theme_minimal() + 
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        panel.border = element_rect(fill = NA, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) 

# ggsave('~/Dropbox/1current/evidence-synthesis-heterogeneity/figures/sd-sigma-grain.pdf',
#          width = 100, height = 100, units = 'mm')

bind_rows(sd_linear_decay %>% 
            select(-starts_with('.')) %>% 
            rename(intercept = b_logphi_Intercept,
                   slope = b_logphi_ln_Grain) %>% 
            mutate(model = 'sd'),
          sigma_linear_decay %>% 
            select(-starts_with('.')) %>% 
            rename(intercept = b_sigma_Intercept,
                   slope = b_sigma_ln_Grain) %>% 
            mutate(model = 'sigma')) %>% 
  group_by(model) %>% 
  mutate(int = median(intercept),
         sl = median(slope),
         q95_int = quantile(intercept, probs = 0.95),
         q5_int = quantile(intercept, probs = 0.05),
         q95_sl = quantile(slope, probs = 0.95),
         q5_sl = quantile(slope, probs = 0.05)) %>% 
  ungroup() %>% 
  distinct(model, .keep_all = TRUE) %>% 
  ggplot() +
  geom_point(aes(x = sl, y = model)) +
  geom_linerange(aes(y = model,
                     xmin = q5_sl, xmax = q95_sl))

peng_sd_extent_plot <- ggplot() + 
  geom_violin(data = sd_extent,
              aes(x = .variable, y = .value)) + 
  geom_boxplot(data = sd_extent,
              aes(x = .variable, y = .value),
              width = 0.1) + 
  geom_text(data = dat %>%
              group_by(Extent) %>%
              summarise(n_obs = n_distinct(Study,Case)) %>%
              ungroup() %>%
              mutate(.variable = unique(sd_extent$.variable)),
            aes(x=.variable,
                y=-0.1,
                label=paste('n == ', n_obs)),
            size=2.5,
            parse = T) +
  scale_x_discrete(breaks = c('sd_Study:Case__Intercept:Extent[0-10]',
                              'sd_Study:Case__Intercept:Extent[10-100]',
                              'sd_Study:Case__Intercept:Extent[102-103]',
                              'sd_Study:Case__Intercept:Extent[103-104]',
                              'sd_Study:Case__Intercept:Extent[104-105]',
                              'sd_Study:Case__Intercept:Extent[105-106]',
                              'sd_Study:Case__Intercept:Extent[106-]'),
                   labels = c('(0-10)', '[10-100)', 
                              expression(paste('[',10^2, '-', 10^3, ')')), 
                              expression(paste('[',10^3, '-', 10^4, ')')),
                              expression(paste('[',10^4, '-', 10^5, ')')),
                              expression(paste('[',10^5, '-', 10^6, ')')),
                              expression(paste('[',10^6, '- )')))) +
  labs(x = expression(paste('Extent [', km^2, ']')),
       y = expression(paste('Residual variation (', omega, ')')),
       tag = '(d)') +
  theme_minimal() + 
  theme(panel.border = element_rect(fill = NA, colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text.x = element_text(hjust = 1, angle = 30))


# ggsave('~/Dropbox/1current/evidence-synthesis-heterogeneity/figures/sd-extent.pdf',
#        width = 150, height = 150, units = 'mm')


cowplot::plot_grid(peng_ms_results,
                   cowplot::plot_grid(peng_sd_sigma_plot,
                                      peng_sd_extent_plot,
                                      nrow = 1),
                   nrow = 2)

ggsave('~/Dropbox/1current/evidence-synthesis-heterogeneity/figures/Fig1.pdf',
       width = 200, height = 200, units = 'mm')


# compare parameter estimates
m1_pars <- peng_m1 %>% 
  gather_draws(`b_.*`, `sd_.*`,
               regex = TRUE,
               seed = seed, 
               ndraws = 1000) %>% 
  ungroup() %>% 
  mutate(model = '1.1')

peng_sd_grain_pars <- peng_sd_linear %>% 
  gather_draws(`b_.*`, `sd_.*`,
               regex = TRUE,
               seed = seed, 
               ndraws = 1000) %>% 
  # mutate(.variable = case_when(str_replace()
  ungroup() %>% 
  mutate(model = '1.2')

peng_sd_extent_pars <- peng_sd_extent %>% 
  gather_draws(`b_.*`, `sd_.*`,
               regex = TRUE,
               seed = seed, 
               ndraws = 1000) %>% 
  ungroup() %>% 
  mutate(model = '1.3')

peng_sigma_grain_pars <- peng_sigma_linear %>% 
  gather_draws(`b_.*`, `sd_.*`,
               regex = TRUE,
               seed = seed, 
               ndraws = 1000) %>% 
  ungroup() %>% 
  mutate(model = '1.4')

peng_sigma_study_pars <- peng_sigma_study %>% 
  gather_draws(`b_.*`, `sd_.*`,
               regex = TRUE,
               seed = seed, 
               ndraws = 1000) %>% 
  ungroup() %>% 
  mutate(model = '1.5')

bind_rows(m1_pars,
          peng_sd_grain_pars,
          peng_sd_extent_pars,
          peng_sigma_grain_pars,
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
                                 '1.2' = '#58508d',
                                 '1.3' = '#bc5090',
                                 '1.4' = '#ff6361',
                                 '1.5' = '#ffa600')) +
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

ggsave('~/Dropbox/1current/evidence-synthesis-heterogeneity/figures/FigS1.pdf',
       width = 125, height = 125, units = 'mm')  


# prediction intervals for best model
newdat <- peng_sigma_linear$data %>% 
  as_tibble()

posterior_predict(peng_sigma_linear,
                  newdata = newdat,
                  re_formula = NULL,
                  seed = seed,
                  ndraws = 100) %>% 
  as_tibble() %>% 
  pivot_longer(everything(),
               names_to = 'row',
               values_to = 'yhat') %>% 
  bind_cols(newdat %>% 
              slice(rep(row_number(), 100))) %>% 
  mutate(group = rep(1:100, each = nrow(newdat))) %>% 
  ggplot() +
  geom_point(aes(x = ln_Grain, y = z)) +
  stat_summary(aes(x = ln_Grain, y = yhat),
            geom = 'line',
            fun = 'median')
  
peng_sigma_study$data %>% 
  as_tibble() %>% 
  add_predicted_draws(peng_sigma_study,
                      seed = seed,
                      ndraws = 10) %>% 
  ungroup() %>% 
  ggplot() +
  geom_point(data = . %>% 
               distinct(.row, .keep_all = TRUE),
             aes(x = ln_Grain, y = z)) + 
  stat_lineribbon(aes(x = ln_Grain, y = .prediction, group = .draw), 
                  .width = c(.95, .50), alpha = 1/4) 
