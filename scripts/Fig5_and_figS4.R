library(tidyverse)
library(magrittr)
library(broom)
library(ggpubr)
library(here)

protist_theme <- theme_pubr() +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
theme_set(protist_theme)

###5A###
data_5a <- data.frame(element = factor(c('Full model', 'Full model with interactions', 
                                  'Poly(A) signal', 'Upstream nucleotide', 
                                  'Downstream nucleotide'), 
                                  levels = c('Downstream nucleotide', 'Upstream nucleotide', 
                                             'Poly(A) signal', 'Full model with interactions', 
                                             'Full model'
                                             )),
                      coefficient = c(.51,.85,.58,.26,.39))

pdf(here('./figs/5A_raw.pdf'))
ggbarplot(data_5a, x = 'element', y = 'coefficient', orientation = 'horiz',
          fill = 'grey', color = 'grey', label = TRUE) +
  rremove('ylab')  +
  theme(text = element_text(size = 14)) + 
  geom_vline(xintercept = 4.5, linetype = 2)
dev.off()

gm_df <- read_tsv(here('./data/scores_8mers/gmur_8mers.tsv.gz'))
gm_4_models <- gm_df %>%
  mutate(signal = factor(case_when(grepl('AGTAAA', mer) ~ 'AGUAAA',
                                   grepl('TGTAAA', mer) ~ 'UGUAAA',
                                   grepl('AGTGAA', mer) ~ 'AGUGAA',
                                   TRUE ~ 'none'),
                         levels = c('none', 'AGUAAA', 'UGUAAA', 'AGUGAA')),
         upstream = factor(case_when(grepl('A[A|T]GT[A|G]AA', mer) ~'A',
                                     grepl('G[A|T]GT[A|G]AA', mer) ~'G',
                                     grepl('C[A|T]GT[A|G]AA', mer) ~'C',
                                     grepl('T[A|T]GT[A|G]AA', mer) ~'U',
                                     TRUE ~ 'none'), levels = c('none', 'A', 'G', 'C', 'U')),
         downstream = factor(case_when(grepl('[A|T]GT[A|G]AAA', mer) ~'A',
                                       grepl('[A|T]GT[A|G]AAG', mer) ~'G',
                                       grepl('[A|T]GT[A|G]AAC', mer) ~'C',
                                       grepl('[A|T]GT[A|G]AAT', mer) ~'U',
                                       TRUE ~ 'none'), levels = c('none', 'A', 'G', 'C', 'U')))
model_full_gm <- lm(score ~signal+upstream+downstream, 
                    data = gm_4_models)
model_full_r_sq_gm <- model_full_gm %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_signal_gm <- lm(score ~upstream+downstream, 
                                 data = gm_4_models)
model_reduced_no_signal_r_sq_gm <- model_reduced_no_signal_gm %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_up_gm <- lm(score ~signal+downstream, 
                             data = gm_4_models)
model_reduced_no_up_r_sq_gm <- model_reduced_no_up_gm %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_down_gm <- lm(score ~signal+upstream, 
                               data = gm_4_models)
model_reduced_no_down_r_sq_gm <- model_reduced_no_down_gm %>%
  glance() %>%
  dplyr::select('adj.r.squared')
r_sq_drop_gm <- tibble(term = c('full model', 'signal sequence', 'upstream nt', 'downstream nt'),
                       adj_rsq = c(as.numeric(model_full_r_sq_gm),
                                   as.numeric(model_reduced_no_signal_r_sq_gm),
                                   as.numeric(model_reduced_no_up_r_sq_gm),
                                   as.numeric(model_reduced_no_down_r_sq_gm)),
                       dataset = 'G. muris') %>%
  mutate(explained_variance = round(ifelse(term == 'full model', adj_rsq, 
                                           adj_rsq[1] -adj_rsq),digits = 2))
pdf(here('./figs/5A_1_raw.pdf'), width = 4, height = 5)
ggbarplot(r_sq_drop_gm, x = 'term', y = 'explained_variance', fill = 'gray',
          label = TRUE, lab.col = "black", lab.pos = "in", 
          lab.hjust = 0.5, lab.vjust = 1, lab.size = 5) +
  rremove('xlab') +
  ylab('adjusted R^2') +
  ylim(0,1) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.y = element_text(size = 14))
dev.off()

model_full_gm_int <- lm(score ~signal*upstream*downstream, 
                        data = gm_4_models)
model_full_r_sq_gm_int <- model_full_gm_int %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_signal_gm_int <- lm(score ~upstream*downstream, 
                                     data = gm_4_models)
model_reduced_no_signal_r_sq_gm_int <- model_reduced_no_signal_gm_int %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_up_gm_int <- lm(score ~signal*downstream, 
                                 data = gm_4_models)
model_reduced_no_up_r_sq_gm_int <- model_reduced_no_up_gm_int %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_down_gm_int <- lm(score ~signal*upstream, 
                                   data = gm_4_models)
model_reduced_no_down_r_sq_gm_int <- model_reduced_no_down_gm_int %>%
  glance() %>%
  dplyr::select('adj.r.squared')
r_sq_drop_gm_int <- tibble(term = c('full model', 'signal sequence', 'upstream nt', 'downstream nt'),
                           adj_rsq = c(as.numeric(model_full_r_sq_gm_int),
                                       as.numeric(model_reduced_no_signal_r_sq_gm_int),
                                       as.numeric(model_reduced_no_up_r_sq_gm_int),
                                       as.numeric(model_reduced_no_down_r_sq_gm_int)),
                           dataset = 'G. muris') %>%
  mutate(explained_variance = round(ifelse(term == 'full model', adj_rsq, 
                                           adj_rsq[1] -adj_rsq),digits = 2))
pdf(here('./figs/5A_2_raw.pdf'), width = 4, height = 5)
ggbarplot(r_sq_drop_gm_int, x = 'term', y = 'explained_variance', fill = 'gray',
          label = TRUE, lab.col = "black", lab.pos = "in", 
          lab.hjust = 0.5, lab.vjust = 1, lab.size = 5) +
  rremove('xlab') +
  ylab('adjusted R^2') +
  ylim(0,1) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.y = element_text(size = 14))
dev.off()


###S4A###
gm_df_full <- model_full_gm %>%
  tidy() %>%
  mutate(estimate = ifelse(term == '(Intercept)', estimate, estimate[1]+estimate)) %>%
  filter(!term == '(Intercept)') %>%
  mutate(dataset = factor(case_when(grepl('signal', term) ~ 'signal',
                                    grepl('upstream', term) ~ 'upstream',
                                    TRUE ~ 'downstream'), 
                          levels = c('upstream', 'signal', 'downstream')),
         term = factor(gsub('signal|upstream|downstream', '', term),
                       levels = c('A', 'G', 'C', 'U', 'AGUAAA', 'UGUAAA', 'AGUGAA')),
         color = factor(case_when(estimate > 0 & p.value < .05 ~ 1,
                                  estimate < 0 & p.value < .05 ~ 2,
                                  TRUE ~ 3)))

pdf(here('./figs/S4A_raw.pdf'), width = 6, height = 4) 
ggplot(gm_df_full , aes(x = term, y = estimate)) +
  geom_segment(aes(xend = term), yend = 0 ) +
  geom_point(size = 2, aes(color = color)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  facet_grid( . ~ dataset, scales = "free_x", space = "free") +
  rremove('xlab') +
  rremove('legend') +
  ylim(-2,2) +
  scale_color_manual(values = c('blue', 'orange', 'gray')) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()


###5B###
fit <- lm(score ~ signal*upstream*downstream, data = gm_4_models)
plot_data_AGUGAA <- fit$model %>%
  filter(signal == 'AGUGAA' & !upstream == 'none' & !downstream == 'none') %>%
  mutate(upstream = factor(upstream, levels = c('A', 'G', 'C', 'U')),
         downstream = factor(downstream, levels = c('A', 'G', 'C', 'U')))
pdf(here('figs/5B_raw.pdf'), width = 4, height = 3)
ggbarplot(plot_data_AGUGAA, x = 'downstream', y = 'score', 
          facet.by = 'upstream', fill = 'gray', nrow = 1,
          xlab = 'downstream nt', ylab = 'coefficient', main = 'AGUGAA') +
  geom_hline(yintercept = 0, linetype = 1) +
  theme(text = element_text(size = 14)) +
  ylim(-3.5,1)
dev.off()


###S4B###
plot_data_AGUAAA <- fit$model %>%
  filter(signal == 'AGUAAA' & !upstream == 'none' & !downstream == 'none') %>%
  mutate(upstream = factor(upstream, levels = c('A', 'G', 'C', 'U')),
         downstream = factor(downstream, levels = c('A', 'G', 'C', 'U')))
pdf(here('figs/S4B_raw.pdf'), width = 4, height = 3)
ggbarplot(plot_data_AGUAAA, x = 'downstream', y = 'score', 
          facet.by = 'upstream', fill = 'gray', nrow = 1,
          xlab = 'downstream nt', ylab = 'coefficient', main = 'AGUAAA') +
  geom_hline(yintercept = 0, linetype = 1) +
  theme(text = element_text(size = 14)) +
  ylim(-3.5,1)
dev.off()


###S4C###
plot_data_UGUAAA <- fit$model %>%
  filter(signal == 'UGUAAA' & !upstream == 'none' & !downstream == 'none') %>%
  mutate(upstream = factor(upstream, levels = c('A', 'G', 'C', 'U')),
         downstream = factor(downstream, levels = c('A', 'G', 'C', 'U')))
pdf(here('figs/S4C_raw.pdf'), width = 4, height = 3)
ggbarplot(plot_data_UGUAAA, x = 'downstream', y = 'score', 
          facet.by = 'upstream', fill = 'gray', nrow = 1,
          xlab = 'downstream nt', ylab = 'coefficient', main = 'UGUAAA') +
  geom_hline(yintercept = 0, linetype = 1) +
  theme(text = element_text(size = 14)) +
  ylim(-3.5,1)
dev.off()


###5C###
data_5c <- read_tsv(here('./data/prem_cl/gmur_prem_cl_4_fig5c.tsv.gz')) %>%
  mutate(signal = factor(signal, levels = c('AGUAAA', 'UGUAAA', 'AGUGAA')),
         prem_cl = factor(prem_cl))

pdf(here('./figs/5C_raw.pdf'))
ggbarplot(data_5c, x = 'signal', y = 'counts', 
          fill = 'prem_cl', palette = c('red', 'gray'),
          position = position_dodge()) +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  rremove('legend.title') +
  rremove('xlab') +
  ylim(0,80)
dev.off()
