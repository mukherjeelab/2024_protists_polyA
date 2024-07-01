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

ng_df <- read_tsv(here('./data/scores_8mers/neg_8mers.tsv.gz'))
tf_df <- read_tsv(here('./data/scores_8mers/trif_8mers.tsv.gz'))
tv_df <- read_tsv(here('./data/scores_8mers/trich_8mers.tsv.gz'))

ng_4_models <- ng_df %>%
  mutate(signal = factor(case_when(grepl('AATAAA', mer) ~ 'AAUAAA',
                                   grepl('ATTAAA', mer) ~ 'AUUAAA',
                                   TRUE ~ 'none'),
                         levels = c('none', 'AAUAAA', 'AUUAAA')),
         upstream = factor(case_when(grepl('AA[A|T]TAAA', mer) ~'A',
                                     grepl('GA[A|T]TAAA', mer) ~'G',
                                     grepl('CA[A|T]TAAA', mer) ~'C',
                                     grepl('TA[A|T]TAAA', mer) ~'U',
                                     TRUE ~ 'none'), levels = c('none', 'A', 'G', 'C', 'U')),
         downstream = factor(case_when(grepl('A[A|T]TAAAA', mer) ~'A',
                                       grepl('A[A|T]TAAAG', mer) ~'G',
                                       grepl('A[A|T]TAAAC', mer) ~'C',
                                       grepl('A[A|T]TAAAT', mer) ~'U',
                                       TRUE ~ 'none'), levels = c('none', 'A', 'G', 'C', 'U')))
tf_4_models <- tf_df %>%
  mutate(signal = factor(case_when(grepl('AATAAA', mer) ~ 'AAUAAA',
                                   grepl('ATTAAA', mer) ~ 'AUUAAA',
                                   TRUE ~ 'none'),
                         levels = c('none', 'AAUAAA', 'AUUAAA')),
         upstream = factor(case_when(grepl('AA[A|T]TAAA', mer) ~'A',
                                     grepl('GA[A|T]TAAA', mer) ~'G',
                                     grepl('CA[A|T]TAAA', mer) ~'C',
                                     grepl('TA[A|T]TAAA', mer) ~'U',
                                     TRUE ~ 'none'), levels = c('none', 'A', 'G', 'C', 'U')),
         downstream = factor(case_when(grepl('A[A|T]TAAAA', mer) ~'A',
                                       grepl('A[A|T]TAAAG', mer) ~'G',
                                       grepl('A[A|T]TAAAC', mer) ~'C',
                                       grepl('A[A|T]TAAAT', mer) ~'U',
                                       TRUE ~ 'none'), levels = c('none', 'A', 'G', 'C', 'U')))
tv_4_models <- tv_df %>%
  mutate(signal = factor(case_when(grepl('AATAAA', mer) ~ 'AAUAAA',
                                   grepl('ATTAAA', mer) ~ 'AUUAAA',
                                   TRUE ~ 'none'),
                         levels = c('none', 'AAUAAA', 'AUUAAA')),
         upstream = factor(case_when(grepl('AA[A|T]TAAA', mer) ~'A',
                                     grepl('GA[A|T]TAAA', mer) ~'G',
                                     grepl('CA[A|T]TAAA', mer) ~'C',
                                     grepl('TA[A|T]TAAA', mer) ~'U',
                                     TRUE ~ 'none'), levels = c('none', 'A', 'G', 'C', 'U')),
         downstream = factor(case_when(grepl('A[A|T]TAAAA', mer) ~'A',
                                       grepl('A[A|T]TAAAG', mer) ~'G',
                                       grepl('A[A|T]TAAAC', mer) ~'C',
                                       grepl('A[A|T]TAAAT', mer) ~'U',
                                       TRUE ~ 'none'), levels = c('none', 'A', 'G', 'C', 'U')))
model_full_ng <- lm(score ~signal+upstream+downstream, 
                    data = ng_4_models)
model_full_r_sq_ng <- model_full_ng %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_signal_ng <- lm(score ~upstream+downstream, 
                                 data = ng_4_models)
model_reduced_no_signal_r_sq_ng <- model_reduced_no_signal_ng %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_up_ng <- lm(score ~signal+downstream, 
                             data = ng_4_models)
model_reduced_no_up_r_sq_ng <- model_reduced_no_up_ng %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_down_ng <- lm(score ~signal+upstream, 
                               data = ng_4_models)
model_reduced_no_down_r_sq_ng <- model_reduced_no_down_ng %>%
  glance() %>%
  dplyr::select('adj.r.squared')
r_sq_drop_ng <- tibble(term = c('full model', 'signal sequence', 'upstream nt', 'downstream nt'),
                       adj_rsq = c(as.numeric(model_full_r_sq_ng),
                                   as.numeric(model_reduced_no_signal_r_sq_ng),
                                   as.numeric(model_reduced_no_up_r_sq_ng),
                                   as.numeric(model_reduced_no_down_r_sq_ng)),
                       dataset = 'N. gruberi') %>%
  mutate(explained_variance = ifelse(term == 'full model', adj_rsq, adj_rsq[1]-adj_rsq))
model_full_tf <- lm(score ~signal+upstream+downstream, 
                    data = tf_4_models)
model_full_r_sq_tf <- model_full_tf %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_signal_tf <- lm(score ~upstream+downstream, 
                                 data = tf_4_models)
model_reduced_no_signal_r_sq_tf <- model_reduced_no_signal_tf %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_up_tf <- lm(score ~signal+downstream, 
                             data = tf_4_models)
model_reduced_no_up_r_sq_tf <- model_reduced_no_up_tf %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_down_tf <- lm(score ~signal+upstream, 
                               data = tf_4_models)
model_reduced_no_down_r_sq_tf <- model_reduced_no_down_tf %>%
  glance() %>%
  dplyr::select('adj.r.squared')
r_sq_drop_tf <- tibble(term = c('full model', 'signal sequence', 'upstream nt', 'downstream nt'),
                       adj_rsq = c(as.numeric(model_full_r_sq_tf),
                                   as.numeric(model_reduced_no_signal_r_sq_tf),
                                   as.numeric(model_reduced_no_up_r_sq_tf),
                                   as.numeric(model_reduced_no_down_r_sq_tf)),
                       dataset = 'T. foetus')  %>%
  mutate(explained_variance = ifelse(term == 'full model', adj_rsq, adj_rsq[1]-adj_rsq))
model_full_tv <- lm(score ~signal+upstream+downstream, 
                    data = tv_4_models)
model_full_r_sq_tv <- model_full_tv %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_signal_tv <- lm(score ~upstream+downstream, 
                                 data = tv_4_models)
model_reduced_no_signal_r_sq_tv <- model_reduced_no_signal_tv %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_up_tv <- lm(score ~signal+downstream, 
                             data = tv_4_models)
model_reduced_no_up_r_sq_tv <- model_reduced_no_up_tv %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_down_tv <- lm(score ~signal+upstream, 
                               data = tv_4_models)
model_reduced_no_down_r_sq_tv <- model_reduced_no_down_tv %>%
  glance() %>%
  dplyr::select('adj.r.squared')
r_sq_drop_tv <- tibble(term = c('full model', 'signal sequence', 'upstream nt', 'downstream nt'),
                       adj_rsq = c(as.numeric(model_full_r_sq_tv),
                                   as.numeric(model_reduced_no_signal_r_sq_tv),
                                   as.numeric(model_reduced_no_up_r_sq_tv),
                                   as.numeric(model_reduced_no_down_r_sq_tv)),
                       dataset = 'T. vaginalis')  %>%
  mutate(explained_variance = ifelse(term == 'full model', adj_rsq, adj_rsq[1]-adj_rsq))

r_sq_drop_all <- rbind(r_sq_drop_ng, r_sq_drop_tf, r_sq_drop_tv) %>%
  mutate(term = factor(term, 
                       levels = c('full model','signal sequence', 'upstream nt','downstream nt')),
         explained_variance = round(explained_variance, digits = 2))

pdf(here('./figs/6A_raw.pdf'), width = 10, height = 5)
ggbarplot(r_sq_drop_all, x = 'term', y = 'explained_variance', fill = 'gray',
          facet.by = 'dataset', label = TRUE,
          lab.col = "black", lab.pos = "in", lab.hjust = 0.5, lab.vjust = 1) +
  rremove('xlab') +
  ylab('adjusted R^2') +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14))
dev.off()

ng_df_full <- model_full_ng %>%
  tidy() %>%
  mutate(estimate = ifelse(term == '(Intercept)', estimate, estimate[1]+estimate)) %>%
  filter(!term == '(Intercept)') %>%
  mutate(dataset = factor(case_when(grepl('signal', term) ~ 'signal',
                                    grepl('upstream', term) ~ 'upstream',
                                    TRUE ~ 'downstream'), 
                          levels = c('upstream', 'signal', 'downstream')),
         term = factor(gsub('signal|upstream|downstream', '', term),
                       levels = c('A', 'G', 'C', 'U', 'AAUAAA', 'AUUAAA')),
         color = factor(case_when(estimate > 0 & p.value < .05 ~ 1,
                                  estimate < 0 & p.value < .05 ~ 2,
                                  TRUE ~ 3)))

tf_df_full <- model_full_tf %>%
  tidy() %>%
  mutate(estimate = ifelse(term == '(Intercept)', estimate, estimate[1]+estimate)) %>%
  filter(!term == '(Intercept)') %>%
  mutate(dataset = factor(case_when(grepl('signal', term) ~ 'signal',
                                    grepl('upstream', term) ~ 'upstream',
                                    TRUE ~ 'downstream'), 
                          levels = c('upstream', 'signal', 'downstream')),
         term = factor(gsub('signal|upstream|downstream', '', term),
                       levels = c('A', 'G', 'C', 'U', 'AAUAAA', 'AUUAAA')),
         color = factor(case_when(estimate > 0 & p.value < .05 ~ 1,
                                  estimate < 0 & p.value < .05 ~ 2,
                                  TRUE ~ 3)))

tv_df_full <- model_full_tv %>%
  tidy() %>%
  mutate(estimate = ifelse(term == '(Intercept)', estimate, estimate[1]+estimate)) %>%
  filter(!term == '(Intercept)') %>%
  mutate(dataset = factor(case_when(grepl('signal', term) ~ 'signal',
                                    grepl('upstream', term) ~ 'upstream',
                                    TRUE ~ 'downstream'), 
                          levels = c('upstream', 'signal', 'downstream')),
         term = factor(gsub('signal|upstream|downstream', '', term),
                       levels = c('A', 'G', 'C', 'U', 'AAUAAA', 'AUUAAA')),
         color = factor(case_when(estimate > 0 & p.value < .05 ~ 1,
                                  estimate < 0 & p.value < .05 ~ 2,
                                  TRUE ~ 3)))

ng_plt <- ggplot(ng_df_full , aes(x = term, y = estimate)) +
  geom_segment(aes(xend = term), yend = 0 ) +
  geom_point(size = 2, aes(color = color)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  facet_grid( . ~ dataset, scales = "free_x", space = "free") +
  rremove('xlab') +
  rremove('legend') +
  ylim(-1,3) +
  ggtitle('N. gruberi') +
  theme(title = element_text(size = 14)) +
  scale_color_manual(values = c('blue', 'orange', 'gray')) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
tf_plt <- ggplot(tf_df_full , aes(x = term, y = estimate)) +
  geom_segment(aes(xend = term), yend = 0 ) +
  geom_point(size = 2, aes(color = color)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  facet_grid( . ~ dataset, scales = "free_x", space = "free") +
  rremove('xlab') +
  rremove('legend') +
  ylim(-1,3) +
  ggtitle('T. foetus') +
  theme(title = element_text(size = 14)) +
  scale_color_manual(values = c('blue', 'orange', 'gray')) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
tv_plt <- ggplot(tv_df_full , aes(x = term, y = estimate)) +
  geom_segment(aes(xend = term), yend = 0 ) +
  geom_point(size = 2, aes(color = color)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  facet_grid( . ~ dataset, scales = "free_x", space = "free") +
  rremove('xlab') +
  rremove('legend') +
  ylim(-1,3) +
  ggtitle('T. vaginalis') +
  theme(title = element_text(size = 14)) +
  scale_color_manual(values = c('blue', 'orange', 'gray')) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
pdf(here('figs/6B_raw.pdf'), width = 6, height = 12) 
ggarrange(ng_plt, tf_plt, tv_plt, nrow =3)
dev.off()
