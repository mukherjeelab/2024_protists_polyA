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

gl_df <- read_tsv(here('./data/scores_8mers/glamb_8mers.tsv.gz'))
gs_df <- read_tsv(here('./data/scores_8mers/gsm_8mers.tsv.gz'))
gm_df <- read_tsv(here('./data/scores_8mers/gmur_8mers.tsv.gz'))

gl_4_models <- gl_df %>%
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
gs_4_models <- gs_df %>%
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
model_full_gl <- lm(score ~signal+upstream+downstream, 
                 data = gl_4_models)
model_full_r_sq_gl <- model_full_gl %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_signal_gl <- lm(score ~upstream+downstream, 
                                 data = gl_4_models)
model_reduced_no_signal_r_sq_gl <- model_reduced_no_signal_gl %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_up_gl <- lm(score ~signal+downstream, 
                                 data = gl_4_models)
model_reduced_no_up_r_sq_gl <- model_reduced_no_up_gl %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_down_gl <- lm(score ~signal+upstream, 
                             data = gl_4_models)
model_reduced_no_down_r_sq_gl <- model_reduced_no_down_gl %>%
  glance() %>%
  dplyr::select('adj.r.squared')
r_sq_drop_gl <- tibble(term = c('full model', 'signal sequence', 'upstream nt', 'downstream nt'),
                       adj_rsq = c(as.numeric(model_full_r_sq_gl),
                                        as.numeric(model_reduced_no_signal_r_sq_gl),
                                        as.numeric(model_reduced_no_up_r_sq_gl),
                                        as.numeric(model_reduced_no_down_r_sq_gl)),
                       dataset = 'G. lamblia') %>%
  mutate(explained_variance = round(ifelse(term == 'full model', adj_rsq, 
                                     adj_rsq[1] -adj_rsq),digits = 2))
pdf(here('./figs/4D_raw.pdf'), width = 4, height = 5)
ggbarplot(r_sq_drop_gl, x = 'term', y = 'explained_variance', fill = 'gray',
          label = TRUE, lab.col = "black", lab.pos = "in", 
          lab.hjust = 0.5, lab.vjust = 1, lab.size = 5) +
  rremove('xlab') +
  ylab('adjusted R^2') +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.y = element_text(size = 14))
dev.off()

model_full_gs <- lm(score ~signal+upstream+downstream, 
                    data = gs_4_models)
model_full_r_sq_gs <- model_full_gs %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_signal_gs <- lm(score ~upstream+downstream, 
                                 data = gs_4_models)
model_reduced_no_signal_r_sq_gs <- model_reduced_no_signal_gs %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_up_gs <- lm(score ~signal+downstream, 
                             data = gs_4_models)
model_reduced_no_up_r_sq_gs <- model_reduced_no_up_gs %>%
  glance() %>%
  dplyr::select('adj.r.squared')
model_reduced_no_down_gs <- lm(score ~signal+upstream, 
                               data = gs_4_models)
model_reduced_no_down_r_sq_gs <- model_reduced_no_down_gs %>%
  glance() %>%
  dplyr::select('adj.r.squared')
r_sq_drop_gs <- tibble(term = c('full model', 'signal sequence', 'upstream nt', 'downstream nt'),
                       adj_rsq = c(as.numeric(model_full_r_sq_gs),
                                   as.numeric(model_reduced_no_signal_r_sq_gs),
                                   as.numeric(model_reduced_no_up_r_sq_gs),
                                   as.numeric(model_reduced_no_down_r_sq_gs)),
                       dataset = 'G. lamblia B-GS') %>%
  mutate(explained_variance = round(ifelse(term == 'full model', adj_rsq, 
                                           adj_rsq[1] -adj_rsq),digits = 2))
pdf(here('./figs/S3C_raw.pdf'), width = 4, height = 5)
ggbarplot(r_sq_drop_gs, x = 'term', y = 'explained_variance', fill = 'gray',
          label = TRUE, lab.col = "black", lab.pos = "in", 
          lab.hjust = 0.5, lab.vjust = 1, lab.size = 5) +
  rremove('xlab') +
  ylab('adjusted R^2') +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.y = element_text(size = 14))
dev.off()

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

gl_df_full <- model_full_gl %>%
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

pdf(here('./figs/4E_raw.pdf'), width = 6, height = 4) 
ggplot(gl_df_full , aes(x = term, y = estimate)) +
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

gs_df_full <- model_full_gs %>%
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

pdf(here('./figs/S3B_raw.pdf'), width = 6, height = 4) 
ggplot(gs_df_full , aes(x = term, y = estimate)) +
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
