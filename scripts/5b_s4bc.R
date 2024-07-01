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

gmur_df <- read_tsv(here('data/scores_8mers/gmur_8mers.tsv.gz')) %>%
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
fit <- lm(score ~ signal*upstream*downstream, data = gmur_df)

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