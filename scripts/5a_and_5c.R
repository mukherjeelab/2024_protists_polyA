library(tidyverse)
library(magrittr)
library(ggpubr)
library(here)

protist_theme <- theme_pubr() +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
theme_set(protist_theme)

data_5a <- data.frame(element = factor(c('Full model', 'Full model with interactions', 
                                  'Poly(A) signal', 'Upstream nucleotide', 
                                  'Downstream nucleotide'), 
                                  levels = c('Downstream nucleotide', 'Upstream nucleotide', 
                                             'Poly(A) signal', 'Full model with interactions', 
                                             'Full model'
                                             )),
                      coefficient = c(.51,.85,.58,.26,.39))

pdf(here('./figs/gmur_R2_drop_5a_new.pdf'))
ggbarplot(data_5a, x = 'element', y = 'coefficient', orientation = 'horiz',
          fill = 'grey', color = 'grey', label = TRUE) +
  rremove('ylab')  +
  theme(text = element_text(size = 14)) + 
  geom_vline(xintercept = 4.5, linetype = 2)
dev.off()

data_5c <- read_tsv(here('./data/prem_cl/gmur_prem_cl_4_fig5c.tsv.gz')) %>%
  mutate(signal = factor(signal, levels = c('AGUAAA', 'UGUAAA', 'AGUGAA')),
         prem_cl = factor(prem_cl))

pdf(here('./figs/gmur_prem_cl_5c_new.pdf'))
ggbarplot(data_5c, x = 'signal', y = 'counts', 
          fill = 'prem_cl', palette = c('black', 'grey'),
          position = position_dodge()) +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  rremove('legend.title') +
  rremove('xlab') +
  ylim(0,80)
dev.off()
