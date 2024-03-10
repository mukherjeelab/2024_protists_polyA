library(tidyverse)
library(magrittr)
library(matrixStats)
library(Biostrings)
library(ggpubr)
library(here)

protist_theme <- theme_pubr() +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
theme_set(protist_theme)

glNT <- readDNAStringSet(here('./data/sequences/glamb_cds_ext.fa.gz')) 
glNT_end <- glNT %>%
  subseq(., start = width(.) -9, end = width(.))
gl_dual_A <- glNT_end %>%
  .[grepl('AATAAA',.)]
gl_dual_G1 <- glNT_end %>%
  .[grepl('AGT[A|G]AA',.)]
gl_dual_G2 <- glNT_end %>%
  .[grepl('TGTAAA',.)]
gl_dual_G <- c(gl_dual_G1, gl_dual_G2)
gsNT <- readDNAStringSet(here('./data/sequences/gsm_cds_ext.fa.gz')) 
gsNT_end <- gsNT %>%
  subseq(., start = width(.) -9, end = width(.))
gs_dual_A <- gsNT_end %>%
  .[grepl('AATAAA',.)]
gs_dual_G1 <- gsNT_end %>%
  .[grepl('AGT[A|G]AA',.)]
gs_dual_G2 <- gsNT_end %>%
  .[grepl('TGTAAA',.)]
gs_dual_G <- c(gs_dual_G1, gs_dual_G2)
gmNT <- readDNAStringSet(here('./data/sequences/gmur_cds_ext.fa.gz')) %>%
  .[width(.) > 9]
gmNT_end <- gmNT %>%
  subseq(., start = width(.) -9, end = width(.))
gm_dual_A <- gmNT_end %>%
  .[grepl('AATAAA',.)]
gm_dual_G1 <- gmNT_end %>%
  .[grepl('AGT[A|G]AA',.)]
gm_dual_G2 <- gmNT_end %>%
  .[grepl('TGTAAA',.)]
gm_dual_G <- c(gm_dual_G1, gm_dual_G2)

du <- data.frame(du = factor(rep(c('no dual', 'AAUAAA', 'WGURAA'), 3), 
                             levels = c('no dual', 'AAUAAA', 'WGURAA')),
                 counts = c((length(glNT)-length(gl_dual_A)-length(gl_dual_G)),
                            length(gl_dual_A), length(gl_dual_G),
                            (length(gsNT)-length(gs_dual_A)-length(gs_dual_G)),
                            length(gs_dual_A), length(gs_dual_G),
                            (length(gmNT)-length(gm_dual_A)-length(gm_dual_G)),
                            length(gm_dual_A), length(gm_dual_G)),
                 dataset = rep(c('G. lamblia', 'G. lamblia B-GS', 'G. muris'), each = 3)) %>%
  group_by(dataset) %>%
  mutate(freq = counts/sum(counts))

pdf(here('paper01/fig_elements/dual_occurence.pdf'), width = 4, height = 3)
ggbarplot(du, x = 'dataset', y = 'freq', fill = 'du', ylab = 'frequency',
          palette = c('gray', 'blue', 'red'))+
  rremove('legend.title') +
  rremove('xlab') +
  theme(legend.position="right") +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

gl_utr3 <- readDNAStringSet(here('./data/sequences/glamb_utr3_ext.fa.gz'))
gl_utr3_start <- gl_utr3 %>%
  subseq(., start = 1, end = ifelse(width(.) > 10, 10, width(.)))
gl_utr3_dualA <- gl_utr3_start %>%
   .[grepl('AATAAA',.)]
gl_utr3_dualG1 <- gl_utr3_start %>%
   .[grepl('AGT[A|G]AA',.)]
gl_utr3_dualG2 <- gl_utr3_start %>%
   .[grepl('TGTAAA',.)]
gl_utr3_dualG <- c(gl_utr3_dualG1, gl_utr3_dualG2)
gl_df_A <- data.frame(id = names(gl_utr3_dualA)) %>%
  mutate(utr_len = width(gl_utr3)[match(id, names(gl_utr3))]-5)
gl_df_A_close <- gl_df_A %>%
  filter(utr_len < 21)
gl_df_A_close_no_signal <- gl_utr3 %>%
  .[names(.) %in% gl_df_A_close$id] %>%
  subseq(., start = 6, end = width(.)) %>%
  .[!grepl('AATAAA|AGTAAA|AGTGAA|TGTAAA',.)]
gl_df_G <- data.frame(id = names(gl_utr3_dualG)) %>%
  mutate(utr_len = width(gl_utr3)[match(id, names(gl_utr3))]-5)
gl_df_G_close <- gl_df_G %>%
  filter(utr_len < 21)
gl_df_G_close_no_signal <- gl_utr3 %>%
  .[names(.) %in% gl_df_G_close$id] %>%
  subseq(., start = 6, end = width(.)) %>%
  .[!grepl('AATAAA|AGTAAA|AGTGAA|TGTAAA',.)] 
gs_utr3 <- readDNAStringSet(here('./data/sequences/gsm_utr3_ext.fa.gz'))
gs_utr3_start <- gs_utr3 %>%
  subseq(., start = 1, end = ifelse(width(.) > 10, 10, width(.)))
gs_utr3_dualA <- gs_utr3_start %>%
   .[grepl('AATAAA',.)]
gs_utr3_dualG1 <- gs_utr3_start %>%
   .[grepl('AGT[A|G]AA',.)]
gs_utr3_dualG2 <- gs_utr3_start %>%
   .[grepl('TGTAAA',.)]
gs_utr3_dualG <- c(gs_utr3_dualG1, gs_utr3_dualG2)
gs_df_A <- data.frame(id = names(gs_utr3_dualA)) %>%
  mutate(utr_len = width(gs_utr3)[match(id, names(gs_utr3))]-5)
gs_df_A_close <- gs_df_A %>%
  filter(utr_len < 21)
gs_df_A_close_no_signal <- gs_utr3 %>%
  .[names(.) %in% gs_df_A_close$id] %>%
  subseq(., start = 6, end = width(.)) %>%
  .[!grepl('AATAAA|AGTAAA|AGTGAA|TGTAAA',.)]
gs_df_G <- data.frame(id = names(gs_utr3_dualG)) %>%
  mutate(utr_len = width(gs_utr3)[match(id, names(gs_utr3))]-5)
gs_df_G_close <- gs_df_G %>%
  filter(utr_len < 21)
gs_df_G_close_no_signal <- gs_utr3 %>%
  .[names(.) %in% gs_df_G_close$id] %>%
  subseq(., start = 6, end = width(.)) %>%
  .[!grepl('AATAAA|AGTAAA|AGTGAA|TGTAAA',.)]
gm_utr3 <- readDNAStringSet(here('./data/sequences/gmur_utr3_ext.fa.gz'))
gm_utr3_start <- gm_utr3 %>%
  subseq(., start = 1, end = ifelse(width(.) > 10, 10, width(.)))
gm_utr3_dualA <- gm_utr3_start %>%
   .[grepl('AATAAA',.)]
gm_utr3_dualG1 <- gm_utr3_start %>%
   .[grepl('AGT[A|G]AA',.)]
gm_utr3_dualG2 <- gm_utr3_start %>%
   .[grepl('TGTAAA',.)]
gm_utr3_dualG <- c(gm_utr3_dualG1, gm_utr3_dualG2)
gm_df_A <- data.frame(id = names(gm_utr3_dualA)) %>%
  mutate(utr_len = width(gm_utr3)[match(id, names(gm_utr3))]-5)
gm_df_A_close <- gm_df_A %>%
  filter(utr_len < 21)
gm_df_A_close_no_signal <- gm_utr3 %>%
  .[names(.) %in% gm_df_A_close$id] %>%
  subseq(., start = 6, end = width(.)) %>%
  .[!grepl('AATAAA|AGTAAA|AGTGAA|TGTAAA',.)]
gm_df_G <- data.frame(id = names(gm_utr3_dualG)) %>%
  mutate(utr_len = width(gm_utr3)[match(id, names(gm_utr3))]-5)
gm_df_G_close <- gm_df_G %>%
  filter(utr_len < 21)
gm_df_G_close_no_signal <- gm_utr3 %>%
  .[names(.) %in% gm_df_G_close$id] %>%
  subseq(., start = 6, end = width(.)) %>%
  .[!grepl('AATAAA|AGTAAA|AGTGAA|TGTAAA',.)]
  
usage_df <- data.frame(dataset = rep(c('G. lamblia', 'G. lamblia B-GS', 'G. muris'), each = 4),
                       signal = rep(c('AAUAAA', 'WGURAA'), 6),
                       usage = rep(rep(c('not used', 'used'), 3), each = 2),
                       freq = c(((nrow(gl_df_A) - length(gl_df_A_close_no_signal))/nrow(gl_df_A)),
                                ((nrow(gl_df_G) - length(gl_df_G_close_no_signal))/nrow(gl_df_G)),
                                (length(gl_df_A_close_no_signal)/nrow(gl_df_A)),
                                (length(gl_df_G_close_no_signal)/nrow(gl_df_G)),
                                ((nrow(gs_df_A) - length(gs_df_A_close_no_signal))/nrow(gs_df_A)),
                                ((nrow(gs_df_G) - length(gs_df_G_close_no_signal))/nrow(gs_df_G)),
                                (length(gs_df_A_close_no_signal)/nrow(gs_df_A)),
                                (length(gs_df_G_close_no_signal)/nrow(gs_df_G)),
                                ((nrow(gm_df_A) - length(gm_df_A_close_no_signal))/nrow(gm_df_A)),
                                ((nrow(gm_df_G) - length(gm_df_G_close_no_signal))/nrow(gm_df_G)),
                                (length(gm_df_A_close_no_signal)/nrow(gm_df_A)),
                                (length(gm_df_G_close_no_signal)/nrow(gm_df_G)))
                       )

pdf(here('paper01/fig_elements/dual_usage.pdf'), width = 6, height = 3)
ggbarplot(usage_df, x = 'dataset', y = 'freq', fill = 'usage',
          ylab = 'frequency', palette = c('gray', 'red'), facet.by = 'signal') +
  rremove('legend.title') +
  rremove('xlab') +
  theme(legend.position="right") +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

usage_df_G_counts <- data.frame(dataset = rep(c('G. lamblia', 'G. lamblia B-GS', 'G. muris'),
                                              each = 2),
                       usage = factor(rep(c('used','not used'), 3), 
                                      levels = c('used', 'not used')),
                       counts = c(length(gl_df_G_close_no_signal),
                                  nrow(gl_df_G)-length(gl_df_G_close_no_signal),
                                  length(gs_df_G_close_no_signal),
                                  nrow(gs_df_G)-length(gs_df_G_close_no_signal),
                                  length(gm_df_G_close_no_signal),
                                  nrow(gm_df_G)-length(gm_df_G_close_no_signal))
)
pdf(here('paper01/fig_elements/2c_new_dual_usage_counts.pdf'), width = 4, height = 3)
ggbarplot(usage_df_G_counts, x = 'usage', y = 'counts', fill = 'usage',
          palette = c('red', 'gray'), facet.by = 'dataset', main = 'WGURAA') +
  rremove('legend') +
  rremove('xlab') +
  theme(text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

usage_df_A_counts <- data.frame(dataset = rep(c('G. lamblia', 'G. lamblia B-GS', 'G. muris'),
                                              each = 2),
                       usage = factor(rep(c('used','not used'), 3), 
                                      levels = c('used', 'not used')),
                       counts = c(length(gl_df_A_close_no_signal),
                                  nrow(gl_df_A)-length(gl_df_A_close_no_signal),
                                  length(gs_df_A_close_no_signal),
                                  nrow(gs_df_A)-length(gs_df_A_close_no_signal),
                                  length(gm_df_A_close_no_signal),
                                  nrow(gm_df_A)-length(gm_df_A_close_no_signal))
)
pdf(here('paper01/fig_elements/s2a_dual_usage_counts_AAUAAA.pdf'), width = 4, height = 3)
ggbarplot(usage_df_A_counts, x = 'usage', y = 'counts', fill = 'usage',
          palette = c('red', 'gray'), facet.by = 'dataset', main = 'AAUAAA') +
  rremove('legend') +
  rremove('xlab') +
  theme(text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

len_df_gl <- data.frame(id = names(gl_utr3), 
                        length = width(gl_utr3)-5,
                        dataset = 'G. lamblia') %>%
  mutate(du1 = factor(case_when(id %in% gl_df_A$id ~ 'AAUAAA',
                                id %in% gl_df_G$id ~ 'WGURAA',
                                TRUE ~ 'no dual usage'), 
                      levels = c('no dual usage', 'WGURAA', 'AAUAAA')),
         du2 = factor(ifelse(id %in% gl_df_A$id | id %in% gl_df_G$id, 'dual usage', 'no dual usage')),
         du3 = factor(ifelse(id %in% gl_df_G$id, 'WGURAA', 'no dual usage'))
         )
len_df_gs <- data.frame(id = names(gs_utr3) , 
                        length = width(gs_utr3)-5,
                        dataset = 'G. lamblia B-GS') %>%
  mutate(du1 = factor(case_when(id %in% gs_df_A$id ~ 'AAUAAA',
                                id %in% gs_df_G$id ~ 'WGURAA',
                                TRUE ~ 'no dual usage'), 
                      levels = c('no dual usage', 'WGURAA', 'AAUAAA')),
         du2 = factor(ifelse(id %in% gs_df_A$id | id %in% gs_df_G$id, 'dual usage', 'no dual usage')),
         du3 = factor(ifelse(id %in% gs_df_G$id, 'WGURAA', 'no dual usage'))
         )
len_df_gm <- data.frame(id = names(gm_utr3) , 
                        length = width(gm_utr3)-5,
                        dataset = 'G. muris') %>%
  mutate(du1 = factor(case_when(id %in% gm_df_A$id ~ 'AAUAAA',
                                id %in% gm_df_G$id ~ 'WGURAA',
                                TRUE ~ 'no dual usage'), 
                      levels = c('no dual usage', 'WGURAA', 'AAUAAA')),
         du2 = factor(ifelse(id %in% gm_df_A$id | id %in% gm_df_G$id, 'dual usage', 'no dual usage')),
         du3 = factor(ifelse(id %in% gm_df_G$id, 'WGURAA', 'no dual usage'))
         )

len_df <- rbind(len_df_gl, len_df_gs, len_df_gm) %>%
  mutate(dataset = factor(dataset))

pdf(here('paper01/fig_elements/2d_length.pdf'), width = 6, height = 4)
ggviolin(len_df, x = 'du1', y = 'length', 
         facet.by = 'dataset', color = 'gray', fill = 'gray') +
  geom_boxplot(width = .05, outlier.size = .001) +
  rremove('xlab') +
  theme(text = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1)) +
  scale_y_log10(limits = c(1,1000))
dev.off()
