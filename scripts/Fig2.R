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


###2B###
#N. gruberi
ngNT <- readDNAStringSet(here('./data/sequences/neg_cds_ext.fa.gz')) %>%
  .[width(.) > 9]
ngNT_end <- ngNT %>%
  subseq(., start = width(.) -9, end = width(.))
ng_dual_A <- ngNT_end %>%
  .[grepl('AATAAA',.)]
ng_dual_G1 <- ngNT_end %>%
  .[grepl('AGT[A|G]AA',.)]
ng_dual_G2 <- ngNT_end %>%
  .[grepl('TGTAAA',.)]
ng_dual_G <- c(ng_dual_G1, ng_dual_G2)

#T. vaginalis
tvNT <- readDNAStringSet(here('./data/sequences/trich_cds_ext.fa.gz')) 
tvNT_end <- tvNT %>%
  subseq(., start = width(.) -9, end = width(.))
tv_dual_A <- tvNT_end %>%
  .[grepl('AATAAA',.)]
tv_dual_G1 <- tvNT_end %>%
  .[grepl('AGT[A|G]AA',.)]
tv_dual_G2 <- tvNT_end %>%
  .[grepl('TGTAAA',.)]
tv_dual_G <- c(tv_dual_G1, tv_dual_G2)

#T. foetus
tfNT <- readDNAStringSet(here('./data/sequences/trif_cds_ext.fa.gz')) 
tfNT_end <- tfNT %>%
  subseq(., start = width(.) -9, end = width(.))
tf_dual_A <- tfNT_end %>%
  .[grepl('AATAAA',.)]
tf_dual_G1 <- tfNT_end %>%
  .[grepl('AGT[A|G]AA',.)]
tf_dual_G2 <- tfNT_end %>%
  .[grepl('TGTAAA',.)]
tf_dual_G <- c(tf_dual_G1, tf_dual_G2)

#G. lamblia A
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

#G. lamblia B
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

#G. muris
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
                 counts = c((length(ngNT)-length(ng_dual_A)-length(ng_dual_G)),
                            length(ng_dual_A), length(ng_dual_G),
                            (length(tvNT)-length(tv_dual_A)-length(tv_dual_G)),
                            length(tv_dual_A), length(tv_dual_G),
                            (length(tfNT)-length(tf_dual_A)-length(tf_dual_G)),
                            length(tf_dual_A), length(tf_dual_G),
                            (length(glNT)-length(gl_dual_A)-length(gl_dual_G)),
                            length(gl_dual_A), length(gl_dual_G),
                            (length(gsNT)-length(gs_dual_A)-length(gs_dual_G)),
                            length(gs_dual_A), length(gs_dual_G),
                            (length(gmNT)-length(gm_dual_A)-length(gm_dual_G)),
                            length(gm_dual_A), length(gm_dual_G)),
                 dataset = rep(c('N. gruberi', 'T. vaginalis', 'T. foetus', 'G. lamblia A', 'G. lamblia B', 'G. muris'), each = 3)) %>%
  group_by(dataset) %>%
  mutate(freq = counts/sum(counts))

pdf(here('figs/2B_raw.pdf'), width = 4, height = 3)
ggbarplot(du, x = 'dataset', y = 'freq', fill = 'du', ylab = 'fraction',
          palette = c('gray', '#00a875ff', '#0071bcff'))+
  rremove('legend.title') +
  rremove('xlab') +
  theme(legend.position="right") +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


###2C###
#N. gruberi
ng_utr3 <- readDNAStringSet(here('./data/sequences/neg_utr3_ext.fa.gz'))
ng_utr3_start <- ng_utr3 %>%
  subseq(., start = 1, end = ifelse(width(.) > 10, 10, width(.)))
ng_utr3_dualA <- ng_utr3_start %>%
  .[grepl('AATAAA',.)]
ng_utr3_dualG1 <- ng_utr3_start %>%
  .[grepl('AGT[A|G]AA',.)]
ng_utr3_dualG2 <- ng_utr3_start %>%
  .[grepl('TGTAAA',.)]
ng_utr3_dualG <- c(ng_utr3_dualG1, ng_utr3_dualG2)
ng_df_A <- data.frame(id = names(ng_utr3_dualA)) %>%
  mutate(utr_len = width(ng_utr3)[match(id, names(ng_utr3))]-5)
ng_df_G <- data.frame(id = names(ng_utr3_dualG)) %>%
  mutate(utr_len = width(ng_utr3)[match(id, names(ng_utr3))]-5)

len_df_ng <- data.frame(id = names(ng_utr3), 
                        length = width(ng_utr3)-5,
                        dataset = 'N. gruberi') %>%
  mutate(du1 = factor(case_when(id %in% ng_df_A$id ~ 'AAUAAA',
                                id %in% ng_df_G$id ~ 'WGURAA',
                                TRUE ~ 'no dual usage'), 
                      levels = c('no dual usage', 'WGURAA', 'AAUAAA')),
         du2 = factor(ifelse(id %in% ng_df_A$id | id %in% ng_df_G$id, 'dual usage', 'no dual usage')),
         du3 = factor(ifelse(id %in% ng_df_G$id, 'WGURAA', 'no dual usage'))
  )

#G. lamblia A
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
gl_df_G <- data.frame(id = names(gl_utr3_dualG)) %>%
  mutate(utr_len = width(gl_utr3)[match(id, names(gl_utr3))]-5)

len_df_gl <- data.frame(id = names(gl_utr3), 
                        length = width(gl_utr3)-5,
                        dataset = 'G. lamblia A') %>%
  mutate(du1 = factor(case_when(id %in% gl_df_A$id ~ 'AAUAAA',
                                id %in% gl_df_G$id ~ 'WGURAA',
                                TRUE ~ 'no dual usage'), 
                      levels = c('no dual usage', 'WGURAA', 'AAUAAA')),
         du2 = factor(ifelse(id %in% gl_df_A$id | id %in% gl_df_G$id, 'dual usage', 'no dual usage')),
         du3 = factor(ifelse(id %in% gl_df_G$id, 'WGURAA', 'no dual usage'))
         )

len_df <- rbind(len_df_ng, len_df_gl) %>%
  mutate(dataset = factor(dataset, levels = c("N. gruberi", "G. lamblia A")))

pdf(here('figs/2C_raw.pdf'), width = 6, height = 4)
ggviolin(len_df, x = 'du1', y = 'length', 
         facet.by = 'dataset', color = 'gray', fill = 'du1', palette = c('gray','#0071bcff', '#00a875ff')) +
  geom_boxplot(width = .05, outlier.size = .001) +
  rremove('xlab') +
  theme(text = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1)) +
  scale_y_log10(limits = c(1,1000)) +
  rremove("legend")
dev.off()
