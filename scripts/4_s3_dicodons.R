library(tidyverse)
library(magrittr)
library(Biostrings)
library(ggpubr)
library(here)

protist_theme <- theme_pubr() +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
theme_set(protist_theme)

glamb_AA <- readAAStringSet(here('./data/sequences/glamb_AA.fa.gz'))
glamb_DNA <- readDNAStringSet(here('./data/sequences/glamb_cds_DNA.fa.gz')) %>%
  .[names(.) %in% names(glamb_AA)]
gsm_AA <- readAAStringSet(here('./data_sequences/gsm_AA.fa.gz'))
gsm_DNA <- readDNAStringSet(here('./data/sequences/gsm_cds_DNA.fa.gz')) %>%
  .[names(.) %in% names(gsm_AA)]
gmur_AA <- readAAStringSet(here('./data/sequences/gmur_AA.fa.gz'))
gmur_DNA <- readDNAStringSet(here('./data/sequences/gmur_cds_DNA.fa.gz')) %>%
  .[names(.) %in% names(gmur_AA)]
trich_AA <- readAAStringSet(here('./data/sequences/trich_AA.fa.gz'))
trich_DNA <- readDNAStringSet(here('./data/sequences/trich_cds_DNA.fa.gz')) %>%
  .[names(.) %in% names(trich_AA)]
trif_AA <- readAAStringSet(here('./data/sequences/trif_AA.fa.gz'))
trif_DNA <- readDNAStringSet(here('./data/sequences/trif_cds_DNA.fa.gz')) %>%
  .[names(.) %in% names(trif_AA)]

gl_hex <- oligonucleotideFrequency(glamb_DNA, width = 6, step = 3) %>%
  t() %>%
  as.data.frame() %>%
  mutate(cnt = rowSums(.),
         dataset = 'G. lamblia') %>%
  dplyr::select(c('cnt', 'dataset')) %>%
  rownames_to_column('hex')  %>%
  mutate(hex = gsub('T', 'U', hex))
gs_hex <- oligonucleotideFrequency(gsm_DNA, width = 6, step = 3) %>%
  t() %>%
  as.data.frame() %>%
  mutate(cnt = rowSums(.),
         dataset = 'G. lamblia B-GS') %>%
  dplyr::select(c('cnt', 'dataset')) %>%
  rownames_to_column('hex')  %>%
  mutate(hex = gsub('T', 'U', hex))
gm_hex <- oligonucleotideFrequency(gmur_DNA, width = 6, step = 3) %>%
  t() %>%
  as.data.frame() %>%
  mutate(cnt = rowSums(.),
         dataset = 'G. muris') %>%
  dplyr::select(c('cnt', 'dataset')) %>%
  rownames_to_column('hex') %>%
  mutate(hex = gsub('T', 'U', hex))
tv_hex <- oligonucleotideFrequency(trich_DNA, width = 6, step = 3) %>%
  t() %>%
  as.data.frame() %>%
  mutate(cnt = rowSums(.),
         dataset = 'T. vaginalis') %>%
  dplyr::select(c('cnt', 'dataset')) %>%
  rownames_to_column('hex') %>%
  mutate(hex = gsub('T', 'U', hex))
tf_hex <- oligonucleotideFrequency(trif_DNA, width = 6, step = 3) %>%
  t() %>%
  as.data.frame() %>%
  mutate(cnt = rowSums(.),
         dataset = 'T. foetus') %>%
  dplyr::select(c('cnt', 'dataset')) %>%
  rownames_to_column('hex') %>%
  mutate(hex = gsub('T', 'U', hex))
hex_df <- rbind(gl_hex, gs_hex, gm_hex, tv_hex, tf_hex) %>%
  mutate(dataset = factor(dataset, 
                          levels = c('G. lamblia', 'G. lamblia B-GS', 'G. muris', 
                                     'T. foetus', 'T. vaginalis')))

CK_df <-  hex_df %>%
  filter(hex %in% c('UGUAAA','UGUAAG','UGCAAA','UGCAAG')) %>%
  mutate(codon = factor(if_else(hex == 'UGUAAA', 'UGU AAA', 'other CK codons'), 
                        levels = c('UGU AAA', 'other CK codons')))
CK_df_summary <- CK_df %>% 
  group_by(dataset, codon) %>%
  summarise(
    sd = sd(cnt),
    cnt = mean(cnt)
  ) %>%
  replace(is.na(.), 0)

pdf(here('./figs/CK_codons_bar.pdf'), width = 7.5, height = 3)
ggplot(CK_df, aes(codon, cnt, color = codon)) +
  geom_bar(stat = "identity", data = CK_df_summary,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2)) +
  scale_color_manual(values = c('red', 'black')) +
  facet_wrap(vars(dataset), nrow = 1, , scales="free_y") +
  theme_pubr() +
  ylab('counts') +
  rremove('xlab') +
  rremove('legend') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

KC_df <-  hex_df %>%
  filter(hex %in% c('AAAUGU','AAGUGU','AAAUGC','AAGUGC')) %>%
  mutate(codon = factor(if_else(hex == 'AAAUGU', 'AAA UGU', 'other KC codons'), 
                        levels = c('AAA UGU', 'other KC codons')))
KC_df_summary <- KC_df %>% 
  group_by(dataset, codon) %>%
  summarise(
    sd = sd(cnt),
    cnt = mean(cnt)
  ) %>%
  replace(is.na(.), 0)

pdf(here('./figs/KC_codons_bar.pdf'), width = 7.5, height = 3)
ggplot(KC_df, aes(codon, cnt, color = codon)) +
  geom_bar(stat = "identity", data = KC_df_summary,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2)) +
  scale_color_manual(values = c('red', 'black')) +
  facet_wrap(vars(dataset), nrow = 1, , scales="free_y") +
  theme_pubr() +
  ylab('counts') +
  rremove('xlab') +
  rremove('legend') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


SK_df <-  hex_df %>%
  filter(hex %in% c('AGCAAA', 'AGUAAA','UCAAAA','UCCAAA','UCGAAA','UCUAAA',
                    'AGCAAG', 'AGUAAG','UCAAAG','UCCAAG','UCGAAG','UCUAAG')) %>%
  mutate(codon = factor(if_else(hex == 'AGUAAA', 'AGU AAA', 'other SK codons'), 
                        levels = c('AGU AAA', 'other SK codons')))
SK_df_summary <- SK_df %>% 
  group_by(dataset, codon) %>%
  summarise(
    sd = sd(cnt),
    cnt = mean(cnt)
  ) %>%
  replace(is.na(.), 0)

pdf(here('./figs/SK_codons_bar.pdf'), width = 7.5, height = 3)
ggplot(SK_df, aes(codon, cnt, color = codon)) +
  geom_bar(stat = "identity", data = SK_df_summary,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2)) +
  scale_color_manual(values = c('red', 'black')) +
  facet_wrap(vars(dataset), nrow = 1, , scales="free_y") +
  theme_pubr() +
  ylab('counts') +
  rremove('xlab') +
  rremove('legend') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

KS_df <-  hex_df %>%
  filter(hex %in% c('AAAAGC', 'AAAAGU','AAAUCA','AAAUCC','AAAUCG','AAAUCU',
                    'AAGAGC', 'AAGAGU','AAGUCA','AAGUCC','AAGUCG','AAGUCU')) %>%
  mutate(codon = factor(if_else(hex == 'AAAAGU', 'AAA AGU', 'other KS codons'), 
                        levels = c('AAA AGU', 'other KS codons')))
KS_df_summary <- KS_df %>% 
  group_by(dataset, codon) %>%
  summarise(
    sd = sd(cnt),
    cnt = mean(cnt)
  ) %>%
  replace(is.na(.), 0)

pdf(here('./figs/KS_codons_bar.pdf'), width = 7.5, height = 3)
ggplot(KS_df, aes(codon, cnt, color = codon)) +
  geom_bar(stat = "identity", data = KS_df_summary,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2)) +
  scale_color_manual(values = c('red', 'black')) +
  facet_wrap(vars(dataset), nrow = 1, , scales="free_y") +
  theme_pubr() +
  ylab('counts') +
  rremove('xlab') +
  rremove('legend') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

SE_df <-  hex_df %>%
  filter(hex %in% c('AGCGAA', 'AGUGAA','UCAGAA','UCCGAA','UCGGAA','UCUGAA',
                    'AGCAAG', 'AGUAAG','UCAGAG','UCCGAG','UCGGAG','UCUGAG')) %>%
  mutate(codon = factor(if_else(hex == 'AGUGAA', 'AGU GAA', 'other SE codons'), 
                        levels = c('AGU GAA', 'other SE codons')))
SE_df_summary <- SE_df %>% 
  group_by(dataset, codon) %>%
  summarise(
    sd = sd(cnt),
    cnt = mean(cnt)
  ) %>%
  replace(is.na(.), 0)

pdf(here('./figs/SE_codons_bar.pdf'), width = 7.5, height = 3)
ggplot(SE_df, aes(codon, cnt, color = codon)) +
  geom_bar(stat = "identity", data = SE_df_summary,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2)) +
  scale_color_manual(values = c('red', 'black')) +
  facet_wrap(vars(dataset), nrow = 1, , scales="free_y") +
  theme_pubr() +
  ylab('counts') +
  rremove('xlab') +
  rremove('legend') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

ES_df <-  hex_df %>%
  filter(hex %in% c('GAAAGC', 'GAAAGU','GAAUCA','GAAUCC','GAAUCG','GAAUCU',
                    'AAGAGC', 'AAGAGU','GAGUCA','GAGUCC','GAGUCG','GAGUCU')) %>%
  mutate(codon = factor(if_else(hex == 'GAAAGU', 'GAA AGU', 'other ES codons'), 
                        levels = c('GAA AGU', 'other ES codons')))
ES_df_summary <- ES_df %>% 
  group_by(dataset, codon) %>%
  summarise(
    sd = sd(cnt),
    cnt = mean(cnt)
  ) %>%
  replace(is.na(.), 0)

pdf(here('./figs/ES_codons_bar.pdf'), width = 7.5, height = 3)
ggplot(ES_df, aes(codon, cnt, color = codon)) +
  geom_bar(stat = "identity", data = ES_df_summary,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2)) +
  scale_color_manual(values = c('red', 'black')) +
  facet_wrap(vars(dataset), nrow = 1, , scales="free_y") +
  theme_pubr() +
  ylab('counts') +
  rremove('xlab') +
  rremove('legend') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
