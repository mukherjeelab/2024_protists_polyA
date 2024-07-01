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


###1C###
ent <- read_tsv(here('./data/pA_signals/ent_pA.tsv.gz')) %>%
  mutate(dataset = 'E. histolitica')
neg <- read_tsv(here('./data/pA_signals/neg_pA.tsv.gz')) %>%
  mutate(dataset = 'N. gruberi')
trif <- read_tsv(here('./data/pA_signals/trif_pA.tsv.gz')) %>%
  mutate(dataset = 'T. foetus')
trich <- read_tsv(here('./data/pA_signals/trich_pA.tsv.gz')) %>%
  mutate(dataset = 'T. vaginalis')
glamb <- read_tsv(here('./data/pA_signals/glamb_pA.tsv.gz')) %>%
  mutate(dataset = 'G. lamblia')
gsm <- read_tsv(here('./data/pA_signals/gsm_pA.tsv.gz')) %>%
  mutate(dataset = 'G. lamblia B-GS')
gmur <- read_tsv(here('./data/pA_signals/gmur_pA.tsv.gz')) %>%
  mutate(dataset = 'G. muris')
hs <- read_tsv('./data/pA_signals/H_sapiens_signals.tsv.gz') %>%
  set_colnames(colnames(ent))
pa_df <- rbind(hs, ent, neg, trif, trich, glamb, gsm, gmur) %>%
  mutate(color = factor(case_when(polyA_signals %in% c('UGUAAA', 'AGUAAA', 'AGUGAA') ~ 'red',
                                  polyA_signals %in% c('AUUAAA', 'AAUAAA') ~ 'blue',
                                  TRUE ~ 'black'), levels = c('red', 'blue', 'black')),
         polyA_signals = factor(polyA_signals, levels = c('AGUAAA', 'UGUAAA', 'AGUGAA', 'AAUUGA', 'AAUUAA', 'UUAAUU', 'AUUAAA', 'AAUAAA')),
         dataset = factor(dataset, levels = c('H. sapiens','E. histolitica', 'N. gruberi', 'T. foetus', 'T. vaginalis', 'G. lamblia', 'G. lamblia B-GS', 'G. muris')))

pdf(here('./figs/1C_raw.pdf'), width = 14, height = 5)
ggbarplot(pa_df, x = 'polyA_signals', y = 'frequency', color = 'black', fill = 'color',
          palette = c('#0071bcff', '#00a875ff', 'black'), orientation = 'horizontal', facet.by = 'dataset',
          nrow = 1) +
  rremove('ylab') +
  rremove('legend') +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 14))
dev.off()


###1D###
eh <- readDNAStringSet(here('./data/sequences/ent_3UTR_combined.fa.gz'))
eh_motif1 <- eh %>%
  .[grepl('TTAATT', .)] %>%
  vmatchPattern('TTAATT',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(eh)[match(names, names(eh))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
eh_motif2 <- eh %>%
  .[grepl('AATTAA', .)] %>%
  vmatchPattern('AATTAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(eh)[match(names, names(eh))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
eh_motif3 <- eh %>%
  .[grepl('AATTGA', .)] %>%
  vmatchPattern('AATTGA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(eh)[match(names, names(eh))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
eh_motif <- rbind(eh_motif1, eh_motif2, eh_motif3) %>%
  mutate(dataset = 'E. histolitica')

ng <- readDNAStringSet(here('./data/sequences/neg_3UTR_combined.fa.gz'))
ng_motif1 <- ng %>%
  .[grepl('AATAAA', .)] %>%
  vmatchPattern('AATAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(ng)[match(names, names(ng))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
ng_motif2 <- ng %>%
  .[grepl('ATTAAA', .)] %>%
  vmatchPattern('ATTAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(ng)[match(names, names(ng))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
ng_motif <- rbind(ng_motif1, ng_motif2) %>%
  mutate(dataset = 'N. gruberi')

tf <- readDNAStringSet(here('./data/sequences/trif_3UTR_combined.fa.gz'))
tf_motif1 <- tf %>%
  .[grepl('AATAAA', .)] %>%
  vmatchPattern('AATAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(tf)[match(names, names(tf))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
tf_motif2 <- tf %>%
  .[grepl('ATTAAA', .)] %>%
  vmatchPattern('ATTAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(tf)[match(names, names(tf))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
tf_motif <- rbind(tf_motif1, tf_motif2) %>%
  mutate(dataset = 'T. foetus')


tv <- readDNAStringSet(here('./data/sequences/trich_3UTR_combined.fa.gz'))
tv_motif1 <- tv %>%
  .[grepl('AATAAA', .)] %>%
  vmatchPattern('AATAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(tv)[match(names, names(tv))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
tv_motif2 <- tv %>%
  .[grepl('ATTAAA', .)] %>%
  vmatchPattern('ATTAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(tv)[match(names, names(tv))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
tv_motif <- rbind(tv_motif1, tv_motif2) %>%
  mutate(dataset = 'T. vaginalis')

gl <- readDNAStringSet(here('./data/sequences/glamb_3UTR_combined.fa.gz'))
gl_motif1 <- gl %>%
  .[grepl('AGTAAA', .)] %>%
  vmatchPattern('AGTAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(gl)[match(names, names(gl))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
gl_motif2 <- gl %>%
  .[grepl('AGTGAA', .)] %>%
  vmatchPattern('AGTGAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(gl)[match(names, names(gl))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
gl_motif3 <- gl %>%
  .[grepl('TGTAAA', .)] %>%
  vmatchPattern('TGTAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(gl)[match(names, names(gl))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
gl_motif <- rbind(gl_motif1, gl_motif2, gl_motif3) %>%
  mutate(dataset = 'G. lamblia')

gm <- readDNAStringSet(here('./data/sequences/gmur_3UTR_combined.fa.gz'))
gm_motif1 <- gm %>%
  .[grepl('AGTAAA', .)] %>%
  vmatchPattern('AGTAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(gm)[match(names, names(gm))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
gm_motif2 <- gm %>%
  .[grepl('AGTGAA', .)] %>%
  vmatchPattern('AGTGAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(gm)[match(names, names(gm))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
gm_motif3 <- gm %>%
  .[grepl('TGTAAA', .)] %>%
  vmatchPattern('TGTAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(gm)[match(names, names(gm))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
gm_motif <- rbind(gm_motif1, gm_motif2, gm_motif3) %>%
  mutate(dataset = 'G. muris')

gs <- readDNAStringSet(here('./data/sequences/gsm_3UTR_combined.fa.gz'))
gs_motif1 <- gs %>%
  .[grepl('AGTAAA', .)] %>%
  vmatchPattern('AGTAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(gs)[match(names, names(gs))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
gs_motif2 <- gs %>%
  .[grepl('AGTGAA', .)] %>%
  vmatchPattern('AGTGAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(gs)[match(names, names(gs))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
gs_motif3 <- gs %>%
  .[grepl('TGTGAA', .)] %>%
  vmatchPattern('TGTGAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(gs)[match(names, names(gs))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
gs_motif <- rbind(gs_motif1, gs_motif2, gs_motif3) %>%
  mutate(dataset = 'G. lamblia B-GS')

hs <- readDNAStringSet(here('./data/sequences/hs_upstream40.fa.gz'))
hs_motif1 <- hs %>%
  .[grepl('AATAAA', .)] %>%
  vmatchPattern('AATAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(hs)[match(names, names(hs))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
hs_motif2 <- hs %>%
  .[grepl('ATTAAA', .)] %>%
  vmatchPattern('ATTAAA',.) %>%
  unlist(., use.names = TRUE) %>%
  as.data.frame() %>%
  mutate(width = width(hs)[match(names, names(hs))],
         strand = ifelse(grepl('\\+', names), '+', '-'),
         dist = ifelse(strand == '+', width-end, start)) %>%
  group_by(names) %>%
  summarise(min_dist = min(dist)) %>%
  filter(min_dist <= 100)
hs_motif <- rbind(hs_motif1, hs_motif2) %>%
  mutate(dataset = 'H. sapiens')

motif_dist <- rbind(hs_motif, eh_motif, ng_motif, tf_motif, tv_motif, gl_motif, gs_motif, gm_motif) %>%
  mutate(dataset = factor(dataset, levels = c('H. sapiens', 'E. histolitica', 'N. gruberi', 'T. foetus', 'T. vaginalis', 'G. lamblia', 'G. lamblia B-GS', 'G. muris')),
         signal = factor(case_when(dataset %in% c('G. lamblia', 'G. lamblia B-GS', 'G. muris') ~ 'AGURAA', 
                                   dataset %in% c('H. sapiens', 'N. gruberi', 'T. vaginalis', 'T. foetus') ~ 'AWUAAA',
                                   TRUE ~ 'other'), levels = c('AWUAAA', 'AGURAA', 'other')),
         min_dist = -min_dist)

pdf(here('./figs/1D_raw.pdf'), width = 12, height = 3)
ggdensity(motif_dist, x = 'min_dist', color = 'signal', palette = c('#00a875ff','#0071bcff', 'black'),
          facet.by = 'dataset', nrow = 1, xlab = 'distance from cleavage site') +
  xlim(c(-40,0)) +
  theme(legend.position="right") +
  theme(legend.title = element_text(size = 14)) +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14))
dev.off()
