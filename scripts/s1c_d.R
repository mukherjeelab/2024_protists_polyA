library(tidyverse)
library(magrittr)
library(ggpubr)
library(here)

protist_theme <- theme_pubr() +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
theme_set(protist_theme)

#S1C
library(tidyverse)
library(magrittr)
library(ggpubr)
library(here)

protist_theme <- theme_pubr() +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
theme_set(protist_theme)

#glamb
polyA_glamb1 <- read_tsv(here('data/nanopolish_res/glamb-1_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_glamb2 <- read_tsv(here('data/nanopolish_res/glamb-2_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_glamb3 <- read_tsv(here('data/nanopolish_res/glamb-3_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_glamb <- rbind(polyA_glamb1, polyA_glamb2, polyA_glamb3) %>%
  group_by(output) %>%
  summarise(cts = dplyr::n()) %>%
  mutate(percent = cts*100/colSums(.[,2])) %>%
  mutate(organism = 'G. lamblia A')

#gsm
polyA_gsm1 <- read_tsv(here('data/nanopolish_res/gsm-1_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_gsm2 <- read_tsv(here('data/nanopolish_res/gsm-2_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_gsm3 <- read_tsv(here('data/nanopolish_res/gsm-3_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_gsm <- rbind(polyA_gsm1, polyA_gsm2, polyA_gsm3) %>%
  group_by(output) %>%
  summarise(cts = dplyr::n()) %>%
  mutate(percent = cts*100/colSums(.[,2])) %>%
  mutate(organism = 'G. lamblia B')

#gmur
polyA_gmur1 <- read_tsv(here('data/nanopolish_res/gmur-1_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_gmur2 <- read_tsv(here('data/nanopolish_res/gmur-2_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_gmur <- rbind(polyA_gmur1, polyA_gmur2) %>%
  group_by(output) %>%
  summarise(cts = dplyr::n()) %>%
  mutate(percent = cts*100/colSums(.[,2])) %>%
  mutate(organism = 'G. muris')

#ent
polyA_ent1 <- read_tsv(here('data/nanopolish_res/ent-1_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_ent2 <- read_tsv(here('data/nanopolish_res/ent-2_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_ent3 <- read_tsv(here('data/nanopolish_res/ent-3_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_ent <- rbind(polyA_ent1, polyA_ent2, polyA_ent3) %>%
  group_by(output) %>%
  summarise(cts = dplyr::n()) %>%
  mutate(percent = cts*100/colSums(.[,2])) %>%
  mutate(organism = 'E. histolitica')

#neg
polyA_neg1 <- read_tsv(here('data/nanopolish_res/neg-1_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 14 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 15 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_neg2 <- read_tsv(here('data/nanopolish_res/neg-2_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 14 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 15 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_neg3 <- read_tsv(here('data/nanopolish_res/neg-3_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 14 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 15 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_neg <- rbind(polyA_neg1, polyA_neg2, polyA_neg3) %>%
  group_by(output) %>%
  summarise(cts = dplyr::n()) %>%
  mutate(percent = cts*100/colSums(.[,2])) %>%
  mutate(organism = 'N. gruberi')

#trich
polyA_trich1 <- read_tsv(here('data/nanopolish_res/trich-1_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_trich2 <- read_tsv(here('data/nanopolish_res/trich-2_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_trich3 <- read_tsv(here('data/nanopolish_res/trich-3_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_trich <- rbind(polyA_trich1, polyA_trich2, polyA_trich3) %>%
  group_by(output) %>%
  summarise(cts = dplyr::n()) %>%
  mutate(percent = cts*100/colSums(.[,2])) %>%
  mutate(organism = 'T. vaginalis')

#trif
polyA_trif1 <- read_tsv(here('data/nanopolish_res/trif-1_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_trif2 <- read_tsv(here('data/nanopolish_res/trif-2_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_trif3 <- read_tsv(here('data/nanopolish_res/trif-3_polyA.tsv.gz')) %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9,10)) %>%
  mutate(output = case_when(qc_tag == 'PASS' & polya_length > 29 ~ '% As included',
                            qc_tag == 'PASS' & polya_length < 30 ~ '% As too short',
                            TRUE ~ '% no As'))
polyA_trif <- rbind(polyA_trif1, polyA_trif2, polyA_trif3) %>%
  group_by(output) %>%
  summarise(cts = dplyr::n()) %>%
  mutate(percent = cts*100/colSums(.[,2])) %>%
  mutate(organism = 'T. foetus')


polyA_df <- rbind(polyA_glamb, polyA_gsm, polyA_gmur, polyA_neg, polyA_trich, polyA_trif, polyA_ent) %>%
  mutate(output = factor(output, 
                         levels = c('% As included', '% As too short', '% no As')),
         organism = factor(organism, 
                           levels = c('E. histolitica', 'N. gruberi', 'T. foetus', 'T. vaginalis', 'G. lamblia A', 'G. lamblia B', 'G. muris')))

pdf(here('figs/S1C_raw.pdf'))
ggbarplot(polyA_df, x = 'organism', y = 'percent', fill = 'output') +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14)) +
  rremove('xlab') +
  rremove('ylab') +
  rremove('legend.title') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#S1D
annot <- read_tsv(here('data/annotation/3utr_annot.tsv.gz')) %>%
  mutate(organism = factor(organism, 
                           levels = c('E. histolitica', 'N. gruberi', 'T. foetus', 'T. vaginalis', 'G. lamblia', 'G. lamblia B-GS', 'G. muris')),
         annotation = factor(annotation, levels = c('previous', 'current')))

pdf(here('figs/S1D_raw.pdf'))
ggbarplot(annot, x = 'organism', y = 'annotated_3utrs', fill = 'annotation', 
          position = position_dodge(0.9), ylab = "#annotated 3'UTRs",
          label = TRUE )  +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14)) +
  rremove('xlab') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

