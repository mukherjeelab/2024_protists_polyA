library(tidyverse)
library(magrittr)
library(Biostrings)
library(broom)
library(rtracklayer)
library(plyranges)
library(data.table)
library(universalmotif)
library(ggmsa)
library(ggseqlogo)
library(ggcoverage)
library(ggpattern)
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
gsm_AA <- readAAStringSet(here('./data/sequences/gsm_AA.fa.gz'))
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

###4A###
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

pdf(here('./figs/4A_raw.pdf'), width = 7.5, height = 3)
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


###S3A###
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

pdf(here('./figs/S3A_1_raw.pdf'), width = 7.5, height = 3)
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

pdf(here('./figs/S3A_2_raw.pdf'), width = 7.5, height = 3)
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

pdf(here('./figs/S3A_3_raw.pdf'), width = 7.5, height = 3)
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

pdf(here('./figs/S3A_4_raw.pdf'), width = 7.5, height = 3)
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

pdf(here('./figs/S3A_5_raw.pdf'), width = 7.5, height = 3)
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


###4B###
cds_cts <- read_tsv(here("data/prem_cl/cds_counts.tsv.gz"))
pdf(here("figs/4B_raw.pdf"))
ggbarplot(cds_cts, x = 'dataset', y = 'counts', fill = 'gray') +
  rremove("xlab")
dev.off()


###4D###
gl_df <- read_tsv(here('./data/scores_8mers/glamb_8mers.tsv.gz'))
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


###S3C###
gs_df <- read_tsv(here('./data/scores_8mers/gsm_8mers.tsv.gz'))
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


###4E###
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


###S3B###
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


###4F###
gl_gtf <- import.gff(here('data/annotation/glamb_combined.gtf.gz'), format = 'gtf') %>%
  mutate(gene_name = 'GL50803_1890',
         gene_biotype = 'protein_coding') %>%
  plyranges::filter(strand == '+')
gl_bed <- gl_gtf %>%
  filter(gene_id == 'GL50803_1890' & type == 'transcript') %>%
  mutate(score = 1,
         start = start - 1)
peak_df_gl <- readDNAStringSet(here('data/sequences/GL50803_1890.fa.gz')) %>%
  vmatchPattern('AGTRAA', ., fixed = FALSE) %>%
  as.data.frame() %>%
  mutate(chr = 'GLCHR05',
         start = start(gl_bed) -1 + start,
         end = start(gl_bed) -1 + end) %>%
  dplyr::select(c('chr','start','end'))

track_df_glamb <- LoadTrackFile(track.file = here('data/sequences/GL50803_1890.bw'), 
                                format = "bw", extend = 0,
                                region = 'GLCHR05:1003565-1004000', 
                                gtf.gr = gl_gtf,
                                bin.size = 1,
                                gene.name = 'GL50803_1890')
pdf(here('figs/4F_raw.pdf'), height = 4, width = 6)
ggcoverage(track_df_glamb, color = 'gray') +
  geom_peak(peak.df = peak_df_gl) +  
  geom_gene(gtf.gr=gl_gtf, arrow.size = 0)
dev.off()


###S3F###
gs_gtf <- import.gff(here('data/annotation/gsm_nanopore.gtf.gz'), format = 'gtf') %>%
  mutate(gene_name = 'GL50581_3478',
         gene_biotype = 'protein_coding') %>%
  plyranges::filter(strand == '+')
gs_bed <- gs_gtf %>%
  filter(gene_id == 'GL50581_3478' & type == 'transcript') %>%
  mutate(score = 1,
         start = start - 1)
#export.bed(gs_bed, here('paper01/GL50581_3478.bed'))
#system('bedtools getfasta -fi ./new_genomes/gsm.fa -bed ./paper01/GL50581_3478.bed -s -fo ./paper01/GL50581_3478.fa')
peak_df_gs <- readDNAStringSet(here('data/sequences/GL50581_3478.fa.gz')) %>%
  vmatchPattern('WGTRAA', ., fixed = FALSE) %>%
  as.data.frame() %>%
  mutate(chr = 'ACGJ01002891',
         start = start(gs_bed) -1 + start,
         end = start(gs_bed) -1 + end) %>%
  dplyr::select(c('chr','start','end'))
track_df_gsm <- LoadTrackFile(track.file = here('data/sequences/GL50581_3478.bw'), 
                              format = "bw", extend = 0,
                              region = 'ACGJ01002891:25960-26820', 
                              gtf.gr = gs_gtf,
                              bin.size = 10,
                              gene.name = 'GL50581_3478')
pdf(here('figs/S3F_raw.pdf'), height = 4, width = 6)
ggcoverage(track_df_gsm, color = 'gray') +
  geom_peak(peak.df = peak_df_gs) +  
  geom_gene(gtf.gr=gs_gtf, arrow.size = 0)
dev.off()

###S3D###
motifs <- list.files(path = paste0(here(), '/data/predictive_kmers'), pattern = '.csv$',
                     full.names = TRUE)
cons_motifs <- c()
for (i in seq_along(motifs)) {
  tmp <- fread(motifs[i]) %>%
    as.matrix() %>%
    set_rownames(c('A','C','G','T'))
  tmp_mot <- c()
  
  for (j in 1:ncol(tmp)) {
    tmp_col <- get_consensus(tmp[,j])
    tmp_mot <- c(tmp_mot, tmp_col)
  }
  cons_motifs <- c(cons_motifs, paste(tmp_mot, collapse = ''))
}

cons_motifs_final <- cons_motifs %>%
  gsub('T', 'U', .) %>%
  set_names(gsub('^.*predictive_kmers/', '', gsub('\\.csv','', motifs)))

#gl pos
gl_pos <- readDNAStringSet(here('data/gkmSVM/gl/glamb_pos.fa'))
gl_pos_motifs <- motifs[c(1:5,11,13)]

mot_df_gl_pos <- data.frame(id = NA, pos = NA, motif = NA, cons_seq = NA) %>%
  na.omit()

for (i in seq_along(gl_pos_motifs)) {
  motif <- fread(gl_pos_motifs[i]) %>%
    as.matrix() %>%
    set_rownames(c('A','C','G','T'))
  mot_sc <- scan_sequences(motifs = motif, sequences = gl_pos) %>%
    as.data.frame()
  tmp_cts <- data.frame(id = NA, pos = NA) %>%
    na.omit()
  
  for (j in 1:nrow(mot_sc)) {
    tmp <- data.frame(id = rep(mot_sc$sequence[j], (mot_sc$stop[j]-mot_sc$start[j]+1)),
                      pos = round(mean(c(mot_sc$start[j], mot_sc$stop[j]))))
    tmp_cts <- rbind(tmp_cts, tmp)
  }
  
  tmp_cts$motif <- gsub('^.*predictive_kmers/', '', gsub('\\.csv','', gl_pos_motifs[i]))
  tmp_cts$cons_seq <- cons_motifs_final[match(tmp_cts$motif, 
                                              names(cons_motifs_final))]
  mot_df_gl_pos <- rbind(mot_df_gl_pos, tmp_cts)
}

pdf(here("figs/S3D_1_raw.pdf"), width = 7, height = 3)
ggdensity(mot_df_gl_pos %>%
            mutate(pos = pos-40) %>%
            filter(cons_seq %in% c('AGAANAUSWSUK','MWWGUAAAYWNN',
                                   'NURWGURAAYNCUWRRKA','YUWNNUGUAAAYN')), 
          x = 'pos', fill = 'cons_seq') +
  xlim(-40,40) +
  rremove('xlab') +
  geom_vline(xintercept = 0, linetype="dashed") +
  theme(legend.position="right") +
  rremove("legend.title")
dev.off()

pdf(here("figs/S3D_2_raw.pdf"))
ggseqlogo(fread(motifs[3]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3D_3_raw.pdf"))
ggseqlogo(fread(motifs[11]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3D_4_raw.pdf"))
ggseqlogo(fread(motifs[13]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3D_5_raw.pdf"))
ggseqlogo(fread(motifs[1]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

#gl neg
gl_neg <- readDNAStringSet(here('data/gkmSVM/gl/glamb_neg.fa'))
gl_neg_motifs <- motifs[c(9,18,19)]

mot_df_gl_neg <- data.frame(id = NA, pos = NA, motif = NA, cons_seq = NA) %>%
  na.omit()

for (i in seq_along(gl_neg_motifs)) {
  motif <- fread(gl_neg_motifs[i]) %>%
    as.matrix() %>%
    set_rownames(c('A','C','G','T'))
  mot_sc <- scan_sequences(motifs = motif, sequences = gl_neg) %>%
    as.data.frame()
  tmp_cts <- data.frame(id = NA, pos = NA) %>%
    na.omit()
  
  for (j in 1:nrow(mot_sc)) {
    tmp <- data.frame(id = rep(mot_sc$sequence[j], (mot_sc$stop[j]-mot_sc$start[j]+1)),
                      pos = round(mean(c(mot_sc$start[j], mot_sc$stop[j]))))
    tmp_cts <- rbind(tmp_cts, tmp)
  }
  
  tmp_cts$motif <- gsub('^.*predictive_kmers/', '', gsub('\\.csv','', gl_neg_motifs[i]))
  tmp_cts$cons_seq <- cons_motifs_final[match(tmp_cts$motif, 
                                              names(cons_motifs_final))]
  mot_df_gl_neg <- rbind(mot_df_gl_neg, tmp_cts)
}

pdf(here("figs/S3D_6_raw.pdf"), width = 7, height = 3)
ggdensity(mot_df_gl_neg %>%
            mutate(pos = pos-40) %>%
            filter(cons_seq %in% c('NWCCAGUGAAGAAN','GWARARNVDWKK')), 
          x = 'pos', fill = 'cons_seq') +
  xlim(-40,40) +
  rremove('xlab') +
  geom_vline(xintercept = 0, linetype="dashed") +
  theme(legend.position="right") +
  rremove("legend.title")
dev.off()

pdf(here("figs/S3D_7_raw.pdf"))
ggseqlogo(fread(motifs[19]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3D_8_raw.pdf"))
ggseqlogo(fread(motifs[9]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()


###S3E###
#gs pos
gs_pos <- readDNAStringSet(here('data/gkmSVM/gs/gsm_pos.fa'))
gs_pos_motifs <- motifs[c(1,11:17)]

mot_df_gs_pos <- data.frame(id = NA, pos = NA, motif = NA, cons_seq = NA) %>%
  na.omit()

for (i in seq_along(gs_pos_motifs)) {
  motif <- fread(gs_pos_motifs[i]) %>%
    as.matrix() %>%
    set_rownames(c('A','C','G','T'))
  mot_sc <- scan_sequences(motifs = motif, sequences = gs_pos) %>%
    as.data.frame()
  tmp_cts <- data.frame(id = NA, pos = NA) %>%
    na.omit()
  
  for (j in 1:nrow(mot_sc)) {
    tmp <- data.frame(id = rep(mot_sc$sequence[j], (mot_sc$stop[j]-mot_sc$start[j]+1)),
                      pos = round(mean(c(mot_sc$start[j], mot_sc$stop[j]))))
    tmp_cts <- rbind(tmp_cts, tmp)
  }
  
  tmp_cts$motif <- gsub('^.*predictive_kmers/', '', gsub('\\.csv','', gs_pos_motifs[i]))
  tmp_cts$cons_seq <- cons_motifs_final[match(tmp_cts$motif, 
                                              names(cons_motifs_final))]
  mot_df_gs_pos <- rbind(mot_df_gs_pos, tmp_cts)
}

pdf(here("figs/S3E_1_raw.pdf"), width = 7, height = 3)
ggdensity(mot_df_gs_pos %>%
            mutate(pos = pos-40) %>%
            filter(cons_seq %in% c('NAAGUGAAYWDA','MWWGUAAAYWNN',
                                   'NURWGURAAYNCUWRRKA','YUWNNUGUAAAYN',
                                   'YNBYWGUAAAYN','SRMNAKNWSUKADM')), 
          x = 'pos', fill = 'cons_seq') +
  xlim(-40,40) +
  rremove('xlab') +
  geom_vline(xintercept = 0, linetype="dashed") +
  theme(legend.position="right") +
  rremove("legend.title")
dev.off()

pdf(here("figs/S3E_2_raw.pdf"))
ggseqlogo(fread(motifs[1]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3E_3_raw.pdf"))
ggseqlogo(fread(motifs[11]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3E_4_raw.pdf"))
ggseqlogo(fread(motifs[12]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3E_5_raw.pdf"))
ggseqlogo(fread(motifs[13]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3E_6_raw.pdf"))
ggseqlogo(fread(motifs[14]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3E_7_raw.pdf"))
ggseqlogo(fread(motifs[15]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

#gs neg
gs_neg <- readDNAStringSet(here('data/gkmSVM/gs/gsm_neg.fa'))
gs_neg_motifs <- motifs[c(8,9)]

mot_df_gs_neg <- data.frame(id = NA, pos = NA, motif = NA, cons_seq = NA) %>%
  na.omit()

for (i in seq_along(gs_neg_motifs)) {
  motif <- fread(gs_neg_motifs[i]) %>%
    as.matrix() %>%
    set_rownames(c('A','C','G','T'))
  mot_sc <- scan_sequences(motifs = motif, sequences = gs_neg) %>%
    as.data.frame()
  tmp_cts <- data.frame(id = NA, pos = NA) %>%
    na.omit()
  
  for (j in 1:nrow(mot_sc)) {
    tmp <- data.frame(id = rep(mot_sc$sequence[j], (mot_sc$stop[j]-mot_sc$start[j]+1)),
                      pos = round(mean(c(mot_sc$start[j], mot_sc$stop[j]))))
    tmp_cts <- rbind(tmp_cts, tmp)
  }
  
  tmp_cts$motif <- gsub('^.*predictive_kmers/', '', gsub('\\.csv','', gs_neg_motifs[i]))
  tmp_cts$cons_seq <- cons_motifs_final[match(tmp_cts$motif, 
                                              names(cons_motifs_final))]
  mot_df_gs_neg <- rbind(mot_df_gs_neg, tmp_cts)
}

pdf(here("figs/S3E_8_raw.pdf"), width = 7, height = 3)
ggdensity(mot_df_gs_neg %>%
            mutate(pos = pos-40) %>%
            filter(cons_seq %in% c('NWCCAGUGAAGAAN','RAGUGAAGSGMWG')), 
          x = 'pos', fill = 'cons_seq') +
  xlim(-40,40) +
  rremove('xlab') +
  geom_vline(xintercept = 0, linetype="dashed") +
  theme(legend.position="right") +
  rremove("legend.title")
dev.off()

pdf(here("figs/S3E_9_raw.pdf"))
ggseqlogo(fread(motifs[8]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

pdf(here("figs/S3E_10_raw.pdf"))
ggseqlogo(fread(motifs[9]) %>%
            as.matrix() %>%
            set_rownames(c('A','C','G','U')))
dev.off()

###S3G###
msa1 <- here('./data/prem_cl/gl_dup.fa')
pdf(here('./figs/S3G_raw.pdf'))
ggmsa(msa1, none_bg = TRUE, font = "helvetical", char_width = .8, seq_name = TRUE) 
dev.off()


###S3H###
msa2 <- here('./data/prem_cl/gs_dup.fa')
pdf(here('./figs/S3H_raw.pdf'))
ggmsa(msa2, none_bg = TRUE, font = "helvetical", char_width = .8, seq_name = TRUE) 
dev.off()