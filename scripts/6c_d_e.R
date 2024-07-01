library(tidyverse)
library(ggpubr)
library(kmer)
library(ape)
library(Biostrings)
library(universalmotif)
library(ggseqlogo)
library(here)

protist_theme <- theme_pubr() +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
theme_set(protist_theme)

cts <- read_tsv(here('data/aux_elements/neg_hexamers_CDS_cts.tsv.gz'))
pdf(here('figs/6C_raw.pdf'), width = 3, height = 4)
ggbarplot(cts, x = 'hexamer', y = 'counts', fill = 'gray', label = TRUE) +
  theme(text = element_text(size = 14)) +
  rremove('ylab') +
  rremove('xlab')
dev.off()

neg_8mers <- read_tsv(here('data/gkmSVM/ng/neg_8mers_AAUAAA.out.gz'), col_names = FALSE) %>%
  arrange(desc(X2))
top_50 <- neg_8mers %>%
  head(n = 50)
seq_top50 <- DNAStringSet(top_50$X1) %>%
  set_names(top_50$X1)
seq_top50_bin <- as.DNAbin(seq_top50)
top_50_clust <- cluster(seq_top50_bin, k = 6)
#plot(top_50_clust)

pwm1 <- seq_top50 %>%
  .[names(.) %in% c('TTTGTAAA', 'TTTTGTAA', 'ATTGTAAA','AATTGTAA')] %>%
  PWM()
pwm2 <- seq_top50 %>%
  .[names(.) %in% c('TTTTTTTT', 'ATTTTTTT', 'CTTTTTTT','GTTTTTTT','TTTTTTTG')] %>%
  PWM()

neg_pos <- readDNAStringSet('data/gkmSVM/ng/neg_pos_AAUAAA.fa.gz')
neg_neg <- readDNAStringSet('data/gkmSVM/ng/neg_neg_AAUAAA.fa.gz')
sc1_p <- scan_sequences(pwm1, neg_pos) %>%
  as.data.frame() %>%
  mutate(region = "3'UTR")
sc2_p <- scan_sequences(pwm2, neg_pos) %>%
  as.data.frame() %>%
  mutate(region = "3'UTR")
sc1_n <- scan_sequences(pwm1, neg_neg) %>%
  as.data.frame() %>%
  mutate(region = "CDS")
sc1_f <- rbind(sc1_p, sc1_n) %>%
  mutate(region = factor(region))
sc2_n <- scan_sequences(pwm2, neg_neg) %>%
  as.data.frame() %>%
  mutate(region = "CDS")
sc2_f <- rbind(sc2_p, sc2_n) %>%
  mutate(region = factor(region))
pdf(here('figs/6D_1_raw.pdf'), width = 5, height = 3)
ggdensity(sc1_f %>% mutate(start = start -40), x ='start', y = '..count..', xlab = 'position',
          fill = 'region') +
  xlim(-40,40) +
  rremove('xlab') +
  geom_vline(xintercept = 0, linetype="dashed") +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
dev.off()

pdf(here('figs/6E_1_raw.pdf'), width = 5, height = 3)
ggdensity(sc2_f %>% mutate(start = start -40), x ='start', y = '..count..', xlab = 'position',
          fill = 'region') +
  xlim(-40,40) +
  rremove('xlab') +
  geom_vline(xintercept = 0, linetype="dashed") +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
dev.off()

pwm1_rna <- pwm1 %>%
  set_rownames(c('A','C','G','U'))
pwm2_rna <- pwm2 %>%
  set_rownames(c('A','C','G','U'))

pdf(here('figs/6D_2_raw.pdf'), width = 4, height = 3)
ggseqlogo(pwm1_rna)
dev.off()

pdf(here('figs/6E_2_raw.pdf'), width = 4, height = 3)
ggseqlogo(pwm2_rna)
dev.off()