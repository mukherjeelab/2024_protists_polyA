library(tidyverse)
library(dplyr)
library(magrittr)
library(data.table)
library(Biostrings)
library(universalmotif)
library(ggpubr)
library(ggseqlogo)
library(here)

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