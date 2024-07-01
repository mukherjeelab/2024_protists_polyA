library(dplyr)
library(tidyverse)
library(data.table)
library(magrittr)
library(matrixStats)
library(ggpubr)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)

#S1B
#gl
annot_glamb1 <- fread(here('data/featureCounts/glamb1_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_glamb1 <- fread(here('data/nanopolish_res/glamb-1_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_glamb1, by = c('readname' = 'V1')) %>%
  na.omit()
annot_glamb2 <- fread(here('data/featureCounts/glamb2_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_glamb2 <- fread(here('data/nanopolish_res/glamb-2_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_glamb2, by = c('readname' = 'V1')) %>%
  na.omit()
annot_glamb3 <- fread(here('data/featureCounts/glamb3_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_glamb3 <- fread(here('data/nanopolish_res/glamb-3_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_glamb3, by = c('readname' = 'V1')) %>%
  na.omit()
glamb_polyA_distrib <- rbind(polyA_glamb1,polyA_glamb2,polyA_glamb3) %>%
  dplyr::rename('GeneID' = 'V4') %>%
  group_by(GeneID) %>%
  mutate(gene_counts = dplyr::n()) %>%
  ungroup() %>%
  mutate(`3'UTR annotation` = factor(if_else(polya_length < 30, 'discarded','kept'))) %>%
  filter(gene_counts >= 10, polya_length <= 200) %>%
  mutate(dataset = 'G. lamblia A')
#gs
annot_gsm1 <- fread(here('data/featureCounts/gsm1_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_gsm1 <- fread(here('data/nanopolish_res/gsm-1_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_gsm1, by = c('readname' = 'V1')) %>%
  na.omit()
annot_gsm2 <- fread(here('data/featureCounts/gsm2_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_gsm2 <- fread(here('data/nanopolish_res/gsm-2_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_gsm2, by = c('readname' = 'V1')) %>%
  na.omit()
annot_gsm3 <- fread(here('data/featureCounts/gsm3_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_gsm3 <- fread(here('data/nanopolish_res/gsm-3_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_gsm3, by = c('readname' = 'V1')) %>%
  na.omit()
gsm_polyA_distrib <- rbind(polyA_gsm1,polyA_gsm2,polyA_gsm3) %>%
  dplyr::rename('GeneID' = 'V4') %>%
  group_by(GeneID) %>%
  mutate(gene_counts = dplyr::n()) %>%
  ungroup() %>%
  mutate(`3'UTR annotation` = factor(if_else(polya_length < 30, 'discarded','kept'))) %>%
  filter(gene_counts >= 10, polya_length <= 200) %>%
  mutate(dataset = 'G. lamblia B')
#gm
annot_gmur1 <- fread(here('data/featureCounts/gmur1_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_gmur1 <- fread(here('data/nanopolish_res/gmur-1_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_gmur1, by = c('readname' = 'V1')) %>%
  na.omit()
annot_gmur2 <- fread(here('data/featureCounts/gmur2_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_gmur2 <- fread(here('data/nanopolish_res/gmur-2_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_gmur2, by = c('readname' = 'V1')) %>%
  na.omit()
gmur_polyA_distrib <- rbind(polyA_gmur1,polyA_gmur2) %>%
  dplyr::rename('GeneID' = 'V4') %>%
  group_by(GeneID) %>%
  mutate(gene_counts = dplyr::n()) %>%
  ungroup() %>%
  mutate(`3'UTR annotation` = factor(if_else(polya_length < 30, 'discarded','kept'))) %>%
  filter(gene_counts >= 10, polya_length <= 200) %>%
  mutate(dataset = 'G. muris')
#eh
annot_ent1 <- fread(here('data/featureCounts/ent1_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_ent1 <- fread(here('data/nanopolish_res/ent-1_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_ent1, by = c('readname' = 'V1')) %>%
  na.omit()
annot_ent2 <- fread(here('data/featureCounts/ent2_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_ent2 <- fread(here('data/nanopolish_res/ent-2_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_ent2, by = c('readname' = 'V1')) %>%
  na.omit()
annot_ent3 <- fread(here('data/featureCounts/ent3_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_ent3 <- fread(here('data/nanopolish_res/ent-3_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_ent3, by = c('readname' = 'V1')) %>%
  na.omit()
ent_polyA_distrib <- rbind(polyA_ent1,polyA_ent2,polyA_ent3) %>%
  dplyr::rename('GeneID' = 'V4') %>%
  group_by(GeneID) %>%
  mutate(gene_counts = dplyr::n()) %>%
  ungroup() %>%
  mutate(`3'UTR annotation` = factor(if_else(polya_length < 30, 'discarded','kept'))) %>%
  filter(gene_counts >= 10, polya_length <= 200) %>%
  mutate(dataset = 'E. histolytica')
#neg
annot_neg1 <- fread(here('data/featureCounts/neg1_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_neg1 <- fread(here('data/nanopolish_res/neg-1_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_neg1, by = c('readname' = 'V1')) %>%
  na.omit()
annot_neg2 <- fread(here('data/featureCounts/neg2_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_neg2 <- fread(here('data/nanopolish_res/neg-2_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_neg2, by = c('readname' = 'V1')) %>%
  na.omit()
annot_neg3 <- fread(here('data/featureCounts/neg3_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_neg3 <- fread(here('data/nanopolish_res/neg-3_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_neg3, by = c('readname' = 'V1')) %>%
  na.omit()
neg_polyA_distrib <- rbind(polyA_neg1,polyA_neg2,polyA_neg3) %>%
  dplyr::rename('GeneID' = 'V4') %>%
  group_by(GeneID) %>%
  mutate(gene_counts = dplyr::n()) %>%
  ungroup() %>%
  mutate(`3'UTR annotation` = factor(if_else(polya_length < 15, 'discarded','kept'))) %>%
  filter(gene_counts >= 10, polya_length <= 200) %>%
  mutate(dataset = 'N. gruberi')
#tv
annot_trich1 <- fread(here('data/featureCounts/trich1_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_trich1 <- fread(here('data/nanopolish_res/trich-1_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_trich1, by = c('readname' = 'V1')) %>%
  na.omit()
annot_trich2 <- fread(here('data/featureCounts/trich2_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_trich2 <- fread(here('data/nanopolish_res/trich-2_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_trich2, by = c('readname' = 'V1')) %>%
  na.omit()
annot_trich3 <- fread(here('data/featureCounts/trich3_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_trich3 <- fread(here('data/nanopolish_res/trich-3_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_trich3, by = c('readname' = 'V1')) %>%
  na.omit()
trich_polyA_distrib <- rbind(polyA_trich1,polyA_trich2,polyA_trich3) %>%
  dplyr::rename('GeneID' = 'V4') %>%
  group_by(GeneID) %>%
  mutate(gene_counts = dplyr::n()) %>%
  ungroup() %>%
  mutate(`3'UTR annotation` = factor(if_else(polya_length < 30, 'discarded','kept'))) %>%
  filter(gene_counts >= 10, polya_length <= 200) %>%
  mutate(dataset = 'T. vaginalis')
#tf
annot_trif1 <- fread(here('data/featureCounts/trif1_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_trif1 <- fread(here('data/nanopolish_res/trif-1_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_trif1, by = c('readname' = 'V1')) %>%
  na.omit()
annot_trif2 <- fread(here('data/featureCounts/trif2_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_trif2 <- fread(here('data/nanopolish_res/trif-2_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_trif2, by = c('readname' = 'V1')) %>%
  na.omit()
annot_trif3 <- fread(here('data/featureCounts/trif3_feature_counts_annot_filtered.tsv.gz')) %>%
  dplyr::select(c(1,4))
polyA_trif3 <- fread(here('data/nanopolish_res/trif-3_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS') %>%
  distinct(readname, .keep_all = T) %>%
  dplyr::select(c(1,9)) %>%
  left_join(., annot_trif3, by = c('readname' = 'V1')) %>%
  na.omit()
trif_polyA_distrib <- rbind(polyA_trif1,polyA_trif2,polyA_trif3) %>%
  dplyr::rename('GeneID' = 'V4') %>%
  group_by(GeneID) %>%
  mutate(gene_counts = dplyr::n()) %>%
  ungroup() %>%
  mutate(`3'UTR annotation` = factor(if_else(polya_length < 30, 'discarded','kept'))) %>%
  filter(gene_counts >= 10, polya_length <= 200) %>%
  mutate(dataset = 'T. foetus')

distr_df <- rbind(glamb_polyA_distrib,gsm_polyA_distrib,gmur_polyA_distrib,
                  ent_polyA_distrib, neg_polyA_distrib, trich_polyA_distrib,
                  trif_polyA_distrib) %>%
  mutate(dataset = factor(dataset, levels = c('E. histolytica','N. gruberi','T. foetus', 'T. vaginalis', 'G. lamblia A','G. lamblia B','G. muris')))

pdf(here('figs/S1B_raw.pdf'))
gghistogram(distr_df, x = 'polya_length', add = 'median', bins = 20,  fill = "3'UTR annotation", palette = c('red', 'lightgray'), xlab = 'polyA length', facet.by = 'dataset', ncol = 4)+
  rremove('ylab') + 
  theme(axis.title.x=element_text(size=14))+
  theme(strip.text.x=element_text(size=14))+
  theme(axis.text.x=element_text(size=14))+
  theme(axis.text.y=element_text(size=14))+
  theme(strip.text.x = element_text(size = 14))
dev.off()

#S1E
file_list <- list.files(here('data/sequences/'), pattern = '*_3UTR_combined.fa')
creatures <- c('E. histolitica','G. lamblia A','G. muris','G. lamblia B','N. gruberi','T. vaginalis','T. foetus')
len_df <- data_frame(utr3_length = NA, dataset = NA)


for(i in 1:length(file_list)) {
  seqs <- readDNAStringSet(here('data/sequences/',file_list[i]))
  len <- unlist(lapply(seqs, length), use.names = T) %>%
    as.data.frame() %>%
    mutate(dataset = creatures[i]) %>%
    dplyr::rename('utr3_length'='.')
  len_df <- bind_rows(len_df, len)
}

#adding human data
if (file.exists("data/annotation/txdb.gencode26.sqlite.gz")) {
  system("gunzip data/annotation/txdb.gencode26.sqlite.gz")
}

txdb <- loadDb(here('data/annotation/txdb.gencode26.sqlite'))
utr_3 <- threeUTRsByTranscript(txdb, use.names=T)
utr_seqs <- getSeq(Hsapiens, unlist(utr_3))
len_hs <- unlist(lapply(utr_seqs, length), use.names = T) %>%
  as.data.frame() %>%
  mutate(dataset = 'H. sapiens') %>%
  dplyr::rename('utr3_length'='.')

len_df_final <- bind_rows(len_df, len_hs) %>%
  na.omit() %>%
  mutate(dataset = factor(dataset, levels = c('E. histolitica','N. gruberi','T. foetus', 'T. vaginalis', 'G. lamblia A','G. lamblia B','G. muris','H. sapiens')),
         log10_length = log10(utr3_length))

pdf(here('figs/S1E_raw.pdf'), width = 4, height = 4)
ggviolin(len_df_final, x = 'dataset', y = 'log10_length',
         ylab = 'log10(length)', add = 'boxplot') +
  rremove('xlab') +
  theme(axis.title.y=element_text(size=14))+
  theme(axis.text.x=element_text(size=14, angle = 45, hjust = 1))+
  theme(axis.text.y=element_text(size=14))
dev.off()
