library(tidyverse)
library(magrittr)
library(data.table)
library(matrixStats)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(pheatmap)
library(viridis)
library(ggpubr)
library(here)

protist_theme <- theme_pubr() +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
theme_set(protist_theme)

###S1A###
breakslist = seq(0, 1, by = .1)

ent_expr_cor <- read_tsv(here('./data/featureCounts/ent_feature_counts_cleaned.tsv.gz')) %>%
  filter(ent.1 > 0 & ent.2 > 0 & ent.3 > 0) %>%
  dplyr::select(c(1:3)) %>%
  set_colnames(c('rep1','rep2','rep3')) %>%
  cor(., method = c('spearman'))
ent_expr_htmp <- pheatmap(ent_expr_cor, 
                            clustering_distance_rows = "euclidean",
                            clustering_distance_cols = "euclidean", 
                            clustering_method = "complete", 
                            color = viridis(10, option = "E", direction = -1),
                            breaks = breakslist,
                            border_color = "black", 
                            main = 'E. histolitica')
pdf(here('./figs/S1A_raw_eh.pdf'), width = 3, height = 3.2)
ent_expr_htmp
dev.off()
neg_expr_cor <- read_tsv(here('./data/featureCounts/neg_feature_counts_cleaned.tsv.gz')) %>%
  filter(neg.1 > 0 & neg.2 > 0 & neg.3 > 0) %>%
  dplyr::select(c(1:3)) %>%
  set_colnames(c('rep1','rep2','rep3')) %>%
  cor(., method = c('spearman'))
neg_expr_htmp <- pheatmap(neg_expr_cor, 
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean", 
                          clustering_method = "complete", 
                          color = viridis(10, option = "E", direction = -1),
                          breaks = breakslist,
                          border_color = "black", 
                          main = 'N. gruberi')
pdf(here('./figs/S1A_raw_ng.pdf'), width = 3, height = 3.2)
neg_expr_htmp
dev.off()

trich_expr_cor <- read_tsv(here('./data/featureCounts/trich_feature_counts_cleaned.tsv.gz')) %>%
  filter(trich.1 > 0 & trich.2 > 0 & trich.3 > 0) %>%
  dplyr::select(c(1:3)) %>%
  set_colnames(c('rep1','rep2','rep3')) %>%
  cor(., method = c('spearman'))
trich_expr_htmp <- pheatmap(trich_expr_cor, 
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean", 
                          clustering_method = "complete", 
                          color = viridis(10, option = "E", direction = -1),
                          breaks = breakslist,
                          border_color = "black", 
                          main = 'T. vaginalis')
pdf(here('./figs/S1A_raw_tv.pdf'), width = 3, height = 3.2)
trich_expr_htmp
dev.off()

trif_expr_cor <- read_tsv(here('./data/featureCounts/trif_feature_counts_cleaned.tsv.gz')) %>%
  filter(trif.1 > 0 & trif.2 > 0 & trif.3 > 0) %>%
  dplyr::select(c(1:3)) %>%
  set_colnames(c('rep1','rep2','rep3')) %>%
  cor(., method = c('spearman'))
trif_expr_htmp <- pheatmap(trif_expr_cor, 
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean", 
                          clustering_method = "complete", 
                          color = viridis(10, option = "E", direction = -1),
                          breaks = breakslist,
                          border_color = "black", 
                          main = 'T. foetus')
pdf(here('./figs/S1A_raw_tf.pdf'), width = 3, height = 3.2)
trif_expr_htmp
dev.off()

glamb_expr_cor <- read_tsv(here('./data/featureCounts/glamb_feature_counts_cleaned.tsv.gz')) %>%
  filter(glamb.1 > 0 & glamb.2 > 0 & glamb.3 > 0) %>%
  dplyr::select(c(1:3)) %>%
  set_colnames(c('rep1','rep2','rep3')) %>%
  cor(., method = c('spearman'))
glamb_expr_htmp <- pheatmap(glamb_expr_cor, 
                            clustering_distance_rows = "euclidean",
                            clustering_distance_cols = "euclidean", 
                            clustering_method = "complete", 
                            color = viridis(10, option = "E", direction = -1),
                            breaks = breakslist,
                            border_color = "black", 
                            main = 'G. lamblia A')
pdf(here('./figs/S1A_raw_gl.pdf'), width = 3, height = 3.2)
glamb_expr_htmp
dev.off()

gsm_expr_cor <- read_tsv(here('./data/featureCounts/gsm_feature_counts_cleaned.tsv.gz')) %>%
  filter(gsm.1 > 0 & gsm.2 > 0 & gsm.3 > 0) %>%
  dplyr::select(c(1:3)) %>%
  set_colnames(c('rep1','rep2','rep3')) %>%
  cor(., method = c('spearman'))
gsm_expr_htmp <- pheatmap(gsm_expr_cor, 
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean", 
                          clustering_method = "complete", 
                          color = viridis(10, option = "E", direction = -1),
                          breaks = breakslist,
                          border_color = "black", 
                          main = 'G. lamblia B')
pdf(here('./figs/S1A_raw_gs.pdf'), width = 3, height = 3.2)
gsm_expr_htmp
dev.off()

gmur_expr_cor <- read_tsv(here('./data/featureCounts/gmur_feature_counts_cleaned.tsv.gz')) %>%
  filter(gmur.1 > 0 & gmur.2 > 0) %>%
  dplyr::select(c(1:2)) %>%
  set_colnames(c('rep1','rep2')) %>%
  cor(., method = c('spearman'))
gmur_expr_htmp <- pheatmap(gmur_expr_cor, 
                          clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean", 
                          clustering_method = "complete", 
                          color = viridis(10, option = "E", direction = -1),
                          breaks = breakslist,
                          border_color = "black", 
                          main = 'G. muris')
pdf(here('./figs/S1A_raw_gm.pdf'), width = 3, height = 3.2)
gmur_expr_htmp
dev.off()


###S1B###
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


###S1C###
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


###S1D###
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


###S1E###
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
