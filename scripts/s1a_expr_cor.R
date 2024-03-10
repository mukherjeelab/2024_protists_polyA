library(tidyverse)
library(magrittr)
library(pheatmap)
library(viridis)
library(here)

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
pdf(here('./figs/s1a_eh.pdf'), width = 3, height = 3.2)
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
pdf(here('./figs/s1a_ng.pdf'), width = 3, height = 3.2)
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
pdf(here('./figs/s1a_tv.pdf'), width = 3, height = 3.2)
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
pdf(here('./figs/s1a_tf.pdf'), width = 3, height = 3.2)
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
                            main = 'G. lamblia')
pdf(here('./figs/s1a_gl.pdf'), width = 3, height = 3.2)
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
                          main = 'G. lamblia B-GS')
pdf(here('./figs/s1a_gs.pdf'), width = 3, height = 3.2)
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
pdf(here('./figs/s1a_gm.pdf'), width = 3, height = 3.2)
gmur_expr_htmp
dev.off()
