library(tidyverse)
library(magrittr)
library(broom)
library(Biostrings)
library(ggpubr)
library(here)

protist_theme <- theme_pubr() +
  theme(text = element_text(size = 14))+
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14))
theme_set(protist_theme)

gl <- read_tsv(here('./data/aux_elements/ugua_df_glamb.tsv.gz')) %>%
  mutate(dataset = 'G. lamblia')
gs <- read_tsv(here('./data/aux_elements/ugua_df_gsm.tsv.gz')) %>%
  mutate(dataset = 'G. lamblia B-GS')
gm <- read_tsv(here('./data/aux_elements/ugua_df_gmur.tsv.gz')) %>%
  mutate(dataset = 'G. muris')
#eh <- read_tsv(here('./data/ugua_df_ent.tsv.gz')) %>%
#  mutate(dataset = 'E. histolitica')
ng <- read_tsv(here('./data/aux_elements/ugua_df_neg.tsv.gz')) %>%
  mutate(dataset = 'N. gruberi')
tf <- read_tsv(here('./data/aux_elements/ugua_df_trif.tsv.gz')) %>%
  mutate(dataset = 'T. foetus')
tv <- read_tsv(here('./data/aux_elements/ugua_df_trich.tsv.gz')) %>%
  mutate(dataset = 'T. vaginalis')
hs <- read_tsv(here('./data/aux_elements/ugua_df_hs.tsv.gz')) %>%
  set_colnames(colnames(gl))

ugua_df <- rbind(hs, gl, gs, gm, ng, tf, tv) %>%
  mutate(dataset = factor(dataset, levels = c('H. sapiens', 'N. gruberi', 'T. foetus', 'T. vaginalis', 'G. lamblia', 'G. lamblia B-GS', 'G. muris')),
         GU_signals = factor(GU_signals, levels = c("UUAG", "UAGU", "GUUA", "GAUU", "AGUU", "GUAU", "UGAU", "AUUG", "UAUG", "UUGA", "AUGU", "UGUA")),
         signal = factor(if_else(GU_signals == 'UGUA', 'UGUA', 'shuffled'), levels = c('UGUA', 'shuffled')))

ugua_df_summary <- ugua_df %>%
  group_by(dataset, signal) %>%
  summarise(
    sd = sd(occurence),
    occurence = round(mean(occurence))
  ) %>%
  replace(is.na(.), 0) 

ugua_df_summary2 <- ugua_df %>%
  group_by(dataset, signal) %>%
  summarise(
    occurence = sum(occurence)
  ) %>%
  replace(is.na(.), 0) %>%
  group_by(dataset) %>%
  mutate(freq = occurence/sum(occurence),
         label_y = c(0.95, 0.5),
         signal = signal)

pdf(here('./figs/ugua_shuf_all.pdf'), width = 10.5, height = 3)  
ggplot(./data/ugua_df, aes(signal, occurence, color = color)) +
  geom_bar(stat = "identity", data = ./data/ugua_df_summary,
            fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2)) +
  scale_color_manual(values = c('black', '#0054ffff')) +
  facet_wrap(vars(dataset), nrow = 1, , scales="free_y") +
  rremove('xlab') +
  rremove('legend') +
  ylab('counts') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

 
stat_summary_ugua <- ./data/ugua_df_summary2 %>%
  pivot_wider(., id_cols = 'signal', names_from = 'dataset', values_from = 'occurence')
p_val_df_ugua <- data.frame(matrix(nrow = 0, ncol = 2)) %>%
  set_colnames(c('dataset', 'p.value'))
  
for (i in 3:ncol(stat_summary_ugua)) {
  ct <- stat_summary_ugua %>%
    dplyr::select(c(1,2,i)) %>%
    column_to_rownames('signal')
  tmp <- data.frame(dataset = colnames(stat_summary_ugua)[i],
                    p.value = chisq.test(ct) %>%
                      tidy %>%
                      dplyr::select('p.value'))
  p_val_df_ugua <- rbind(p_val_df_ugua, tmp)
}
#write_tsv(p_val_df_ugua, here('paper01/ugua_vs_shuf_pvals.tsv.gz'))
stat_summary_ugua_4plot <- p_val_df_ugua %>%
  mutate(group1 = 'H. sapiens',
         group2 = dataset,
         pval = c('0.047', '3.65e-100', '4.48e-98', '2.11e-19', '5.27e-25', '4.29e-22'))
theme_set(protist_theme)

pdf(here('./figs/mosaic_shuf_ugua.pdf'), height = 5) 
ggbarplot(./data/ugua_df_summary2, x = 'dataset', y = 'freq', fill = 'signal',
          palette = c('blue', 'gray')) +
  geom_text(
    aes(y = label_y, label = occurence)
  ) +
  stat_pvalue_manual(stat_summary_ugua_4plot, label = 'pval', 
                     y.position = c(1.1,1.25,1.4,1.55,1.7,1.85)) +
  ylim(0,1.9) +
  rremove('xlab') +
  rremove('ylab') +
  rremove('legend.title') +
  theme(legend.position = 'right') +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.line = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1))
dev.off()

gl2 <- read_tsv(here('./data/aux_elements/gu-rich_df_glamb.tsv.gz')) %>%
  group_by(GU_motifs) %>%
  summarise(ratio = occurence[1]/occurence[2]) %>%
  mutate(dataset = 'G. lamblia')
gs2 <- read_tsv(here('./data/aux_elements/gu-rich_df_gsm.tsv.gz')) %>%
  group_by(GU_motifs) %>%
  summarise(ratio = occurence[1]/occurence[2]) %>%
  mutate(dataset = './data/aux_elements/G. lamblia B-GS')
gm2 <- read_tsv(here('./data/aux_elements/gu-rich_df_gmur.tsv.gz')) %>%
  group_by(GU_motifs) %>%
  summarise(ratio = occurence[1]/occurence[2]) %>%
  mutate(dataset = 'G. muris')
ng2 <- read_tsv(here('./data/aux_elements/gu-rich_df_neg.tsv.gz')) %>%
  group_by(GU_motifs) %>%
  summarise(ratio = occurence[1]/occurence[2]) %>%
  mutate(dataset = 'N. gruberi')
tf2 <- read_tsv(here('./data/aux_elements/gu-rich_df_trif.tsv.gz')) %>%
  group_by(GU_motifs) %>%
  summarise(ratio = occurence[1]/occurence[2]) %>%
  mutate(dataset = 'T. foetus')
tv2 <- read_tsv(here('./data/aux_elements/gu-rich_df_trich.tsv.gz')) %>%
  group_by(GU_motifs) %>%
  summarise(ratio = occurence[1]/occurence[2]) %>%
  mutate(dataset = 'T. vaginalis')
to_cnt2 <- gsub('U', 'T', gl2$GU_motifs)
up_seq <- readDNAStringSet(here('./paper01/RefSeq_human/upstream30.fa'))
down_seq <- readDNAStringSet(here('./paper01/RefSeq_human/downstream30.fa'))
hs2 <- data.frame(GU_motifs = gl2$GU_motifs, ratio = NA, dataset = 'H. sapiens')
for (i in seq_along(to_cnt2)) {
  hs2$ratio[i] = sum(vcountPattern(to_cnt2[i], down_seq))/sum(vcountPattern(to_cnt2[i], up_seq))
}

ug_rich <- rbind(gl2, gs2, gm2, ng2, tf2, tv2, hs2) %>%
  mutate(dataset = factor(dataset, 
                          levels = c('H. sapiens', 'N. gruberi', 'T. foetus', 
                                     'T. vaginalis', 'G. lamblia', 'G. lamblia B-GS', 
                                     'G. muris')))
'''pdf(here('./figs/gu_rich_all.pdf'), width = 10.5, height = 4.5)
ggviolin(ug_rich, x = 'dataset', y = 'ratio', ylab = 'log2(downstream/upstream) ratio',
         add = 'jitter', main = 'GU- and U-rich elements') +
  rremove('xlab') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

ug_rich_summary <- ug_rich %>%
  mutate(ratio = ifelse(is.finite(ratio) == TRUE, ratio, 0)) %>%
  group_by(dataset) %>%
  summarise(
    sd = sd(ratio),
    ratio = mean(ratio)
  ) %>%
  replace(is.na(.), 0)
pdf(here('./paper01/figs/new/gu_rich_barplot_all.pdf'), width = 10.5, height = 3)  
ggplot(ug_rich, aes(dataset, ratio)) +
  geom_bar(stat = "identity", data = ug_rich_summary,
           fill = NA, color = "black") +
  geom_jitter(position = position_jitter(0.2)) +
  scale_color_manual(values = c('black')) +
  theme_pubr() +
  rremove('xlab') +
  rremove('legend') +
  ylab('log2(downstream/upstream) ratio') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text = element_text(size = 8))
dev.off()'''


ug_rich_stat <- ug_rich %>% 
  pivot_wider(id_cols = 'GU_motifs', names_from = 'dataset', values_from = 'ratio') %>%
  dplyr::select(c(1,5:7,2:4,8))
ug_rich_stat[sapply(ug_rich_stat, is.infinite)] <- 0

p_val_df_ug_rich <- data.frame(matrix(nrow = 0, ncol = 2)) %>%
  set_colnames(c('dataset', 'p.value'))
for (i in 2:(ncol(ug_rich_stat)-1)) {
  x <- ug_rich_stat[, (ncol(ug_rich_stat))] %>%
    deframe()
  y <- ug_rich_stat[, i] %>%
    deframe()
  tmp <- data.frame(dataset = colnames(ug_rich_stat)[i],
                    p.value = wilcox.test(x = x, y = y, 
                                          paired = FALSE, alternative = 'g') %>%
                      tidy() %>%
                      dplyr::select(c('p.value'))
                    )
  p_val_df_ug_rich <- rbind(p_val_df_ug_rich, tmp)                  
}

p_val_df_ug_rich_4plot <- p_val_df_ug_rich %>%
  mutate(group1 = 'H. sapiens',
         group2 = dataset,
         pval = c('3.59e-5','2.54e-3','1.83e-6', '2.88e-6', '1.11e-4', '1.14e-4'))
pdf(here('./figs/gu_rich_boxplot_all.pdf'), height = 4.5)
ggplot(ug_rich, aes(x = dataset, y = ratio)) +
  geom_boxplot() +
  geom_jitter() +
  rremove('xlab') +
  ylab('log2(downstream/upstream) ratio') +
  ggtitle('GU- and U-rich elements') +
  theme(title = element_text(size = 14)) +
  ylim(0,13) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_pvalue_manual(p_val_df_ug_rich_4plot, label = 'pval', 
                     y.position = c(8,8.9,9.8,10.7,11.6,12.5))
dev.off()
