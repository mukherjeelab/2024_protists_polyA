library(tidyverse)
library(data.table)
library(magrittr)
library(matrixStats))
library(here)

gsm_ref <- fread(here('data/annotation/gsm_reference_annot.tsv.gz'))

gsm_ref_3utr <- gsm_ref %>%
  dplyr::select(c('seq_id','primary_tag','start','end','gene_id')) %>%
  filter(primary_tag == 'three_prime_utr')%>%
  dplyr::rename(c('utr3_annot_start'='start','utr3_annot_end'='end','GeneID'='gene_id','Chr'='seq_id')) %>%
  dplyr::select(-c('primary_tag'))

gsm_ref_5utr <- gsm_ref %>%
  dplyr::select(c('seq_id','primary_tag','start','end','gene_id')) %>%
  filter(primary_tag == 'five_prime_utr') %>%
  dplyr::rename(c('utr5_annot_start'='start','utr5_annot_end'='end','GeneID'='gene_id','Chr'='seq_id'))%>%
  dplyr::select(-c('primary_tag'))

fc_gsm <- fread(here('data/featureCounts/gsm_feature_counts_cleaned.tsv.gz')) %>%
  left_join(.,gsm_ref_3utr) %>%
  left_join(.,gsm_ref_5utr) %>%
  mutate(Start = case_when(is.na(utr3_annot_start) == F & Strand == '-' ~ as.integer(utr3_annot_end + 1),
                           is.na(utr5_annot_start) == F & Strand == '+' ~ as.integer(utr5_annot_end + 1),
                           TRUE ~ Start)) %>%
  mutate(End = case_when(is.na(utr3_annot_start) == F & Strand == '+' ~ as.integer(utr3_annot_start - 1),
                         is.na(utr5_annot_start) == F & Strand == '-' ~ as.integer(utr5_annot_start - 1),
                         TRUE ~ End)) %>%
  dplyr::select(-c('utr3_annot_start','utr3_annot_end','utr5_annot_start','utr5_annot_end'))

#rep1
fc_gsm1 <- fc_gsm %>%
  filter(.[,1] > 0) %>%
  dplyr::select(c(1,'GeneID':'Length'))
polyA_gsm1 <- fread(here('data/nanopolish_res/gsm-1_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS' & polya_length >= 30)
annot_gsm1 <- fread(here('data/featureCounts/gsm1_feature_counts_annot_filtered.tsv.gz')) %>%
  filter(V2 == 'Assigned' & V1 %in% polyA_gsm1$readname) %>%
  distinct(.,V1, .keep_all = T) %>%
  dplyr::select(c(1,4)) %>%
  set_colnames(c('readName','GeneID')) %>%
  left_join(.,fc_gsm1)
compare_gsm1 <- fread(here('data/alignment_coords/gsm-1.sorted.bed.gz')) %>%
  distinct(V4, .keep_all = T) %>%
  dplyr::select(-c(5)) %>%
  set_colnames(c('Chr','readStart','readEnd','readName','Strand')) %>% 
  left_join(annot_gsm1,.) %>%
  dplyr::select(c(1,2,4,5,9,6,10,7,8)) %>%
  na.omit() %>%
  filter(if_else(Strand == '+', readEnd > End & readStart < End, readStart < Start & readEnd > Start)) %>%
  mutate(dataset = 'gsm-1')

#rep2
fc_gsm2 <- fc_gsm %>%
  filter(.[,2] > 0) %>%
  dplyr::select(c(1,'GeneID':'Length'))
polyA_gsm2 <- fread(here('data/nanopolish_res/gsm-2_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS' & polya_length >= 30)
annot_gsm2 <- fread(here('data/featureCounts/gsm2_feature_counts_annot_filtered.tsv.gz')) %>%
  filter(V2 == 'Assigned' & V1 %in% polyA_gsm2$readname) %>%
  distinct(.,V1, .keep_all = T) %>%
  dplyr::select(c(1,4)) %>%
  set_colnames(c('readName','GeneID')) %>%
  left_join(.,fc_gsm2)
compare_gsm2 <- fread(here('data/alignment_coords/gsm-2.sorted.bed.gz')) %>%
  distinct(V4, .keep_all = T) %>%
  dplyr::select(-c(5)) %>%
  set_colnames(c('Chr','readStart','readEnd','readName','Strand')) %>% 
  left_join(annot_gsm2,.) %>%
  dplyr::select(c(1,2,4,5,9,6,10,7,8)) %>%
  na.omit() %>%
  filter(if_else(Strand == '+', readEnd > End & readStart < End, readStart < Start & readEnd > Start)) %>%
  mutate(dataset = 'gsm-2')

#rep3
fc_gsm3 <- fc_gsm %>%
  filter(.[,3] > 0) %>%
  dplyr::select(c(1,'GeneID':'Length'))
polyA_gsm3 <- fread(here('data/nanopolish_res/gsm-3_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS' & polya_length >= 30)
annot_gsm3 <- fread(here('data/featureCounts/gsm3_feature_counts_annot_filtered.tsv.gz')) %>%
  filter(V2 == 'Assigned' & V1 %in% polyA_gsm3$readname) %>%
  distinct(.,V1, .keep_all = T) %>%
  dplyr::select(c(1,4)) %>%
  set_colnames(c('readName','GeneID')) %>%
  left_join(.,fc_gsm3)
compare_gsm3 <- fread(here('data/alignment_coords/gsm-3.sorted.bed.gz')) %>%
  distinct(V4, .keep_all = T) %>%
  dplyr::select(-c(5)) %>%
  set_colnames(c('Chr','readStart','readEnd','readName','Strand')) %>% 
  left_join(annot_gsm3,.) %>%
  dplyr::select(c(1,2,4,5,9,6,10,7,8)) %>%
  na.omit() %>%
  filter(if_else(Strand == '+', readEnd > End & readStart < End, readStart < Start & readEnd > Start)) %>%
  mutate(dataset = 'gsm-3')

#combining data and calculating isoforms frequencies
gsm_all <- rbind(compare_gsm1, compare_gsm2, compare_gsm3) %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  mutate(utr3_length = if_else(Strand == '+', (readEnd-End)+1, (Start-readStart)+1)) %>%
  arrange(GeneID)

gsm_groupped <- gsm_all %>%
  mutate(group_var = paste0(GeneID,';',utr3_length)) %>%
  group_by(group_var) %>%
  arrange(group_var) %>%
  summarise(counts = n()) %>%
  separate(group_var, into = c('GeneID','utr3_length'), sep = ';') %>%
  mutate(utr3_length = as.numeric(utr3_length))

gsm_iso_freq_tmp <- gsm_groupped %>% 
  filter(counts >= 10) %>%
  group_by(GeneID) %>%
  arrange(utr3_length, .by_group = T) %>%
  mutate(iso_freq = counts/sum(counts)) %>%
  filter(iso_freq >= 0.1) %>%
  mutate(diff = utr3_length - min(utr3_length)) %>%
  mutate(temp = if_else(diff <= 20, 'min', 'others')) %>%
  mutate(group_var = paste0(GeneID,'_',temp)) %>%
  group_by(group_var) %>%
  filter(temp == 'others'|temp == 'min'&counts==max(counts)) %>%
  group_by(GeneID) %>%
  mutate(diff = utr3_length - min(utr3_length)) %>%
  mutate(temp = if_else(diff <= 20, 1, diff%/%20+1)) %>%
  mutate(group_var = paste0(GeneID,'_',temp)) %>%
  group_by(group_var) %>%
  filter(counts == max(counts)) %>%
  mutate(max_len = max(utr3_length)) %>%
  distinct(max_len, .keep_all = T) %>%
  group_by(GeneID) %>%
  mutate(cl_site_n = 1:n(), sum_iso = sum(iso_freq), mean_site = mean(cl_site_n)) %>%
  mutate(category = case_when(mean_site == 1 ~ 'singles',
                              mean_site == 1.5 ~ 'doubles',
                              TRUE ~ 'multiples')) %>%
  dplyr::select(c('GeneID','utr3_length','iso_freq','sum_iso','cl_site_n', 'category'))

group_iso_multiples_1 <- gsm_iso_freq_tmp %>%
  filter(category=='multiples' & cl_site_n == 1)

group_iso_multiples_filtered <- gsm_iso_freq_tmp %>%
  filter(category=='multiples' & cl_site_n > 1) %>%
  group_by(GeneID) %>%
  mutate(diff = utr3_length - min(utr3_length), temp = if_else(diff <= 20, '_same','_diff'),
         grouping_var = paste0(GeneID,temp)) %>%
  group_by(grouping_var) %>%
  filter(temp=='_diff'|temp=='_same'&iso_freq==max(iso_freq)) %>%
  filter(utr3_length==max(utr3_length)) %>%
  ungroup() %>%
  dplyr::select(-c(7:9)) %>%
  group_by(GeneID) %>%
  mutate(cl_site_n = (1:n())+1, mean_cl_site = mean(cl_site_n), category = if_else(mean_cl_site == 2, 'doubles', category)) %>%
  dplyr::select(-c(7))

gsm_iso_freq <- rbind(subset(gsm_iso_freq_tmp, !category =='multiples'), group_iso_multiples_1, group_iso_multiples_filtered) %>%
  group_by(GeneID) %>%
  arrange(GeneID) %>%
  mutate(cl_site_n = 1:n(), sum_iso = sum(iso_freq), mean_site = mean(cl_site_n)) %>%
  mutate(category = case_when(mean_site == 1 ~ 'singles',
                              mean_site == 1.5 ~ 'doubles',
                              TRUE ~ 'multiples')) %>%
  dplyr::select(c('GeneID','utr3_length','iso_freq','sum_iso','cl_site_n', 'category'))

#calculating coordinates
gsm_3UTR_coord_full_transcr <- gsm_iso_freq %>%
  dplyr::select(c('GeneID','utr3_length')) %>%
  left_join(., subset(fc_gsm, GeneID %in% .$GeneID, select = c(GeneID:Strand))) %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(cl_site = if_else(Strand == '+', End + utr3_length, Start - utr3_length)) %>%
  mutate(Start2 = if_else(Strand == '+', Start, cl_site), End2 = if_else(Strand == '+', cl_site, End), score = 60) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'GeneID', 'score','Strand')) %>%
  filter(Start2 > 0)

gsm_3UTR_coord <- gsm_iso_freq %>%
  dplyr::select(c('GeneID','utr3_length')) %>%
  left_join(., subset(fc_gsm, GeneID %in% .$GeneID, select = c(GeneID:Strand))) %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(cl_site = if_else(Strand == '+', End + utr3_length, Start - utr3_length)) %>%
  mutate(Start2 = if_else(Strand == '+', End + 1, cl_site), End2 = if_else(Strand == '+', cl_site, Start - 1), score = 60) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'GeneID', 'score','Strand')) %>%
  filter(Start2 > 0)

gsm_5UTR_coord_full_transcr <- gsm_all %>%
  mutate(if_else(readStart < 1, 1, readStart+0)) %>%
  filter(if_else(Strand == '+', readStart < Start, readEnd > End)) %>%
  dplyr::select(c(2:8)) %>%
  mutate(readStart = if_else(Strand == '+', readStart-11, readStart+0), readEnd = ifelse(Strand == '+', readEnd+0, readEnd+11), utr5_length = if_else(Strand == '-', (readEnd-End)+1, (Start-readStart)+1)) %>%
  mutate(group_var = paste0(GeneID,';',utr5_length)) %>%
  group_by(group_var) %>%
  arrange(group_var) %>%
  summarise(counts = n()) %>%
  separate(group_var, into = c('GeneID','utr5_length'), sep = ';') %>%
  mutate(utr5_length = as.numeric(utr5_length)) %>% 
  filter(counts >= 10) %>%
  group_by(GeneID) %>%
  filter(counts == max(counts)) %>%
  filter(utr5_length == max(utr5_length)) %>%
  left_join(., subset(fc_gsm, GeneID %in% .$GeneID, select = c(GeneID:Strand))) %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(utr_5 = if_else(Strand == '+', Start - utr5_length, End + utr5_length)) %>%
  mutate(Start2 = if_else(Strand == '+', utr_5, Start), End2 = if_else(Strand == '+', End, utr_5), score = 60) %>%
  filter(Start2 > 0) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'GeneID', 'score','Strand'))

gsm_5UTR_coord <- gsm_all %>%
  mutate(if_else(readStart < 1, 1, readStart+0)) %>%
  filter(if_else(Strand == '+', readStart < Start, readEnd > End)) %>%
  dplyr::select(c(2:8)) %>%
  mutate(readStart = if_else(Strand == '+', readStart-11, readStart+0), readEnd = ifelse(Strand == '+', readEnd+0, readEnd+11), utr5_length = if_else(Strand == '-', (readEnd-End)+1, (Start-readStart)+1)) %>%
  mutate(group_var = paste0(GeneID,';',utr5_length)) %>%
  group_by(group_var) %>%
  arrange(group_var) %>%
  summarise(counts = n()) %>%
  separate(group_var, into = c('GeneID','utr5_length'), sep = ';') %>%
  mutate(utr5_length = as.numeric(utr5_length)) %>% 
  filter(counts >= 10) %>%
  group_by(GeneID) %>%
  filter(counts == max(counts)) %>%
  filter(utr5_length == max(utr5_length)) %>%
  left_join(., subset(fc_gsm, GeneID %in% .$GeneID, select = c(GeneID:Strand))) %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(utr_5 = if_else(Strand == '+', Start - utr5_length, End + utr5_length)) %>%
  mutate(Start2 = if_else(Strand == '+', utr_5, End + 1), End2 = if_else(Strand == '+', Start - 1, utr_5), score = 60) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'GeneID', 'score','Strand')) %>%
  filter(Start2 > 0) %>%
  na.omit()

full_transcr_coord_4gtf <- gsm_3UTR_coord_full_transcr %>%
  left_join(., gsm_5UTR_coord_full_transcr, by = c('GeneID')) %>%
  dplyr::select(-c('score.x','Chr.y','score.y', 'Strand.y')) %>%
  mutate(Start = case_when(is.na(Start2.y == T) ~ Start2.x,
                           Strand.x == '+' & Start2.y < Start2.x ~ Start2.y,
                           Strand.x == '+' & Start2.y > Start2.x ~ Start2.x,
                           Strand.x == '-' ~ Start2.x)) %>%
  mutate(End = case_when(is.na(End2.y == T) ~ End2.x,
                         Strand.x == '+' ~ End2.x,
                         Strand.x == '-' & End2.y > End2.x ~ End2.y,
                         Strand.x == '-' & End2.y < End2.x ~ End2.x)) %>%
  mutate(filter_var = paste0(Chr.x,':',Start,'-',End)) %>%
  distinct(filter_var, .keep_all = T) %>%
  group_by(GeneID) %>%
  mutate(source = 'nanopore', feature = 'transcript', score = '.', frame = '.', tr_counts = 1:n(),
         gene_id = paste0('gene_id "',GeneID,'"'),
         description = paste0(' description "',gsm_ref$description[match(GeneID, gsm_ref$gene_id)],'";'),
         transcript_id = paste0(GeneID,'_',source,'.',tr_counts),
         Transcript_id = paste0(' transcript_id "',GeneID,'_',source,'.',tr_counts,'"')) %>%
  arrange(GeneID) %>%
  dplyr::select(c('GeneID','transcript_id','Chr.x', 'source','feature','Start','End','score','Strand.x','frame','gene_id','Transcript_id','description')) %>%
  unite(V9, gene_id:description, sep = ';', remove = T, na.rm = T) %>%
  na.omit()

utr3_coord_4gtf_pos <- gsm_3UTR_coord %>%
  filter(Strand == '+') %>%
  dplyr::select(-c('score')) %>%
  set_colnames(c('Chr.x','Start','End','GeneID','Strand.x')) %>%
  left_join(.,subset(full_transcr_coord_4gtf, select = -c(Start))) %>%
  mutate(feature = 'three_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

utr3_coord_4gtf_neg <- gsm_3UTR_coord %>%
  filter(Strand == '-') %>%
  dplyr::select(-c('score')) %>%
  set_colnames(c('Chr.x','Start','End','GeneID','Strand.x')) %>%
  left_join(.,subset(full_transcr_coord_4gtf, select = -c(End))) %>%
  mutate(feature = 'three_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

utr5_coord_4gtf_pos <- gsm_5UTR_coord %>%
  filter(Strand == '+') %>%
  dplyr::select(-c('score')) %>%
  set_colnames(c('Chr.x','Start','End','GeneID','Strand.x')) %>%
  left_join(.,subset(full_transcr_coord_4gtf, select = -c(End))) %>%
  mutate(feature = 'five_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

utr5_coord_4gtf_neg <- gsm_5UTR_coord %>%
  filter(Strand == '-') %>%
  dplyr::select(-c('score')) %>%
  set_colnames(c('Chr.x','Start','End','GeneID','Strand.x')) %>%
  left_join(.,subset(full_transcr_coord_4gtf, select = -c(Start))) %>%
  mutate(feature = 'five_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

cds_coord_4gtf <- full_transcr_coord_4gtf %>%
  dplyr::select(-c('Start','End','frame')) %>%
  left_join(.,subset(gsm_ref, primary_tag == 'CDS', select = c(1,4,5,8,10)), by = c('GeneID'='gene_id', 'Chr.x'='seq_id')) %>%
  mutate(length = end - start) %>%
  group_by(transcript_id) %>%
  filter(length == max(length)) %>%
  dplyr::select(c(1,2,3,4,5,9,10,6,7,11,8)) %>%
  mutate(feature = 'CDS') %>%
  dplyr::rename(c('Start' = 'start','End'='end')) %>%
  na.omit()

exon_coord_4gtf <- full_transcr_coord_4gtf %>%
  dplyr::select(-c('Start','End','frame')) %>%
  left_join(.,subset(gsm_ref, primary_tag == 'exon', select = c(1,4,5,8,10)), by = c('GeneID'='gene_id', 'Chr.x'='seq_id')) %>%
  dplyr::select(c(1,2,3,4,5,9,10,6,7,11,8)) %>%
  mutate(feature = 'exon') %>%
  dplyr::rename(c('Start' = 'start','End'='end')) %>%
  na.omit()

gene_coord_4gtf <- full_transcr_coord_4gtf %>%
  mutate(feature = 'gene')

gsm_gtf_new <- rbind(full_transcr_coord_4gtf, utr3_coord_4gtf_pos,utr3_coord_4gtf_neg, utr5_coord_4gtf_pos,utr5_coord_4gtf_neg,cds_coord_4gtf,exon_coord_4gtf,gene_coord_4gtf) %>%
  arrange(Chr.x,Start) %>%
  ungroup() %>%
  dplyr::select(-c('GeneID','transcript_id'))

terms_to_filter <- bind_rows(utr3_coord_4gtf_pos,utr3_coord_4gtf_neg) %>%
  dplyr::select(c(1:3,6,7)) %>%
  left_join(.,gsm_ref_3utr, by = c('GeneID'='GeneID','Chr.x'='Chr')) %>%
  ungroup() %>%
  mutate(length = End - Start, ref_length = utr3_annot_end - utr3_annot_start) %>%
  na.omit() %>%
  filter(abs(ref_length - length) <= 20) %>%
  dplyr::select('transcript_id') %>%
  deframe()

if(length(terms_to_filter) > 0) {
  filtered_gtf <- gsm_gtf_new %>%
    filter(!grepl(paste0(terms_to_filter, collapse="|"), V9))
  if(nrow(filtered_gtf > 0)) {
    fwrite(filtered_gtf, here('data/annotation/gsm_nanopore_filtered.gtf'), sep = '\t', col.names = F, quote = F)
    system("gzip data/annotation/gsm_nanopore_filtered.gtf")
  }
} else {
  fwrite(gsm_gtf_new, here('data/annotation/gsm_nanopore.gtf'), sep = '\t', col.names = F, quote = F)
  system("gzip data/annotation/gsm_nanopore.gtf")
}