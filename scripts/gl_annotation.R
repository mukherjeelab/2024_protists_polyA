library(tidyverse)
library(data.table)
library(magrittr)
library(matrixStats)
library(readxl)
library(here)

glamb_ref <- fread(here('data/annotation/glamb_reference_annot.tsv.gz'))

glamb_ref_3utr <- glamb_ref %>%
  dplyr::select(c('seq_id','primary_tag','start','end','gene_id')) %>%
  filter(primary_tag == 'three_prime_utr')%>%
  dplyr::rename(c('utr3_annot_start'='start','utr3_annot_end'='end','GeneID'='gene_id','Chr'='seq_id')) %>%
  dplyr::select(-c('primary_tag'))

glamb_ref_5utr <- glamb_ref %>%
  dplyr::select(c('seq_id','primary_tag','start','end','gene_id')) %>%
  filter(primary_tag == 'five_prime_utr') %>%
  dplyr::rename(c('utr5_annot_start'='start','utr5_annot_end'='end','GeneID'='gene_id','Chr'='seq_id'))%>%
  dplyr::select(-c('primary_tag'))

fc_glamb <- fread(here('data/featureCounts/glamb_feature_counts_cleaned.tsv.gz')) %>%
  left_join(.,glamb_ref_3utr) %>%
  left_join(.,glamb_ref_5utr) %>%
  mutate(Start = case_when(is.na(utr3_annot_start) == F & Strand == '-' ~ as.integer(utr3_annot_end + 1),
                           is.na(utr5_annot_start) == F & Strand == '+' ~ as.integer(utr5_annot_end + 1),
                           TRUE ~ Start)) %>%
  mutate(End = case_when(is.na(utr3_annot_start) == F & Strand == '+' ~ as.integer(utr3_annot_start - 1),
                         is.na(utr5_annot_start) == F & Strand == '-' ~ as.integer(utr5_annot_start - 1),
                         TRUE ~ End)) %>%
  dplyr::select(-c('utr3_annot_start','utr3_annot_end','utr5_annot_start','utr5_annot_end'))

#rep1
fc_glamb1 <- fc_glamb %>%
  filter(.[,1] > 0) %>%
  dplyr::select(c(1,'GeneID':'Length'))
polyA_glamb1 <- fread(here('data/nanopolish_res/glamb-1_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS' & polya_length >= 30)
annot_glamb1 <- fread(here('data/featureCounts/glamb1_feature_counts_annot_filtered.tsv.gz')) %>%
  filter(V2 == 'Assigned' & V1 %in% polyA_glamb1$readname) %>%
  distinct(.,V1, .keep_all = T) %>%
  dplyr::select(c(1,4)) %>%
  set_colnames(c('readName','GeneID')) %>%
  left_join(.,fc_glamb1)
compare_glamb1 <- fread(here('data/alignment_coords/glamb-1.sorted.bed.gz')) %>%
  distinct(V4, .keep_all = T) %>%
  dplyr::select(-c(5)) %>%
  set_colnames(c('Chr','readStart','readEnd','readName','Strand')) %>% 
  left_join(annot_glamb1,.) %>%
  dplyr::select(c(1,2,4,5,9,6,10,7,8)) %>%
  na.omit() %>%
  filter(if_else(Strand == '+', readEnd > End & readStart < End, readStart < Start & readEnd > Start)) %>%
  mutate(dataset = 'glamb-1')

#rep2
fc_glamb2 <- fc_glamb %>%
  filter(.[,2] > 0) %>%
  dplyr::select(c(1,'GeneID':'Length'))
polyA_glamb2 <- fread(here('data/nanopolish_res/glamb-2_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS' & polya_length >= 30)
annot_glamb2 <- fread(here('data/featureCounts/glamb2_feature_counts_annot_filtered.tsv.gz')) %>%
  filter(V2 == 'Assigned' & V1 %in% polyA_glamb2$readname) %>%
  distinct(.,V1, .keep_all = T) %>%
  dplyr::select(c(1,4)) %>%
  set_colnames(c('readName','GeneID')) %>%
  left_join(.,fc_glamb2)
compare_glamb2 <- fread(here('data/alignment_coords/glamb-2.sorted.bed.gz')) %>%
  distinct(V4, .keep_all = T) %>%
  dplyr::select(-c(5)) %>%
  set_colnames(c('Chr','readStart','readEnd','readName','Strand')) %>% 
  left_join(annot_glamb2,.) %>%
  dplyr::select(c(1,2,4,5,9,6,10,7,8)) %>%
  na.omit() %>%
  filter(if_else(Strand == '+', readEnd > End & readStart < End, readStart < Start & readEnd > Start)) %>%
  mutate(dataset = 'glamb-2')

#rep3
fc_glamb3 <- fc_glamb %>%
  filter(.[,3] > 0) %>%
  dplyr::select(c(1,'GeneID':'Length'))
polyA_glamb3 <- fread(here('data/nanopolish_res/glamb-3_polyA.tsv.gz')) %>%
  filter(qc_tag == 'PASS' & polya_length >= 30)
annot_glamb3 <- fread(here('data/featureCounts/glamb3_feature_counts_annot_filtered.tsv.gz')) %>%
  filter(V2 == 'Assigned' & V1 %in% polyA_glamb3$readname) %>%
  distinct(.,V1, .keep_all = T) %>%
  dplyr::select(c(1,4)) %>%
  set_colnames(c('readName','GeneID')) %>%
  left_join(.,fc_glamb3)
compare_glamb3 <- fread(here('data/alignment_coords/glamb-3.sorted.bed.gz')) %>%
  distinct(V4, .keep_all = T) %>%
  dplyr::select(-c(5)) %>%
  set_colnames(c('Chr','readStart','readEnd','readName','Strand')) %>% 
  left_join(annot_glamb3,.) %>%
  dplyr::select(c(1,2,4,5,9,6,10,7,8)) %>%
  na.omit() %>%
  filter(if_else(Strand == '+', readEnd > End & readStart < End, readStart < Start & readEnd > Start)) %>%
  mutate(dataset = 'glamb-3')

#combining data and calculating isoforms frequencies
glamb_all <- rbind(compare_glamb1, compare_glamb2, compare_glamb3) %>%
  mutate(Start = as.integer(Start), End = as.integer(End)) %>%
  mutate(utr3_length = if_else(Strand == '+', (readEnd-End)+1, (Start-readStart)+1)) %>%
  arrange(GeneID)

glamb_groupped <- glamb_all %>%
  mutate(group_var = paste0(GeneID,';',utr3_length)) %>%
  group_by(group_var) %>%
  arrange(group_var) %>%
  summarise(counts = n()) %>%
  separate(group_var, into = c('GeneID','utr3_length'), sep = ';') %>%
  mutate(utr3_length = as.numeric(utr3_length))

glamb_iso_freq <- glamb_groupped %>% 
  filter(counts >= 10) %>%
  group_by(GeneID) %>%
  arrange(utr3_length, .by_group = T) %>%
  mutate(iso_freq = counts/sum(counts)) %>%
  filter(iso_freq >= 0.1) %>%
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

#calculating coordinates
glamb_3UTR_coord_full_transcr <- glamb_iso_freq %>%
  dplyr::select(c('GeneID','utr3_length')) %>%
  left_join(., subset(fc_glamb, GeneID %in% .$GeneID, select = c(GeneID:Strand))) %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(cl_site = if_else(Strand == '+', End + utr3_length, Start - utr3_length)) %>%
  mutate(Start2 = if_else(Strand == '+', Start, cl_site), End2 = if_else(Strand == '+', cl_site, End), score = 60) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'GeneID', 'score','Strand')) %>%
  filter(Start2 > 0)

glamb_3UTR_coord <- glamb_iso_freq %>%
  dplyr::select(c('GeneID','utr3_length')) %>%
  left_join(., subset(fc_glamb, GeneID %in% .$GeneID, select = c(GeneID:Strand))) %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(cl_site = if_else(Strand == '+', End + utr3_length, Start - utr3_length)) %>%
  mutate(Start2 = if_else(Strand == '+', End + 1, cl_site), End2 = if_else(Strand == '+', cl_site, Start - 1), score = 60) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'GeneID', 'score','Strand')) %>%
  filter(Start2 > 0)

glamb_5UTR_coord_full_transcr <- glamb_all %>%
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
  left_join(., subset(fc_glamb, GeneID %in% .$GeneID, select = c(GeneID:Strand))) %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(utr_5 = if_else(Strand == '+', Start - utr5_length, End + utr5_length)) %>%
  mutate(Start2 = if_else(Strand == '+', utr_5, Start), End2 = if_else(Strand == '+', End, utr_5), score = 60) %>%
  filter(Start2 > 0) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'GeneID', 'score','Strand'))

glamb_5UTR_coord <- glamb_all %>%
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
  left_join(., subset(fc_glamb, GeneID %in% .$GeneID, select = c(GeneID:Strand))) %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(utr_5 = if_else(Strand == '+', Start - utr5_length, End + utr5_length)) %>%
  mutate(Start2 = if_else(Strand == '+', utr_5, End + 1), End2 = if_else(Strand == '+', Start - 1, utr_5), score = 60) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'GeneID', 'score','Strand')) %>%
  filter(Start2 > 0) %>%
  na.omit()

full_transcr_coord_4gtf <- glamb_3UTR_coord_full_transcr %>%
  left_join(., glamb_5UTR_coord_full_transcr, by = c('GeneID')) %>%
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
         description = paste0(' description "',glamb_ref$description[match(GeneID, glamb_ref$gene_id)],'";'),
         transcript_id = paste0(GeneID,'_',source,'.',tr_counts),
         Transcript_id = paste0(' transcript_id "',GeneID,'_',source,'.',tr_counts,'"')) %>%
  arrange(GeneID) %>%
  dplyr::select(c('GeneID','transcript_id','Chr.x', 'source','feature','Start','End','score','Strand.x','frame','gene_id','Transcript_id','description')) %>%
  unite(V9, gene_id:description, sep = ';', remove = T, na.rm = T) %>%
  na.omit()

utr3_coord_4gtf_pos <- glamb_3UTR_coord %>%
  filter(Strand == '+') %>%
  dplyr::select(-c('score')) %>%
  set_colnames(c('Chr.x','Start','End','GeneID','Strand.x')) %>%
  left_join(.,subset(full_transcr_coord_4gtf, select = -c(Start))) %>%
  mutate(feature = 'three_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

utr3_coord_4gtf_neg <- glamb_3UTR_coord %>%
  filter(Strand == '-') %>%
  dplyr::select(-c('score')) %>%
  set_colnames(c('Chr.x','Start','End','GeneID','Strand.x')) %>%
  left_join(.,subset(full_transcr_coord_4gtf, select = -c(End))) %>%
  mutate(feature = 'three_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

utr5_coord_4gtf_pos <- glamb_5UTR_coord %>%
  filter(Strand == '+') %>%
  dplyr::select(-c('score')) %>%
  set_colnames(c('Chr.x','Start','End','GeneID','Strand.x')) %>%
  left_join(.,subset(full_transcr_coord_4gtf, select = -c(End))) %>%
  mutate(feature = 'five_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

utr5_coord_4gtf_neg <- glamb_5UTR_coord %>%
  filter(Strand == '-') %>%
  dplyr::select(-c('score')) %>%
  set_colnames(c('Chr.x','Start','End','GeneID','Strand.x')) %>%
  left_join(.,subset(full_transcr_coord_4gtf, select = -c(Start))) %>%
  mutate(feature = 'five_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

cds_coord_4gtf <- full_transcr_coord_4gtf %>%
  dplyr::select(-c('Start','End','frame')) %>%
  left_join(.,subset(glamb_ref, primary_tag == 'CDS', select = c(1,4,5,8,10)), by = c('GeneID'='gene_id', 'Chr.x'='seq_id')) %>%
  mutate(length = end - start) %>%
  group_by(transcript_id) %>%
  filter(length == max(length)) %>%
  dplyr::select(c(1,2,3,4,5,9,10,6,7,11,8)) %>%
  mutate(feature = 'CDS') %>%
  dplyr::rename(c('Start' = 'start','End'='end')) %>%
  na.omit()

exon_coord_4gtf <- full_transcr_coord_4gtf %>%
  dplyr::select(-c('Start','End','frame')) %>%
  left_join(.,subset(glamb_ref, primary_tag == 'exon', select = c(1,4,5,8,10)), by = c('GeneID'='gene_id', 'Chr.x'='seq_id')) %>%
  dplyr::select(c(1,2,3,4,5,9,10,6,7,11,8)) %>%
  mutate(feature = 'exon') %>%
  dplyr::rename(c('Start' = 'start','End'='end')) %>%
  na.omit()

gene_coord_4gtf <- full_transcr_coord_4gtf %>%
  mutate(feature = 'gene')

glamb_gtf_new <- rbind(full_transcr_coord_4gtf, utr3_coord_4gtf_pos,utr3_coord_4gtf_neg, utr5_coord_4gtf_pos,utr5_coord_4gtf_neg,cds_coord_4gtf,exon_coord_4gtf,gene_coord_4gtf) %>%
  arrange(Chr.x,Start) %>%
  ungroup() %>%
  dplyr::select(-c('GeneID','transcript_id'))

#adding data from previous paper
singles <- read_xlsx(here('data/annotation/Supplemental_table_3_UTRLengths.xlsx'), sheet = 'Singles') %>%
  dplyr::select(-c(3)) %>%
  arrange(gene_id)
doubles <- read_xlsx(here('data/annotation/Supplemental_table_3_UTRLengths.xlsx'), sheet = 'Doubles') %>%
  pivot_longer(cols = 2:3, names_to = 'temp', values_to = 'utr_length') %>%
  arrange(gene_id) %>%
  dplyr::select(-c(2,3))
triple <- read_xlsx(here('data/annotation/Supplemental_table_3_UTRLengths.xlsx'), sheet = 'Triple') %>%
  pivot_longer(cols = 2:4, names_to = 'temp', values_to = 'utr_length') %>%
  arrange(gene_id) %>%
  dplyr::select(-c(2,3))

utr_3_Danielle_full_transcr <- rbind(singles,doubles,triple) %>% 
  left_join(., subset(fc_glamb, GeneID %in% .$gene_id, select = c(4:8)), by = c('gene_id' = 'GeneID'))  %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(cl_site = if_else(Strand == '+', End + utr_length, Start - utr_length)) %>%
  mutate(Start2 = if_else(Strand == '+', Start, cl_site), End2 = if_else(Strand == '+', cl_site, End), score = 60) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'gene_id', 'score','Strand')) %>%
  filter(Start2 > 0)

utr_3_Danielle <- rbind(singles,doubles,triple) %>% 
  left_join(., subset(fc_glamb, GeneID %in% .$gene_id, select = c(4:8)), by = c('gene_id' = 'GeneID'))  %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  mutate(cl_site = if_else(Strand == '+', End + utr_length, Start - utr_length)) %>%
  mutate(Start2 = if_else(Strand == '+', End + 1, cl_site), End2 = if_else(Strand == '+', cl_site, Start - 1), score = 60) %>%
  dplyr::select(c('Chr','Start2', 'End2', 'gene_id', 'score','Strand')) %>%
  filter(Start2 > 0)

full_transcr_Danielle_4gtf <- utr_3_Danielle_full_transcr %>%
  dplyr::rename('GeneID'='gene_id') %>%
  mutate(source = 'QuantSeq_and_nanopore', feature = 'transcript', score = '.', frame = '.', tr_counts = 1:n(),
         gene_id = paste0('gene_id "',GeneID,'"'),
         description = paste0(' description "',glamb_ref$description[match(GeneID, glamb_ref$gene_id)],'";'),
         transcript_id = paste0(GeneID,'_',source,'.',tr_counts),
         Transcript_id = paste0(' transcript_id "',GeneID,'_',source,'.',tr_counts,'"')) %>%
  arrange(GeneID) %>%
  dplyr::select(c('GeneID','transcript_id','Chr', 'source','feature','Start2','End2','score','Strand','frame','gene_id','Transcript_id','description')) %>%
  unite(V9, gene_id:description, sep = ';', remove = T, na.rm = T) %>%
  na.omit()

utr_3_Danielle_4gtf_pos <- utr_3_Danielle %>%
  filter(Strand == '+') %>%
  dplyr::select(-c('score')) %>%
  dplyr::rename('GeneID'='gene_id') %>%
  left_join(.,subset(full_transcr_Danielle_4gtf, select = -c(Start2))) %>%
  mutate(feature = 'three_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

utr_3_Danielle_4gtf_neg <- utr_3_Danielle %>%
  filter(Strand == '-') %>%
  dplyr::select(-c('score')) %>%
  dplyr::rename('GeneID'='gene_id') %>%
  left_join(.,subset(full_transcr_Danielle_4gtf, select = -c(End2))) %>%
  mutate(feature = 'three_prime_utr') %>%
  dplyr::select(4,6,1,7,8,2,3,9,5,10,11) %>%
  na.omit()

cds_coord_Danielle_4gtf <- full_transcr_Danielle_4gtf %>%
  dplyr::select(-c('Start2','End2','frame')) %>%
  left_join(.,subset(glamb_ref, primary_tag == 'CDS', select = c(1,4,5,8,10)), by = c('GeneID'='gene_id', 'Chr'='seq_id')) %>%
  mutate(length = end - start) %>%
  group_by(transcript_id) %>%
  filter(length == max(length)) %>%
  dplyr::select(c(1,2,3,4,5,9,10,6,7,11,8)) %>%
  mutate(feature = 'CDS') %>%
  dplyr::rename(c('Start2' = 'start','End2'='end')) %>%
  mutate_if(is.integer, as.numeric) %>%
  na.omit()

exon_coord_Danielle_4gtf <- full_transcr_Danielle_4gtf %>%
  dplyr::select(-c('Start2','End2','frame')) %>%
  left_join(.,subset(glamb_ref, primary_tag == 'exon', select = c(1,4,5,8,10)), by = c('GeneID'='gene_id', 'Chr'='seq_id')) %>%
  dplyr::select(c(1,2,3,4,5,9,10,6,7,11,8)) %>%
  mutate(feature = 'exon') %>%
  dplyr::rename(c('Start2' = 'start','End2'='end')) %>%
  mutate_if(is.integer, as.numeric) %>%
  na.omit()

gene_coord_Danielle_4gtf <- full_transcr_Danielle_4gtf %>%
  mutate(feature = 'gene')

glamb_gtf_new_Danielle <- bind_rows(full_transcr_Danielle_4gtf, utr_3_Danielle_4gtf_pos,utr_3_Danielle_4gtf_neg, cds_coord_Danielle_4gtf,exon_coord_Danielle_4gtf, gene_coord_Danielle_4gtf) %>%
  dplyr::arrange(Chr,Start2) %>%
  ungroup() %>%
  dplyr::select(-c('GeneID','transcript_id'))

singles2remove <- glamb_iso_freq %>%
  filter(category == 'singles') %>%
  filter(GeneID%in%full_transcr_Danielle_4gtf$GeneID) %>%
  mutate(transcript_id = full_transcr_coord_4gtf$transcript_id[match(GeneID, full_transcr_coord_4gtf$GeneID)]) %>%
  dplyr::select(c('transcript_id')) %>%
  deframe()

doubles2remove <- glamb_iso_freq %>%
  filter(category == 'doubles') %>%
  filter(GeneID%in%doubles$gene_id) %>%
  mutate(transcript_id = paste0(GeneID,'_nanopore.',cl_site_n)) %>%
  dplyr::select(c('transcript_id')) %>%
  deframe()

doubles2selectiveremove <- glamb_iso_freq %>%
  filter(category == 'doubles' &!GeneID%in%names(doubles2remove)) %>%
  left_join(., singles, by = c('GeneID'='gene_id')) %>%
  na.omit() %>%
  group_by(GeneID) %>%
  mutate(transcript_id = paste0(GeneID,'_nanopore.',cl_site_n), abs_diff = abs(utr3_length - utr_length)) %>%
  filter(abs_diff == max(abs_diff)) %>%
  dplyr::select(c('transcript_id')) %>%
  deframe()

lines2remove <- c(singles2remove, doubles2remove, doubles2selectiveremove)

full_transcr_coord_4gtf_comb <- full_transcr_Danielle_4gtf %>%
  left_join(., glamb_5UTR_coord_full_transcr, by = c('GeneID')) %>%
  dplyr::select(-c('score.x','Chr.y','score.y', 'Strand.y')) %>%
  mutate(Start = case_when(is.na(Start2.y) == T ~ Start2.x,
                           Strand.x == '+' & Start2.y < Start2.x ~ Start2.y,
                           Strand.x == '+' & Start2.y > Start2.x ~ Start2.x,
                           Strand.x == '-' ~ Start2.x)) %>%
  mutate(End = case_when(is.na(End2.y) == T ~ End2.x,
                         Strand.x == '+' ~ End2.x,
                         Strand.x == '-' & End2.y > End2.x ~ End2.y,
                         Strand.x == '-' & End2.y < End2.x ~ End2.x)) %>%
  mutate(score = '.') %>%
  dplyr::select(1:5,13:15,8:10) %>%
  set_colnames(colnames(full_transcr_Danielle_4gtf)) %>%
  na.omit() %>%
  dplyr::select(-c(1,2))

utr5_coord_4gtf_comb <- bind_rows(utr5_coord_4gtf_pos, utr5_coord_4gtf_neg) %>%
  mutate(source = if_else(GeneID%in%full_transcr_Danielle_4gtf$GeneID,'QuantSeq_and_nanopore',source)) %>%
  set_colnames(colnames(full_transcr_Danielle_4gtf)) %>%
  na.omit() %>%
  ungroup() %>%
  dplyr::select(-c(1,2))

qs_nanopore_gtf <- glamb_gtf_new_Danielle %>%
  filter(!feature == 'transcript')

nanopore_gtf <- glamb_gtf_new %>%
  filter(!feature == 'five_prime_utr') %>%
  filter(!grepl(paste0(lines2remove, collapse="|"), V9)) %>%
  set_colnames(colnames(full_transcr_Danielle_4gtf)[3:11])

combined_gtf <- bind_rows(full_transcr_coord_4gtf_comb, utr5_coord_4gtf_comb, qs_nanopore_gtf, nanopore_gtf) %>%
  group_by(source) %>%
  dplyr::arrange(Chr,Start2,source)

fwrite(combined_gtf, here('data/annotation/glamb_combined.gtf'), sep = '\t', col.names = F, quote = F)
system("gzip data/annotation/glamb_combined.gtf")