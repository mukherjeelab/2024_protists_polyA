library(tidyverse)
library(rtracklayer)
library(ggcoverage)
library(ggpattern)
library(plyranges)
library(Biostrings)
library(here)

#4B
cds_cts <- read_tsv(here("data/prem_cl/cds_counts.tsv.gz"))
pdf(here("figs/4B_raw.pdf"))
ggbarplot(cds_cts, x = 'dataset', y = 'counts', fill = 'gray') +
  rremove("xlab")
dev.off()
#4F
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

#S3F
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