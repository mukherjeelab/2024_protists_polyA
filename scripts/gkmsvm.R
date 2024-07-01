library(gkmSVM)
library(caret)
library(tidyverse)
library(magrittr)
library(data.table)
library(here)

set.seed(109)

#gl
system("gunzip data/gkmSVM/gl/*.fa.gz")
gkmsvm_kernel(here('data/gkmSVM/gl/glamb_pos.fa'),
              here('data/gkmSVM/gl/glamb_neg.fa'), 
              here('data/gkmSVM/gl/glamb_kernel.out'), 
              L = 8, addRC = FALSE)
gkmsvm_trainCV(here('data/gkmSVM/gl/glamb_kernel.out'),
               here('data/gkmSVM/gl/glamb_pos.fa'),
               here('data/gkmSVM/gl/glamb_neg.fa'),
               svmfnprfx=here('data/gkmSVM/gl/glamb'), 
               outputCVpredfn=here('data/gkmSVM/gl/glamb.out'), 
               outputROCfn=here('data/gkmSVM/gl/glamb_roc.out'))
gkmsvm_classify(here('data/gkmSVM/8mers.fa'),
                svmfnprfx=here('data/gkmSVM/gl/glamb'),
                here('data/gkmSVM/gl/glamb_8mers.out'), 
                L = 8, addRC = FALSE)
gkmsvm_classify(here('data/gkmSVM/gl/glamb_all.fa'),
                svmfnprfx=here('data/gkmSVM/gl/glamb'), 
                here('data/gkmSVM/gl/glamb_glamb.out'), 
                L = 8, addRC = FALSE)
glamb_pred_glamb <- fread(here('data/gkmSVM/gl/glamb_glamb.out')) %>%
  arrange(desc(V2)) %>%
  mutate(label = if_else(grepl('pos',V1), 1, 0),
         prediction = if_else(V2 > 0, 1, 0))
glamb_pred_conf_mtx_glamb <- confusionMatrix(data = factor(glamb_pred_glamb$prediction), 
                                             reference = factor(glamb_pred_glamb$label),
                                             positive = '1',
                                             mode = 'everything')

#gs
system("gunzip data/gkmSVM/gs/*.fa.gz")
gkmsvm_kernel(here('data/gkmSVM/gs/gsm_pos.fa'),
              here('data/gkmSVM/gs/gsm_neg.fa'), 
              here('data/gkmSVM/gs/gsm_kernel.out'), 
              L = 8, addRC = FALSE)
gkmsvm_trainCV(here('data/gkmSVM/gs/gsm_kernel.out'),
               here('data/gkmSVM/gs/gsm_pos.fa'),
               here('data/gkmSVM/gs/gsm_neg.fa'),
               svmfnprfx=here('data/gkmSVM/gs/gsm'), 
               outputCVpredfn=here('data/gkmSVM/gs/gsm.out'), 
               outputROCfn=here('data/gkmSVM/gs/gsm_roc.out'))
gkmsvm_classify(here('data/gkmSVM/8mers.fa'),
                svmfnprfx=here('data/gkmSVM/gs/gsm'),
                here('data/gkmSVM/gs/gsm_8mers.out'), 
                L = 8, addRC = FALSE)               
gkmsvm_classify(here('data/gkmSVM/gs/gsm_all.fa'),
                svmfnprfx=here('data/gkmSVM/gs/gsm'), 
                here('data/gkmSVM/gs/gsm_gsm.out'), 
                L = 8, addRC = FALSE)               
gsm_pred_gsm <- fread(here('data/gkmSVM/gs/gsm_gsm.out')) %>%
  arrange(desc(V2)) %>%
  mutate(label = if_else(grepl('pos',V1), 1, 0),
         prediction = if_else(V2 > 0, 1, 0))               
gsm_pred_conf_mtx_gsm <- confusionMatrix(data = factor(gsm_pred_gsm$prediction), 
                                         reference = factor(gsm_pred_gsm$label),
                                         positive = '1',
                                         mode = 'everything')
               
#gm
system("gunzip data/gkmSVM/gm/*.fa.gz")
gkmsvm_kernel(here('data/gkmSVM/gm/gmur_pos.fa'),
              here('data/gkmSVM/gm/gmur_neg.fa'), 
              here('data/gkmSVM/gm/gmur_kernel.out'), 
              L = 8, addRC = FALSE)
gkmsvm_trainCV(here('data/gkmSVM/gm/gmur_kernel.out'),
               here('data/gkmSVM/gm/gmur_pos.fa'),
               here('data/gkmSVM/gm/gmur_neg.fa'),
               svmfnprfx=here('data/gkmSVM/gm/gmur'), 
               outputCVpredfn=here('data/gkmSVM/gm/gmur.out'), 
               outputROCfn=here('data/gkmSVM/gm/gmur_roc.out'))
gkmsvm_classify(here('data/gkmSVM/8mers.fa'),
                svmfnprfx=here('data/gkmSVM/gm/gmur'),
                here('data/gkmSVM/gm/gmur_8mers.out'), 
                L = 8, addRC = FALSE)
gkmsvm_classify(here('data/gkmSVM/gm/gmur_all.fa'),
                svmfnprfx=here('data/gkmSVM/gm/gmur'), 
                here('data/gkmSVM/gm/gmur_gmur.out'), 
                L = 8, addRC = FALSE)
gmur_pred_gmur <- fread(here('data/gkmSVM/gm/gmur_gmur.out')) %>%
  arrange(desc(V2)) %>%
  mutate(label = if_else(grepl('pos',V1), 1, 0),
         prediction = if_else(V2 > 0, 1, 0))               
gmur_pred_conf_mtx_gmur <- confusionMatrix(data = factor(gmur_pred_gmur$prediction), 
                                             reference = factor(gmur_pred_gmur$label),
                                             positive = '1',
                                             mode = 'everything')    

#ng (example for sub1)
system("gunzip data/gkmSVM/ng/*.fa.gz")
gkmsvm_kernel(here('data/gkmSVM/ng/neg_pos.fa'),
              here('data/gkmSVM/ng/neg_neg_sub1.fa'), 
              here('data/gkmSVM/ng/neg_kernel.out'), 
              L = 8, addRC = FALSE)
gkmsvm_trainCV(here('data/gkmSVM/ng/neg_kernel.out'),
               here('data/gkmSVM/ng/neg_pos.fa'),
               here('data/gkmSVM/ng/neg_neg_sub1.fa'),
               svmfnprfx=here('data/gkmSVM/ng/neg'), 
               outputCVpredfn=here('data/gkmSVM/ng/neg.out'), 
               outputROCfn=here('data/gkmSVM/ng/neg_roc.out'))
gkmsvm_classify(here('data/gkmSVM/8mers.fa'),
                svmfnprfx=here('data/gkmSVM/ng/neg'),
                here('data/gkmSVM/ng/neg_8mers.out'), 
                L = 8, addRC = FALSE)
gkmsvm_classify(here('data/gkmSVM/ng/neg_all_sub1.fa'),
                svmfnprfx=here('data/gkmSVM/ng/neg'), 
                here('data/gkmSVM/ng/neg_neg.out'), 
                L = 8, addRC = FALSE)
neg_pred_neg <- fread(here('data/gkmSVM/ng/neg_neg.out')) %>%
  arrange(desc(V2)) %>%
  mutate(label = if_else(grepl('pos',V1), 1, 0),
         prediction = if_else(V2 > 0, 1, 0))
neg_pred_conf_mtx_neg <- confusionMatrix(data = factor(neg_pred_neg$prediction), 
                                             reference = factor(neg_pred_neg$label),
                                             positive = '1',
                                             mode = 'everything')

#ng_AAUAAA
system("gunzip data/gkmSVM/ng_AAUAAA/*.fa.gz")
gkmsvm_kernel(here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_pos.fa'),
              here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_neg.fa'), 
              here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_kernel.out'), 
              L = 8, addRC = FALSE)
gkmsvm_trainCV(here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_kernel.out'),
               here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_pos.fa'),
               here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_neg.fa'),
               svmfnprfx=here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA'), 
               outputCVpredfn=here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA.out'), 
               outputROCfn=here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_roc.out'))
gkmsvm_classify(here('data/gkmSVM/8mers.fa'),
                svmfnprfx=here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA'),
                here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_8mers.out'), 
                L = 8, addRC = FALSE)
gkmsvm_classify(here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_all.fa'),
                svmfnprfx=here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA'), 
                here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_neg_AAUAAA.out'), 
                L = 8, addRC = FALSE)
neg_AAUAAA_pred_neg_AAUAAA <- fread(here('data/gkmSVM/ng_AAUAAA/neg_AAUAAA_neg_AAUAAA.out')) %>%
  arrange(desc(V2)) %>%
  mutate(label = if_else(grepl('pos',V1), 1, 0),
         prediction = if_else(V2 > 0, 1, 0))
neg_AAUAAA_pred_conf_mtx_neg_AAUAAA <- confusionMatrix(data = factor(neg_AAUAAA_pred_neg_AAUAAA$prediction), 
                                             reference = factor(neg_AAUAAA_pred_neg_AAUAAA$label),
                                             positive = '1',
                                             mode = 'everything')

#tv (example for sub1)
system("gunzip data/gkmSVM/tv/*.fa.gz")
gkmsvm_kernel(here('data/gkmSVM/tv/trich_pos.fa'),
              here('data/gkmSVM/tv/trich_neg_sub1.fa'), 
              here('data/gkmSVM/tv/trich_kernel.out'), 
              L = 8, addRC = FALSE)
gkmsvm_trainCV(here('data/gkmSVM/tv/trich_kernel.out'),
               here('data/gkmSVM/tv/trich_pos.fa'),
               here('data/gkmSVM/tv/trich_neg_sub1.fa'),
               svmfnprfx=here('data/gkmSVM/tv/trich'), 
               outputCVpredfn=here('data/gkmSVM/tv/trich.out'), 
               outputROCfn=here('data/gkmSVM/tv/trich_roc.out'))
gkmsvm_classify(here('data/gkmSVM/8mers.fa'),
                svmfnprfx=here('data/gkmSVM/tv/trich'),
                here('data/gkmSVM/tv/trich_8mers.out'), 
                L = 8, addRC = FALSE)
gkmsvm_classify(here('data/gkmSVM/tv/trich_all_sub1.fa'),
                svmfnprfx=here('data/gkmSVM/tv/trich'), 
                here('data/gkmSVM/tv/trich_trich.out'), 
                L = 8, addRC = FALSE)
trich_pred_trich <- fread(here('data/gkmSVM/tv/trich_trich.out')) %>%
  arrange(desc(V2)) %>%
  mutate(label = if_else(grepl('pos',V1), 1, 0),
         prediction = if_else(V2 > 0, 1, 0))
trich_pred_conf_mtx_trich <- confusionMatrix(data = factor(trich_pred_trich$prediction), 
                                             reference = factor(trich_pred_trich$label),
                                             positive = '1',
                                             mode = 'everything')

#tf (example for sub1)
system("gunzip data/gkmSVM/tf/*.fa.gz")
gkmsvm_kernel(here('data/gkmSVM/tf/trif_pos.fa'),
              here('data/gkmSVM/tf/trif_neg_sub1.fa'), 
              here('data/gkmSVM/tf/trif_kernel.out'), 
              L = 8, addRC = FALSE)
gkmsvm_trainCV(here('data/gkmSVM/tf/trif_kernel.out'),
               here('data/gkmSVM/tf/trif_pos.fa'),
               here('data/gkmSVM/tf/trif_neg_sub1.fa'),
               svmfnprfx=here('data/gkmSVM/tf/trif'), 
               outputCVpredfn=here('data/gkmSVM/tf/trif.out'), 
               outputROCfn=here('data/gkmSVM/tf/trif_roc.out'))
gkmsvm_classify(here('data/gkmSVM/8mers.fa'),
                svmfnprfx=here('data/gkmSVM/tf/trif'),
                here('data/gkmSVM/tf/trif_8mers.out'), 
                L = 8, addRC = FALSE)
gkmsvm_classify(here('data/gkmSVM/tf/trif_all_sub1.fa'),
                svmfnprfx=here('data/gkmSVM/tf/trif'), 
                here('data/gkmSVM/tf/trif_trif.out'), 
                L = 8, addRC = FALSE)
trif_pred_trif <- fread(here('data/gkmSVM/tf/trif_trif.out')) %>%
  arrange(desc(V2)) %>%
  mutate(label = if_else(grepl('pos',V1), 1, 0),
         prediction = if_else(V2 > 0, 1, 0))
trif_pred_conf_mtx_trif <- confusionMatrix(data = factor(trif_pred_trif$prediction), 
                                             reference = factor(trif_pred_trif$label),
                                             positive = '1',
                                             mode = 'everything')