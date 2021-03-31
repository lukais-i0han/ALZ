library('stageR')
library('DEXSeq')
library('dplyr')
library('stringr')



Age_Group <- c('4','12','18')

tx2gene <- read.table('stage_archives/refs/Mouse/t2g.txt',header = T)

counts_FAD <- read.csv('stage_archives/counts_FAD.csv',header = T,
                       stringsAsFactors = F,row.names = 1,check.names=FALSE)

metadado_FAD <- read.csv('stage_archives/refs/Mouse/metadado_FAD.csv',header = T, stringsAsFactors = F,
                         row.names = 1)

colnames(counts_FAD) <- metadado_FAD$SpecimenID

FAD4 <- readRDS('stage_archives/refs/Mouse/Month_4.fullAnalysis.rds')
FAD4 <- FAD4$isoformFeatures

FAD4_CCX <- FAD4[FAD4$condition_1 =='Alzheimer_4_CCX' & FAD4$condition_2 == 'Control_4_CCX',]
FAD4_HIP <- FAD4[FAD4$condition_1 =='Alzheimer_4_HIP' & FAD4$condition_2 == 'Control_4_HIP',]
            
FAD12 <- readRDS('stage_archives/refs/Mouse/Month_12.fullAnalysis.rds')
FAD12 <- FAD12$isoformFeatures

FAD12_CCX <- FAD12[FAD12$condition_1 =='Alzheimer_12_CCX' & FAD12$condition_2 == 'Control_12_CCX',]
FAD12_HIP <- FAD12[FAD12$condition_1 =='Alzheimer_12_HIP' & FAD12$condition_2 == 'Control_12_HIP',]


FAD18 <- readRDS('stage_archives/refs/Mouse/Month_18.fullAnalysis.rds')
FAD18 <- FAD18$isoformFeatures

FAD18_CCX <- FAD18[FAD18$condition_1 =='Alzheimer_18_CCX' & FAD18$condition_2 == 'Control_18_CCX',]
FAD18_HIP <- FAD18[FAD18$condition_1 =='Alzheimer_18_HIP' & FAD18$condition_2 == 'Control_18_HIP',]


FAD_all <- rbind(FAD4_CCX,FAD4_HIP,FAD12_CCX,FAD12_HIP,FAD18_CCX,FAD18_HIP)

FAD_all <- FAD_all %>% mutate(Group = case_when(
  str_detect(condition_1,'_4_CCX') ~ 'CCX_4',
  str_detect(condition_1,'_12_CCX') ~ 'CCX_12',
  str_detect(condition_1,'_18_CCX') ~ 'CCX_18',
  str_detect(condition_1,'_4_HIP') ~ 'HIP_4',
  str_detect(condition_1,'_12_HIP') ~ 'HIP_12',
  TRUE ~ 'HIP_18'
))

metadado_FAD$'Group' <- paste(metadado_FAD$Region,metadado_FAD$Month,sep='_')

saveRDS(FAD_all,'stage_archives/refs/Mouse/FAD_all.rds')
write.csv(metadado_FAD,'stage_archives/metadado_FAD.csv')


