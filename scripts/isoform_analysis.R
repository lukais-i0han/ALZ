### DTU analysis
library('IsoformSwitchAnalyzeR')
library('dplyr')
library('stringr')
## Human


mayo_A_ISW <- readRDS('archives/ISW/MAYO/Age_A.fullAnalysis.rds')
stage_ALZ_A <- read.csv('stage_archives/results/Human/stage_ALZ_A.csv',header = T,
              stringsAsFactors = F,row.names = 1)
stage_PSP_A <- read.csv('stage_archives/results/Human/stage_PSP_A.csv',header = T,
                        stringsAsFactors = F,row.names = 1)

mayo_B_ISW <- readRDS('archives/ISW/MAYO/Age_B.fullAnalysis.rds')
stage_ALZ_B <- read.csv('stage_archives/results/Human/stage_ALZ_B.csv',header = T,
                stringsAsFactors = F, row.names = 1)
stage_PSP_B <- read.csv('stage_archives/results/Human/stage_PSP_B.csv',header = T,
                        stringsAsFactors = F,row.names = 1)


mayo_C_ISW <- readRDS('archives/ISW/MAYO/Age_C.fullAnalysis.rds')
stage_ALZ_C <- read.csv('stage_archives/results/Human/stage_ALZ_C.csv',header = T,
                stringsAsFactors = F,row.names = 1)

## Human
mayo_A_features <- mayo_A_ISW$isoformFeatures

mayo_A_analysis <- mayo_A_ISW$isoformSwitchAnalysis
mayo_A_analysis <- mayo_A_analysis[order(mayo_A_analysis$iso_ref),]

mayo_A_features$gene_q_value <- mayo_A_analysis$gene_q_value
mayo_A_features$iso_q_value <- mayo_A_analysis$iso_q_value
mayo_A <- mayo_A_features[,c(3:9,16,22:23,27:28)]
mayo_A$'Age_Group' <- 'A' 

mayo_B_features <- mayo_B_ISW$isoformFeatures

mayo_B_analysis <- mayo_B_ISW$isoformSwitchAnalysis
mayo_B_analysis <- mayo_B_analysis[order(mayo_B_analysis$iso_ref),]

mayo_B_features$gene_q_value <- mayo_B_analysis$gene_q_value
mayo_B_features$iso_q_value <- mayo_B_analysis$iso_q_value
mayo_B <- mayo_B_features[,c(3:9,16,22:23,27:28)]
mayo_B$'Age_Group' <- 'B'

mayo_C_features <- mayo_C_ISW$isoformFeatures

mayo_C_analysis <- mayo_C_ISW$isoformSwitchAnalysis
mayo_C_analysis <- mayo_C_analysis[order(mayo_C_analysis$iso_ref),]

mayo_C_features$gene_q_value <- mayo_C_analysis$gene_q_value
mayo_C_features$iso_q_value <- mayo_C_analysis$iso_q_value
mayo_C <- mayo_C_features[,c(3:9,16,22:23,27:28)]
mayo_C$'Age_Group' <- 'C'

mayo_ALZ <- rbind(mayo_A,mayo_B,mayo_C)

mayo_ALZ <- mayo_ALZ[mayo_ALZ$condition_1 == 'AD' &
                       mayo_ALZ$condition_2 == 'Control',]

stage_ALZ <- c(stage_ALZ_A$txID,stage_ALZ_B$txID,stage_ALZ_C$txID)

mayo_ALZ <- mayo_ALZ[mayo_ALZ$isoform_id %in% stage_ALZ,]


mayo_PSP <- rbind(mayo_A,mayo_B)
mayo_PSP <- mayo_PSP[mayo_PSP$condition_1 == 'Control' &
                     mayo_PSP$condition_2 == 'PSP',]
stage_PSP <- c(stage_PSP_A$txID,stage_PSP_B$txID)
mayo_PSP <- mayo_PSP[mayo_PSP$isoform_id %in% stage_PSP,]



mayo_ISW_PA <- mayo_ISW[mayo_ISW$condition_2 == 'Pathologic_Aging',]
mayo_ISW_PA$dIF <- mayo_ISW_PA$dIF*(-1)
mayo_ISW_PA$condition_1 <- 'PA'
mayo_ISW_PA$condition_2 <- 'Control'


mayo_ISW <- rbind(mayo_ALZ,mayo_PSP)
mayo_ISW <- mayo_ISW  %>% mutate(DTU = case_when(
  abs(mayo_ISW$dIF) > 0.05 & mayo_ISW$iso_q_value < 0.01 ~ 'DTU',
  TRUE ~ 'Non-DTU'
))


mayo_ISW_sig <- mayo_ISW[mayo_ISW$DTU == 'DTU',]


saveRDS(mayo_ISW,'results/isoform_results/mayo_ISW.rds')
saveRDS(mayo_ISW_sig,'results/isoform_results/mayo_ISW_sig.rds')

### Mouse - FAD5X

month_4 <- readRDS('archives/ISW/5XFAD/Month_4.fullAnalysis.rds')

month_12 <- readRDS('archives/ISW/5XFAD/Month_12.fullAnalysis.rds')

month_18 <- readRDS('archives/ISW/5XFAD/Month_18.fullAnalysis.rds')

month_4_features <- month_4$isoformFeatures

month_4_analysis <- month_4$isoformSwitchAnalysis
month_4_analysis <- month_4_analysis[order(month_4_analysis$iso_ref),]

month_4_features$gene_q_value <- month_4_analysis$gene_q_value
month_4_features$iso_q_value <- month_4_analysis$iso_q_value

month_4 <- month_4_features[,c(3:9,16,22:23,27:28)]
month_4 <- month_4[month_4$condition_1 == 'Alzheimer_4_CCX' & month_4$condition_2 == 'Control_4_CCX'|
          month_4$condition_1 == 'Alzheimer_4_HIP' & month_4$condition_2 == 'Control_4_HIP',]

month_4$'Age_Group' <- '4'

month_12_features <- month_12$isoformFeatures

month_12_analysis <- month_12$isoformSwitchAnalysis
month_12_analysis <- month_12_analysis[order(month_12_analysis$iso_ref),]

month_12_features$gene_q_value <- month_12_analysis$gene_q_value
month_12_features$iso_q_value <- month_12_analysis$iso_q_value

month_12 <- month_12_features[,c(3:9,16,22:23,27:28)]
month_12 <- month_12[month_12$condition_1 == 'Alzheimer_12_CCX' & month_12$condition_2 == 'Control_12_CCX'|
                 month_12$condition_1 == 'Alzheimer_12_HIP' & month_12$condition_2 == 'Control_12_HIP',]
month_12$'Age_Group' <- '12'


month_18_features <- month_18$isoformFeatures

month_18_analysis <- month_18$isoformSwitchAnalysis
month_18_analysis <- month_18_analysis[order(month_18_analysis$iso_ref),]

month_18_features <- month_18_features[month_18_features$iso_ref %in% 
                      month_18_analysis$iso_ref,]

month_18_features$gene_q_value <- month_18_analysis$gene_q_value
month_18_features$iso_q_value <- month_18_analysis$iso_q_value

month_18 <- month_18_features[,c(3:9,16,22:23,27:28)]
month_18 <- month_18[month_18$condition_1 == 'Alzheimer_18_CCX' & month_18$condition_2 == 'Control_18_CCX'|
            month_18$condition_1 == 'Alzheimer_18_HIP' & month_18$condition_2 == 'Control_18_HIP',]
month_18$'Age_Group' <- '18'

FAD5X <- rbind(month_4,month_12,month_18)

FAD5X <- FAD5X %>% mutate(DTU = case_when(
  abs(FAD5X$dIF) > 0.05 & FAD5X$iso_q_value < 0.01 ~ 'DTU',
  TRUE ~ 'Non-DTU'
))

FAD_CCX_4 <- read.csv('stage_archives/results/Mouse/stage_FAD_CCX_4.csv',header = T,
             stringsAsFactors = F, row.names = 1)
FAD_CCX_12 <- read.csv('stage_archives/results/Mouse/stage_FAD_CCX_12.csv',header = T,
                      stringsAsFactors = F, row.names = 1)
FAD_CCX_18 <- read.csv('stage_archives/results/Mouse/stage_FAD_CCX_18.csv',header = T,
                      stringsAsFactors = F, row.names = 1)

FAD_HIP_4 <- read.csv('stage_archives/results/Mouse/stage_FAD_HIP_4.csv',header = T,
                      stringsAsFactors = F, row.names = 1)
FAD_HIP_12 <- read.csv('stage_archives/results/Mouse/stage_FAD_HIP_12.csv',header = T,
                       stringsAsFactors = F)
FAD_HIP_18 <- read.csv('stage_archives/results/Mouse/stage_FAD_HIP_18.csv',header = T,
                       stringsAsFactors = F, row.names = 1)



FAD_stage <- c(FAD_CCX_4$txID,FAD_CCX_12$txID,FAD_CCX_18$txID,
               FAD_HIP_4$txID,FAD_HIP_12$txID,FAD_HIP_18$txID)

FAD5X$isoform_id <- str_sub(FAD5X$isoform_id,start = 1,end = 18)

FAD5X <- FAD5X[FAD5X$isoform_id %in% FAD_stage,]

FAD5X <- FAD5X %>% mutate(Region = case_when(
  str_detect(FAD5X$condition_1,'_CCX') ~ 'Cortex',
  TRUE ~ 'Hippocampus'
))


FAD5X_sig <- FAD5X[FAD5X$DTU == 'DTU',]

saveRDS(FAD5X,'results/isoform_results/FAD5X_ISW.rds')
saveRDS(FAD5X_sig,'results/isoform_results/FAD5X_ISW_sig.rds')

### Mouse - Tau

month_4 <- readRDS('archives/ISW/TAU/age_4.fullAnalysis.rds')

month_17 <- readRDS('archives/ISW/TAU/age_17.fullAnalysis.rds')

month_4_features <- month_4$isoformFeatures

month_4_analysis <- month_4$isoformSwitchAnalysis
month_4_analysis <- month_4_analysis[order(month_4_analysis$iso_ref),]

month_4_features$iso_q_value <- month_4_analysis$padj

month_4 <- month_4_features[,c(3:9,16,22:23,27:28)]
month_4$'Age_Group' <- '4'

month_17_features <- month_17$isoformFeatures

month_17_analysis <- month_17$isoformSwitchAnalysis
month_17_analysis <- month_17_analysis[order(month_17_analysis$iso_ref),]

month_17_features$iso_q_value <- month_17_analysis$padj

month_17 <- month_17_features[,c(3:9,16,22:23,27:28)]
month_17$'Age_Group' <- '17'

Tau <- rbind(month_4,month_17)

stage_4 <-Tau_4 = read.csv('stage_archives/results/Mouse/stage_TAU_4.csv',
                                   header = T, stringsAsFactors = F, row.names = 1)
                  

stage_17$transcript <- gsub('"',replacement = '',stage_17$transcript) ### importei direto pelo R

stage_TAU <- c(stage_TAU$Tau_4$txID,stage_17$transcript)

Tau <- Tau[Tau$isoform_id %in% stage_TAU,]

Tau <- Tau %>% mutate(DTU = case_when(
  abs(Tau$dIF) > 0.05 & iso_q_value < 0.01 ~ 'DTU',
  TRUE ~ 'Non-DTU'
  ))

Tau_sig <- Tau[Tau$DTU == 'DTU',]

saveRDS(Tau,'results/isoform_results/Tau_ISW.rds')
saveRDS(Tau_sig,'results/isoform_results/Tau_ISW_sig.rds')
