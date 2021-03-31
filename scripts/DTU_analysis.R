#### Import of Data
library('dplyr')
library('reshape')
library('RColorBrewer')

mayo_A <- readRDS('last_analysis/Age_A.fullAnalysis.rds')
mayo_A <- mayo_A$switchConsequence

mayo_B <- readRDS('last_analysis/Age_B.fullAnalysis.rds')
mayo_B <- mayo_B$switchConsequence

mayo_C <- readRDS('last_analysis/Age_C.fullAnalysis.rds')
mayo_C <- mayo_C$switchConsequence

FAD_4 <- readRDS('last_analysis/Month_4.fullAnalysis.rds')
FAD_4 <- FAD_4$switchConsequence

FAD_12 <- readRDS('last_analysis/Month_12.fullAnalysis.rds')
FAD_12 <- FAD_12$switchConsequence

FAD_18 <- readRDS('last_analysis/Month_18.fullAnalysis.rds')
FAD_18 <- FAD_18$switchConsequence

mayo_DEG_DTU <- readRDS('results/mayo_DEG_DTU.rds')

PA_NDEG <- read.csv('results/PA_NDEG.csv',header = T, stringsAsFactors = F,
                        row.names = 1)

FAD_DEG_DTU <- readRDS('results/FAD_DEG_DTU.rds')

##### Human N DEG

ALZ_NDEG <- mayo_DEG_DTU$ALZ_NDEG
PSP_NDEG <- mayo_DEG_DTU$PSP_NDEG


### 

mayo_A_ALZ <- mayo_A[complete.cases(mayo_A),]
mayo_A_ALZ <- mayo_A_ALZ[mayo_A_ALZ$condition_1 == 'AD' &
  mayo_A_ALZ$condition_2 =='Control',]
mayo_A_ALZ$'Age_Group' <- 'A'

mayo_B_ALZ <- mayo_B[complete.cases(mayo_B),]
mayo_B_ALZ <- mayo_B_ALZ[mayo_B_ALZ$condition_1 == 'AD' &
  mayo_B_ALZ$condition_2 =='Control',]
mayo_B_ALZ$'Age_Group' <- 'B'

mayo_C_ALZ <- mayo_C[complete.cases(mayo_C),]
mayo_C_ALZ <- mayo_C_ALZ[mayo_C_ALZ$condition_1 == 'AD' &
  mayo_C_ALZ$condition_2 =='Control',]
mayo_C_ALZ$'Age_Group' <- 'C'

Alz_all_DTU <- rbind(mayo_A_ALZ,mayo_B_ALZ,mayo_C_ALZ)

names(Alz_all_DTU)[names(Alz_all_DTU) == 'gene_name'] <- 'Gene_Symbol'

ALZ_NDEG_DTU <- merge(Alz_all_DTU,ALZ_NDEG,by='Gene_Symbol')


ALZ_NDEG_DTU <- ALZ_NDEG_DTU[,c(1,3,7,10:12,16,20)]
names(ALZ_NDEG_DTU)[names(ALZ_NDEG_DTU) == 'Age_Group.y'] <- 'Age_Group'


### PA

mayo_A_PA <- mayo_A[complete.cases(mayo_A),]
mayo_A_PA <- mayo_A_PA[mayo_A_PA$condition_1 == 'Control' &
                           mayo_A_PA$condition_2 =='Pathologic_Aging',]
mayo_A_PA$'Age_Group' <- 'A'

mayo_B_PA <- mayo_B[complete.cases(mayo_B),]
mayo_B_PA <- mayo_B_PA[mayo_B_PA$condition_1 == 'Control' &
                           mayo_B_PA$condition_2 =='Pathologic_Aging',]
mayo_B_PA$'Age_Group' <- 'B'

mayo_C_PA <- mayo_C[complete.cases(mayo_C),]
mayo_C_PA <- mayo_C_PA[mayo_C_PA$condition_1 == 'Control' &
                           mayo_C_PA$condition_2 =='Pathologic_Aging',]
mayo_C_PA$'Age_Group' <- 'C'

PA_all_DTU <- rbind(mayo_A_PA,mayo_B_PA,mayo_C_PA)

names(PA_all_DTU)[names(PA_all_DTU) == 'gene_name'] <- 'Gene_Symbol'

PA_NDEG_DTU <- merge(PA_all_DTU,PA_NDEG,by='Gene_Symbol')


PA_NDEG_DTU <- PA_NDEG_DTU[,c(1,3,7,10:12,16,21)]
names(PA_NDEG_DTU)[names(PA_NDEG_DTU) == 'Age_Group.y'] <- 'Age_Group'


### PSP

mayo_A_PSP <- mayo_A[complete.cases(mayo_A),]
mayo_A_PSP <- mayo_A_PSP[mayo_A_PSP$condition_1 == 'Control' &
                         mayo_A_PSP$condition_2 =='PSP',]
mayo_A_PSP$'Age_Group' <- 'A'

mayo_B_PSP <- mayo_B[complete.cases(mayo_B),]
mayo_B_PSP <- mayo_B_PSP[mayo_B_PSP$condition_1 == 'Control' &
                           mayo_B_PSP$condition_2 =='PSP',]
mayo_B_PSP$'Age_Group' <- 'B'

PSP_all_DTU <- rbind(mayo_A_PSP,mayo_B_PSP)

names(PSP_all_DTU)[names(PSP_all_DTU) == 'gene_name'] <- 'Gene_Symbol'

PSP_NDEG_DTU <- merge(PSP_all_DTU,PSP_NDEG,by='Gene_Symbol')


PSP_NDEG_DTU <- PSP_NDEG_DTU[,c(1,4,7,10:12,16,20)]
names(PSP_NDEG_DTU)[names(PSP_NDEG_DTU) == 'Age_Group.y'] <- 'Age_Group'

### Isoform Groups 


### Mouse 

FAD_CCX_NDEG <- FAD_DEG_DTU$FAD_NDEG

FAD_CCX <- FAD_CCX_NDEG[FAD_CCX_NDEG$Area == 'Cortex',]

FAD_CCX <- FAD_CCX_NDEG[FAD_CCX_NDEG$Area == 'Hippocampus',]


FAD_4 <- FAD_4[FAD_4$condition_1 == 'Alzheimer_4_CCX' &
               FAD_4$condition_2 == 'Control_4_CCX' |
               FAD_4$condition_1 == 'Alzheimer_4_HIP' &
               FAD_4$condition_2 == 'Control_4_HIP',]
FAD_4$'Age_Group' <- '4'


FAD_12 <- FAD_12[FAD_12$condition_1 == 'Alzheimer_12_CCX' &
                 FAD_12$condition_2 == 'Control_12_CCX' |
                 FAD_12$condition_1 == 'Alzheimer_12_HIP' &
                 FAD_12$condition_2 == 'Control_12_HIP',]
FAD_12$'Age_Group' <- '12'


FAD_18 <- FAD_18[FAD_18$condition_1 == 'Alzheimer_18_CCX' &
                   FAD_18$condition_2 == 'Control_18_CCX' |
                   FAD_18$condition_1 == 'Alzheimer_18_HIP' &
                   FAD_18$condition_2 == 'Control_18_HIP',]
FAD_18$'Age_Group' <- '18'

FAD_all_DTU <- rbind(FAD_4,FAD_12,FAD_18)
FAD_all_DTU <- FAD_all_DTU[complete.cases(FAD_all_DTU),]

FAD_all_DTU <- FAD_all_DTU %>% mutate(Region = case_when(
  str_detect(condition_1,'CCX') ~ 'Cortex',
  TRUE ~ 'Hippocampus'
))

names(FAD_all_DTU)[names(FAD_all_DTU) == 'gene_name'] <- 'Gene_Symbol'

FAD_NDEG <- rbind(FAD_CCX,FAD_CCX)
FAD_NDEG <- FAD_NDEG[FAD_NDEG$sig !='FDR<0.01',]

FAD_NDEG_DTU <- merge(FAD_all_DTU,FAD_NDEG,by='Gene_Symbol')

FAD_NDEG_DTU <- FAD_NDEG_DTU[,c(1,6,7,10,12,13,14,23,29)]

#####

write.csv(ALZ_NDEG_DTU,'results/isoform_ALZ.csv')
write.csv(PA_NDEG_DTU,'results/isoform_PA.csv')
write.csv(PSP_NDEG_DTU,'results/isoform_PSP.csv')
write.csv(FAD_NDEG_DTU,'results/isoform_FAD.csv')


#### Gráficos


ALZ <- read.csv('results/isoform_ALZ.csv',row.names = 1,stringsAsFactors = F)
ALZ_ABC <- ALZ[ALZ$DTU_column == 'ABC' |
                 ALZ$DTU_column == 'AC'|
                 ALZ$DTU_column == 'C',]

write.csv(ALZ_ABC,'results/tables/ALZ_NDEG_DTU.csv')

ALZ_A <- ALZ[ALZ$Age_Group == 'A',c(6,8)]
ALZ_B <- ALZ[ALZ$Age_Group == 'B',c(6,8)]
ALZ_C <- ALZ[ALZ$Age_Group == 'C',c(6,8)]

ALZ_A <- with(ALZ_A,table(switchConsequence,DTU_column))

ALZ_A <- data.frame(ALZ_A)

ALZ_A <- ALZ_A[ALZ_A$Freq != 0,]

ALZ_A <- data.frame(SWC = ALZ_A$switchConsequence,
                    DTU = ALZ_A$DTU_column,
                    value = ALZ_A$Freq,stringsAsFactors = F, 
                    Age_Group = rep('A',100))

ALZ_B <- with(ALZ_B,table(switchConsequence,DTU_column))

ALZ_B <- data.frame(ALZ_B)

ALZ_B <- ALZ_B[ALZ_B$Freq != 0,]

ALZ_B <- data.frame(SWC = ALZ_B$switchConsequence,
                    DTU = ALZ_B$DTU_column,
                    value = ALZ_B$Freq,stringsAsFactors = F, 
                    Age_Group = rep('B',98))



ALZ_C <- with(ALZ_C,table(switchConsequence,DTU_column))

ALZ_C <- data.frame(ALZ_C)

ALZ_C <- ALZ_C[ALZ_C$Freq != 0,]

ALZ_C <- data.frame(SWC = ALZ_C$switchConsequence,
                    DTU = ALZ_C$DTU_column,
                    value = ALZ_C$Freq,stringsAsFactors = F, 
                    Age_Group = rep('C',100))


ALZ_ABC <- rbind(ALZ_A,ALZ_B,ALZ_C)

ggplot(ALZ_ABC,aes(x=SWC, y=value, fill=DTU))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
   theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
    facet_grid(~Age_Group) 



### PA

PA <- read.csv('results/isoform_PA.csv',row.names = 1,stringsAsFactors = F)


PA_A <- PA[PA$Age_Group == 'A',c(6,8)]
PA_B <- PA[PA$Age_Group == 'B',c(6,8)]
PA_C <- PA[PA$Age_Group == 'C',c(6,8)]

PA_A <- with(PA_A,table(switchConsequence,DTU_column))

PA_A <- data.frame(PA_A)

PA_A <- PA_A[PA_A$Freq != 0,]

PA_A <- data.frame(SWC = PA_A$switchConsequence,
                    DTU = PA_A$DTU_column,
                    value = PA_A$Freq,stringsAsFactors = F, 
                    Age_Group = rep('A',64))

PA_B <- with(PA_B,table(switchConsequence,DTU_column))

PA_B <- data.frame(PA_B)

PA_B <- PA_B[PA_B$Freq != 0,]

PA_B <- data.frame(SWC = PA_B$switchConsequence,
                    DTU = PA_B$DTU_column,
                    value = PA_B$Freq,stringsAsFactors = F, 
                    Age_Group = rep('B',64))



PA_C <- with(PA_C,table(switchConsequence,DTU_column))

PA_C <- data.frame(PA_C)

PA_C <- PA_C[PA_C$Freq != 0,]

PA_C <- data.frame(SWC = PA_C$switchConsequence,
                    DTU = PA_C$DTU_column,
                    value = PA_C$Freq,stringsAsFactors = F, 
                    Age_Group = rep('C',60))


PA_ABC <- rbind(PA_A,PA_B,PA_C)

ggplot(PA_ABC,aes(x=SWC, y=value, fill=DTU))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Age_Group) 



####

PSP <- read.csv('results/isoform_PSP.csv',row.names = 1,stringsAsFactors = F)


PSP_A <- PSP[PSP$Age_Group == 'A',c(6,8)]
PSP_B <- PSP[PSP$Age_Group == 'B',c(6,8)]

PSP_A <- with(PSP_A,table(switchConsequence,DTU_column))

PSP_A <- data.frame(PSP_A)

PSP_A <- PSP_A[PSP_A$Freq != 0,]

PSP_A <- data.frame(SWC = PSP_A$switchConsequence,
                   DTU = PSP_A$DTU_column,
                   value = PSP_A$Freq,stringsAsFactors = F, 
                   Age_Group = rep('A',42))

PSP_B <- with(PSP_B,table(switchConsequence,DTU_column))

PSP_B <- data.frame(PSP_B)

PSP_B <- PSP_B[PSP_B$Freq != 0,]

PSP_B <- data.frame(SWC = PSP_B$switchConsequence,
                   DTU = PSP_B$DTU_column,
                   value = PSP_B$Freq,stringsAsFactors = F, 
                   Age_Group = rep('B',45))


PSP_AB <- rbind(PSP_A,PSP_B)

ggplot(PSP_AB,aes(x=SWC, y=value, fill=DTU))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Age_Group) 

### Mouse





FAD <- read.csv('results/isoform_FAD.csv',header = T, stringsAsFactors = F,
                row.names = 1)

#### Cortex

FAD_CCX <- FAD[FAD$Region == 'Cortex',]
FAD_4_12_18 <- FAD_CCX[FAD_CCX$DTU_column == '4_12_18' |
                       FAD_CCX$DTU_column == '4_18' |
                       FAD_CCX$DTU_column == '18',]

write.csv(FAD_4_12_18,'results/tables/FAD_CCX_4:12:18.csv')

FAD_CCX_4 <- FAD_CCX[FAD_CCX$Month == '4',c(5,9)]
FAD_CCX_12 <- FAD_CCX[FAD_CCX$Month == '12',c(5,9)]
FAD_CCX_18 <- FAD_CCX[FAD_CCX$Month == '18',c(5,9)]


FAD_CCX_4 <- with(FAD_CCX_4,table(switchConsequence,DTU_column))
FAD_CCX_12 <- with(FAD_CCX_12,table(switchConsequence,DTU_column))
FAD_CCX_18 <- with(FAD_CCX_18,table(switchConsequence,DTU_column))

FAD_CCX_4 <- data.frame(FAD_CCX_4)
FAD_CCX_12 <- data.frame(FAD_CCX_12)
FAD_CCX_18 <- data.frame(FAD_CCX_18)

FAD_CCX_4 <- FAD_CCX_4[FAD_CCX_4$Freq != '0',]
FAD_CCX_4$'Age_Group' <- '4'

FAD_CCX_12 <- FAD_CCX_12[FAD_CCX_12$Freq != '0',]
FAD_CCX_12$'Age_Group' <- '12'

FAD_CCX_18 <- FAD_CCX_18[FAD_CCX_18$Freq != '0',]
FAD_CCX_18$'Age_Group' <- '18'

FAD_CCX_SWC <- rbind(FAD_CCX_4,FAD_CCX_12,FAD_CCX_18)

FAD_CCX_SWC <- data.frame(SWC = FAD_CCX_SWC$switchConsequence,
                    DTU = FAD_CCX_SWC$DTU_column,
                    value = FAD_CCX_SWC$Freq,stringsAsFactors = F, 
                    Age_Group = FAD_CCX_SWC$Age_Group)

FAD_CCX_SWC <- FAD_CCX_SWC %>% mutate(DTU = case_when(
  DTU == '4_12_18' ~ '4:12:18',
  DTU == '4_12' ~ '4:12',
  DTU == '12_18' ~ '12:18',
  DTU == '12' ~ '12',
  DTU == '4' ~ '4',
  TRUE ~ '18'
)) 

FAD_CCX_SWC$DTU <- factor(FAD_CCX_SWC$DTU,levels = 
              c('4','4:12','4:12:18','12','12:18','18'))
FAD_CCX_SWC$Age_Group <- factor(FAD_CCX_SWC$Age_Group,levels=c('4','12','18'))


ggplot(FAD_CCX_SWC,aes(x=SWC, y=value, fill=DTU))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Age_Group) 


### Hippocampus

FAD_HIP <- FAD[FAD$Region == 'Hippocampus',]

FAD_4_12_18 <- FAD_HIP[FAD_HIP$DTU_column == '4_12_18' |
                         FAD_HIP$DTU_column == '4_18' |
                         FAD_HIP$DTU_column == '18',]


write.csv(FAD_4_12_18,'results/tables/FAD_HIP_4:12:18.csv')

FAD_HIP_4 <- FAD_HIP[FAD_HIP$Month == '4',c(5,9)]
FAD_HIP_12 <- FAD_HIP[FAD_HIP$Month == '12',c(5,9)]
FAD_HIP_18 <- FAD_HIP[FAD_HIP$Month == '18',c(5,9)]

FAD_HIP_4 <- with(FAD_HIP_4,table(switchConsequence,DTU_column))
FAD_HIP_12 <- with(FAD_HIP_12,table(switchConsequence,DTU_column))
FAD_HIP_18 <- with(FAD_HIP_18,table(switchConsequence,DTU_column))

FAD_HIP_4 <- data.frame(FAD_HIP_4)
FAD_HIP_12 <- data.frame(FAD_HIP_12)
FAD_HIP_18 <- data.frame(FAD_HIP_18)

FAD_HIP_4 <- FAD_HIP_4[FAD_HIP_4$Freq != '0',]
FAD_HIP_4$'Age_Group' <- '4'

FAD_HIP_12 <- FAD_CCX_HIP[HIP_12$Freq != '0',]
FAD_HIP_12$'Age_Group' <- '12'

FAD_HIP_18 <- FAD_HIP_18[FAD_HIP_18$Freq != '0',]
FAD_HIP_18$'Age_Group' <- '18'

FAD_HIP_SWC <- rbind(FAD_HIP_4,FAD_HIP_12,FAD_HIP_18)

FAD_HIP_SWC <- data.frame(SWC = FAD_HIP_SWC$switchConsequence,
                          DTU = FAD_HIP_SWC$DTU_column,
                          value = FAD_HIP_SWC$Freq,stringsAsFactors = F, 
                          Age_Group = FAD_HIP_SWC$Age_Group)

FAD_HIP_SWC <- FAD_HIP_SWC %>% mutate(DTU = case_when(
  DTU == '4_12_18' ~ '4:12:18',
  DTU == '4_12' ~ '4:12',
  DTU == '12_18' ~ '12:18',
  DTU == '12' ~ '12',
  DTU == '4' ~ '4',
  TRUE ~ '18'
)) 

FAD_HIP_SWC$DTU <- factor(FAD_HIP_SWC$DTU,levels = 
                            c('4','4:12','4:12:18','12','12:18','18'))
FAD_HIP_SWC$Age_Group <- factor(FAD_HIP_SWC$Age_Group,levels=c('4','12','18'))


ggplot(FAD_HIP_SWC,aes(x=SWC, y=value, fill=DTU))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Age_Group) 

