### PCA_plots

library('ggfortify')
library('ggplot2')
library('stringr')
library('DESeq2')

### Import of objects

vsd_FAD <- readRDS('archives/vsd_FAD.rds')
vsd_TAU <- readRDS('archives/vsd_TAU.rds')
ALZ_Age <- readRDS('archives/ALZ_Age.rds')
PSP_Age <- readRDS('archives/PSP_Age.rds')

### Import of metadata

FAD_metadata <- read.csv('metadado/FAD_metadata.csv', header = T, 
                         stringsAsFactors = F,
                         row.names =1 )

FAD_metadata$Month <- factor(FAD_metadata$Month, levels = c('4','12','18'))


TAU_metadata <- read.csv('metadado/Tau_metadado.csv',header = T,
                         stringsAsFactors = F,
                         row.names = 1)

mayo_metadado <- read.csv('metadado/mayo_metadado.csv',header = T,
                          stringsAsFactors = F,
                          row.names = 1)
### FAD

FAD_CCX <- vsd_FAD[,FAD_metadata$Region == 'CCX']
FAD_HIP <- vsd_FAD[,FAD_metadata$Region == 'HIP']

pcDat_CCX <- prcomp(t(FAD_CCX), center = T)

autoplot(pcDat_CCX, data = FAD_metadata[FAD_metadata$Region == 'CCX',],
         colour = 'Month' ,
         shape= 'Group', 
         x = 1, y = 3,
         loadings.colour = 'blue',size=4) + 
  scale_x_continuous(limits = c(-0.3, 0.3))+
  scale_y_continuous(limits = c(-0.3, 0.5))+
  theme_bw() + labs(title='FAD5X - Cortex')


pcDat_HIP <- prcomp(t(FAD_HIP),center = T)

autoplot(pcDat_HIP, data = FAD_metadata[FAD_metadata$Region == 'HIP',],
         colour = 'Month' ,
         shape= 'Group', 
         x = 1, y = 2,
         loadings.colour = 'blue',size=4) + 
  scale_x_continuous(limits = c(-0.3, 0.3))+
  scale_y_continuous(limits = c(-0.3, 0.5))+
  theme_bw() + labs(title='FAD5X - Hippocampus')

### TAU

TAU_metadata$Month <- str_replace(TAU_metadata$Month,'_months','')
colnames(TAU_metadata)[colnames(TAU_metadata) == 'mutation'] <- 'Group'
TAU_metadata$Month <- factor(TAU_metadata$Month, levels = c('4','17'))


pcDat_TAU <- prcomp(t(vsd_TAU),center = T)

autoplot(pcDat_TAU, data = TAU_metadata,
         colour = 'Month' ,
         shape= 'Group', 
         x = 1, y = 2,
         loadings.colour = 'blue',size=4) + 
  scale_x_continuous(limits = c(-1, 0.5))+
  scale_y_continuous(limits = c(-1, 0.5))+
  theme_bw() + labs(title='Tau35 - Hippocampus')

### ALZ

ALZ <- vst(ALZ_Age,fitType = 'local',blind = F)
ALZ <- assay(ALZ)

ALZ_metadado <- mayo_metadado[mayo_metadado$Group != 'PSP',]

colnames(ALZ) <- ALZ_metadado$ID

pcaDat_ALZ <- prcomp(t(ALZ),scale = T,center = T)

autoplot(pcaDat_ALZ, data = ALZ_metadado,
         colour = 'Age_Group' ,
         shape= 'Group', 
         x = 1, y = 2,
         loadings.colour = 'blue',size=4) + 
  scale_x_continuous(limits = c(-0.25, 0.25))+
  scale_y_continuous(limits = c(-0.5, 0.5))+
  theme_bw() + labs(title='ALZ - Hippocampus')


### PSP

PSP_metadado <- mayo_metadado[mayo_metadado$Group != 'AD' &
                                mayo_metadado$Age_Group != 'C',]

PSP <- vst(PSP_Age,blind = F,fitType = 'local')
PSP <- assay(PSP)

colnames(PSP) <- PSP_metadado$ID

pcDat_PSP <- prcomp(t(PSP),scale. = T, center = T)

autoplot(pcDat_PSP, data = PSP_metadado,
         colour = 'Age_Group' ,
         shape= 'Group', 
         x = 1, y = 2,
         loadings.colour = 'blue',size=4) + 
  scale_x_continuous(limits = c(-0.25, 0.25))+
  scale_y_continuous(limits = c(-0.5, 0.5))+
  theme_bw() + labs(title='PSP - Hippocampus')

saveRDS(ALZ,'results/vsd_ALZ.rds')
saveRDS(PSP,'results/vsd_PSP.rds')
