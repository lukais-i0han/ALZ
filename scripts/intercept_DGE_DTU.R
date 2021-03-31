### Interception of DEG and DTU
library('stringr')
library('dplyr')
library("UpSetR")
library('tidyr')

### Human DGE

ALZ_sig <- read.csv('archives/ALZ_sig.csv',header = T, stringsAsFactors = F,
                 row.names = 1)

names(ALZ_sig)[names(ALZ_sig) == 'AgeAtDeath'] <- 'Age_Group'


PSP_sig <- read.csv('archives/PSP_sig.csv',header = T, stringsAsFactors = F,
                row.names = 1)

ALZ <- list(ALZ_A = unique(ALZ_A$Gene_Symbol),
            ALZ_B = unique(ALZ_B$Gene_Symbol),
            ALZ_C = unique(ALZ_C$Gene_Symbol))

intercept_ALZ <- venn::venn(ALZ,
                  zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)


PSP <- list(PSP_A = unique(PSP_A$Gene_Symbol),
            PSP_B = unique(PSP_B$Gene_Symbol))

intercept_PSP <- venn::venn(PSP,
  zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)




ALZ_A <- ALZ_sig[ALZ_sig$Age_Group == 'A',]
ALZ_B <- ALZ_sig[ALZ_sig$Age_Group == 'B',]
ALZ_C <- ALZ_sig[ALZ_sig$Age_Group == 'C',]

PSP_A <- PSP_sig[PSP_sig$Age_Group == 'A',]
PSP_B <- PSP_sig[PSP_sig$Age_Group == 'B',]


ALZ_PSP <- list(ALZ_A = unique(ALZ_A$Gene_Symbol),
                  ALZ_B = unique(ALZ_B$Gene_Symbol),
                  ALZ_C = unique(ALZ_C$Gene_Symbol),
                  PSP_A = unique(PSP_A$Gene_Symbol),
                  PSP_B = unique(PSP_B$Gene_Symbol))

intercept_ALZ_PSP <- venn::venn(ALZ_PSP,
      zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)
inters_ALZ_PSP <-attr(intercept_ALZ_PSP,"intersections")



upset(fromList(ALZ_PSP),sets = c('ALZ_A','ALZ_B','ALZ_C',
                                 'PSP_A','PSP_B'),order.by = c('freq','degree'))

inters <- plyr::ldply(inters_ALZ_PSP, rbind)
inters <- t(inters)
inters <- data.frame(inters,stringsAsFactors = F)
colnames(inters) <- inters[1,]
inters <- inters[-1,]
inters[is.na(inters)] <-""

saveRDS(inters,'results/inters_ALZ_PSP.rds')

### Mouse DGE

FAD5X <- read.csv('archives/FAD5X_sig.csv', header = T, stringsAsFactors = F,
                  row.names = 1)
Tau <- read.csv('archives/Tau_sig.csv', header = T, stringsAsFactors = F,
                row.names = 1)

FAD_4CCX <- FAD5X[FAD5X$Month == '4' & FAD5X$Area == 'Cortex',]
FAD_12CCX <- FAD5X[FAD5X$Month == '12' & FAD5X$Area == 'Cortex',]
FAD_18CCX <- FAD5X[FAD5X$Month == '18' & FAD5X$Area == 'Cortex',]

FAD_4HIP <- FAD5X[FAD5X$Month == '4' & FAD5X$Area == 'Hippocampus',]
FAD_12HIP <- FAD5X[FAD5X$Month == '12' & FAD5X$Area == 'Hippocampus',]
FAD_18HIP <- FAD5X[FAD5X$Month == '18' & FAD5X$Area == 'Hippocampus',]

TAU_4 <- Tau[Tau$Month == '4',]
TAU_17 <- Tau[Tau$Month == '17',]

FAD_CCX <- list(FAD_4_CCX = unique(FAD_4CCX$Gene_Symbol),
            FAD_12_CCX = unique(FAD_12CCX$Gene_Symbol),
            FAD_18_CCX = unique(FAD_18CCX$Gene_Symbol))
            
FAD_HIP <- list(FAD_4_HIP = unique(FAD_4HIP$Gene_Symbol),
                 FAD_12_HIP = unique(FAD_12HIP$Gene_Symbol),
                 FAD_18_HIP = unique(FAD_18HIP$Gene_Symbol))


TAU <- list(TAU_4_HIP = unique(TAU_4$Gene_Name),
            TAU_17_HIP = unique(TAU_17$Gene_Name))

intercept_FAD_CCX <- venn::venn(FAD_CCX,
            zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)

intercept_FAD_HIP <- venn::venn(FAD_HIP,
            zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)


intercept_TAU <- venn::venn(TAU,
      zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)


### FAD-TAU
FAD_TAU <- list(FAD_4_CCX = unique(FAD_4CCX$Gene_Symbol),
                FAD_12_CCX = unique(FAD_12CCX$Gene_Symbol),
                FAD_18_CCX = unique(FAD_18CCX$Gene_Symbol),
                FAD_4_HIP = unique(FAD_4HIP$Gene_Symbol),
                FAD_12_HIP = unique(FAD_12HIP$Gene_Symbol),
                FAD_18_HIP = unique(FAD_18HIP$Gene_Symbol),
                TAU_4_HIP = unique(TAU_4$Gene_Name),
                TAU_17_HIP = unique(TAU_17$Gene_Name))

upset(fromList(FAD_TAU),sets = c('FAD_4_CCX','FAD_12_CCX','FAD_18_CCX',
                                 'FAD_4_HIP','FAD_12_HIP','FAD_18_HIP',
                                 'TAU_4_HIP','TAU_17_HIP'),
                               order.by = c('freq','degree'))

intercept_FAD_TAU <- venn::venn(FAD_TAU,
        zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)
inters_FAD_TAU <-attr(intercept_FAD_TAU,"intersections")

inters <- plyr::ldply(inters_FAD_TAU, rbind)
inters <- t(inters)
inters <- data.frame(inters,stringsAsFactors = F)
colnames(inters) <- inters[1,]
inters <- inters[-1,]
inters[is.na(inters)] <-""

saveRDS(inters,'results/inters_FAD_TAU.rds')

### Human DTU

mayo <- readRDS('results/isoform_results/mayo_ISW_sig.rds')

mayo_A <- mayo[mayo$condition_1 == 'AD' & mayo$AgeAtDeath == 'A',]
mayo_B <- mayo[mayo$condition_1 == 'AD' & mayo$AgeAtDeath == 'B',]
mayo_C <- mayo[mayo$condition_1 == 'AD' & mayo$AgeAtDeath == 'C',]

psp_A <- mayo[mayo$condition_1 == 'PSP' & mayo$AgeAtDeath == 'A',]
psp_B <- mayo[mayo$condition_1 == 'PSP' & mayo$AgeAtDeath == 'B',]

mayo_psp_DTU <- list(ALZ_A = unique(mayo_A$gene_name),
                     ALZ_B = unique(mayo_B$gene_name), 
                     ALZ_C = unique(mayo_C$gene_name),
                     PSP_A = unique(psp_A$gene_name),
                     PSP_B = unique(psp_B$gene_name)
                     )

upset(fromList(mayo_psp_DTU),sets =  c('ALZ_A','ALZ_B','ALZ_C',
                    'PSP_A','PSP_B'),order.by = c('freq','degree'))

mayo_psp_DTU <- venn::venn(mayo_psp_DTU,
      zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)
inters_mayo_psp_DTU <-attr(mayo_psp_DTU,"intersections")

inters <- plyr::ldply(inters_mayo_psp_DTU, rbind)
inters <- t(inters)
inters <- data.frame(inters,stringsAsFactors = F)
colnames(inters) <- inters[1,]
inters <- inters[-1,]
inters[is.na(inters)] <-""

saveRDS(inters,'results/inter_mayo_psp_DTU.rds')


### Mouse DTU

FAD5X_DTU <- readRDS('results/isoform_results/FAD5X_ISW_sig.rds')

Tau_DTU <- readRDS('results/isoform_results/Tau_ISW.rds')

FAD_4CCX <- FAD5X_DTU[FAD5X_DTU$Group == '4' & FAD5X_DTU$Region == 'Cortex',]
FAD_12CCX <- FAD5X_DTU[FAD5X_DTU$Group == '12' & FAD5X_DTU$Region == 'Cortex',]
FAD_18CCX <- FAD5X_DTU[FAD5X_DTU$Group == '18' & FAD5X_DTU$Region == 'Cortex',]

FAD_4HIP <- FAD5X_DTU[FAD5X_DTU$Group == '4' & FAD5X_DTU$Region == 'Hippocampus',]
FAD_12HIP <- FAD5X_DTU[FAD5X_DTU$Group == '12' & FAD5X_DTU$Region == 'Hippocampus',]
FAD_18HIP <- FAD5X_DTU[FAD5X_DTU$Group == '4' & FAD5X_DTU$Region == 'Hippocampus',]

TAU_4 <- Tau_DTU[Tau_DTU$Group == '4',]
TAU_17 <- Tau_DTU[Tau_DTU$Group  == '17',]

FAD_TAU_DTU <- list(
                TAU_4_HIP = unique(TAU_4$gene_name),
                TAU_17_HIP = unique(TAU_17$gene_name))

upset(fromList(FAD_TAU_DTU),sets = c('FAD_4_CCX','FAD_12_CCX','FAD_18_CCX',
                                 'FAD_4_HIP','FAD_12_HIP','FAD_18_HIP',
                                 'TAU_4_HIP','TAU_17_HIP'),order.by = "freq")

intercept_FAD_TAU <- venn::venn(FAD_TAU_DTU,
        zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)
inters_FAD_TAU_DTU <-attr(intercept_FAD_TAU,"intersections")

inters <- plyr::ldply(inters_FAD_TAU_DTU, rbind)
inters <- t(inters)
inters <- data.frame(inters,stringsAsFactors = F)
colnames(inters) <- inters[1,]
inters <- inters[-1,]
inters[is.na(inters)] <-""

saveRDS(inters,'results/inters_FAD_TAU_DTU.rds')
