#### Transcriptgramer
library('stringr')
library('transcriptogramer')

### import vsd objects

vsd_ALZ <- readRDS('results/vsd_filtered/ALZ_filtered.rds')
vsd_PSP <- readRDS('results/vsd_filtered/PSP_filtered.rds')
vsd_FAD <- readRDS('results/vsd_filtered/FAD_filtered.rds')
vsd_TAU <- readRDS('results/vsd_filtered/TAU_filtered.rds')

### import of tables with genes of each condition

ALZ <- read.csv('archives/ALZ_all.csv', header = T, stringsAsFactors = F,
                 row.names = 1)



PSP <- read.csv('archives/PSP_all.csv',header = T, stringsAsFactors = F,
                row.names = 1)

### Mayo - ALZ

ALZ$ENSG <-str_sub(ALZ$ENSG,start = 1,end = 15)

ALZ_A <- ALZ[ALZ$Age_Group == 'A',]
ALZ_B <- ALZ[ALZ$Age_Group == 'B',]
ALZ_C <- ALZ[ALZ$Age_Group == 'C',]

mayo_metadado <- read.csv('metadado/mayo_metadado.csv',header = T, stringsAsFactors = F,
                          row.names = 1)

ALZ_metadado <- mayo_metadado[mayo_metadado$Group != 'PSP',]

vsd_ALZA <- vsd_ALZ[,str_detect(ALZ_metadado$Diag_Age,'_A')]
vsd_ALZA <- vsd_ALZA[rownames(vsd_ALZA) %in% ALZ_A$ENSG,]

vsd_ALZB <- vsd_ALZ[,str_detect(ALZ_metadado$Diag_Age,'_B')]
vsd_ALZB <- vsd_ALZB[rownames(vsd_ALZB) %in% ALZ_B$ENSG,]

vsd_ALZC <- vsd_ALZ[,str_detect(ALZ_metadado$Diag_Age,'_C')]
vsd_ALZC <- vsd_ALZC[rownames(vsd_ALZC) %in% ALZ_C$ENSG,]


vsd_ALZ <- vsd_ALZ[rownames(vsd_ALZ) %in% ALZ$ENSG,]

t_human <- transcriptogramPreprocess(association = association, ordering = Hs900,
                               radius = 80)

ALZ_dictionary <- read.csv('results/tables/ALZ_dictionary.csv',header = T, stringsAsFactors = F,
                            row.names = 1)

t_all <- transcriptogramStep1(object = t_human, expression = vsd_ALZ,
                                       dictionary = ALZ_dictionary,nCores = 4)

t_A <- transcriptogramStep1(object = t_human, expression = vsd_ALZA,
                            dictionary = ALZ_dictionary,nCores = 4)

t_B <- transcriptogramStep1(object = t_human, expression = vsd_ALZB,
                            dictionary = ALZ_dictionary,nCores = 4)

t_C <- transcriptogramStep1(object = t_human, expression = vsd_ALZC,
                            dictionary = ALZ_dictionary,nCores = 4)

### Second Step

t_all <- transcriptogramStep2(object = t_all, nCores = 4 )
saveRDS(t_all,'trans_results/t_all_ALZ.rds')

t_A <- transcriptogramStep2(object = t_A, nCores = 4 )
saveRDS(t_A,'trans_results/t_A_ALZ.rds')

t_B <- transcriptogramStep2(object = t_B, nCores = 4 )
saveRDS(t_B,'trans_results/t_B_ALZ.rds')

t_C <- transcriptogramStep2(object = t_C, nCores = 4 )
saveRDS(t_C,'trans_results/t_C_ALZ.rds')

### Third Step
ALZ_A <- ALZ_metadado[ALZ_metadado$Age_Group == 'A',]
ALZ_B <- ALZ_metadado[ALZ_metadado$Age_Group == 'B',]
ALZ_C <- ALZ_metadado[ALZ_metadado$Age_Group == 'C',]


t_all <- differentiallyExpressed(object = t_all,
          levels = ifelse(ALZ_metadado$Group == 'Control',TRUE,FALSE)
          , pValue = 0.01,adjustMethod = 'fdr',
          trend = F, title = "All Groups",DEsymbols)

t_A <- differentiallyExpressed(object = t_A,
        levels = ifelse(ALZ_A$Group == 'Control',TRUE,FALSE), pValue = 0.01,
        trend = F, title = "Group A",adjustMethod = 'fdr',DEsymbols) #### No clusters differentially expressed


t_B <- differentiallyExpressed(object = t_B,
    levels = ifelse(ALZ_B$Group == 'Control',TRUE,FALSE), pValue = 0.01,
    trend = F, title = "Group B",DEsymbols) 

t_C <- differentiallyExpressed(object = t_C,
          levels = ifelse(ALZ_C$Group == 'Control',TRUE,FALSE), pValue = 0.01,
          trend = F, title = "Group C",DEsymbols)


###Cluster enrinchement
t_all_cluster <- clusterEnrichment(object = t_all, species = 'Homo sapiens',
                  pValue = 0.005, nCores = 5)
saveRDS(t_all_cluster,'trans_results/all_cluster.rds')

t_B_cluster <- clusterEnrichment(object = t_B, species = 'Homo sapiens',
                              pValue = 0.005, nCores = 5)
saveRDS(t_B_cluster,'trans_results/t_B_cluster_mayo.rds')

t_C_cluster <- clusterEnrichment(object = t_C, species = 'Homo sapiens',
                                 pValue = 0.005, nCores = 5)
saveRDS(t_C_cluster,'trans_results/t_C_cluster_mayo.rds')

###  PSP

PSP$ENSG <-str_sub(PSP$ENSG,start = 1,end = 15)

PSP_A <- PSP[PSP$Age_Group == 'A',]
PSP_B <- PSP[PSP$Age_Group == 'B',]

PSP_metadado <- read.csv('metadado/mayo_metadado.csv',header = T, stringsAsFactors = F,
                          row.names = 1)

PSP_metadado <- PSP_metadado[PSP_metadado$Group != 'AD' &
                            PSP_metadado$Age_Group != 'C',]

vsd_PSPA <- vsd_PSP[,str_detect(PSP_metadado$Diag_Age,'_A')]
vsd_PSPA <- vsd_PSPA[rownames(vsd_PSPA) %in% PSP_A$ENSG,]

vsd_PSPB <- vsd_PSP[,str_detect(PSP_metadado$Diag_Age,'_B')]
vsd_PSPB <- vsd_PSPB[rownames(vsd_PSPB) %in% PSP_B$ENSG,]

vsd_PSP <- vsd_PSP[rownames(vsd_PSP) %in% PSP$ENSG,]

t_human <- transcriptogramPreprocess(association = association, ordering = Hs900,
                                     radius = 80)

dictionary_PSP <- read.csv('results/tables/PSP_dictionary.csv',header = T, stringsAsFactors = F,
                            row.names = 1)

t_all <- transcriptogramStep1(object = t_human, expression = vsd_PSP,
                              dictionary = dictionary_PSP,nCores = 4)

t_A <- transcriptogramStep1(object = t_human, expression = vsd_PSPA,
                            dictionary = dictionary_PSP,nCores = 4)

t_B <- transcriptogramStep1(object = t_human, expression = vsd_PSPB,
                            dictionary = dictionary_PSP,nCores = 4)

### Second Step

t_all <- transcriptogramStep2(object = t_all, nCores = 4 )
saveRDS(t_all,'trans_results/t_all_PSP.rds')

t_A <- transcriptogramStep2(object = t_A, nCores = 4 )
saveRDS(t_A,'trans_results/t_A_PSP.rds')

t_B <- transcriptogramStep2(object = t_B, nCores = 4 )
saveRDS(t_B,'trans_results/t_B_PSP.rds')


### Third Step
PSP_metadadoA <- PSP_metadado[PSP_metadado$Age_Group == 'A',]
PSP_metadadoB <- PSP_metadado[PSP_metadado$Age_Group == 'B',]


t_all <- differentiallyExpressed(object = t_all,
        levels = ifelse(PSP_metadado$Group == 'Control',TRUE,FALSE)
        , pValue = 0.01,
        trend = F, title = "All Groups",DEsymbols)


t_A <- differentiallyExpressed(object = t_A,
        levels = ifelse(PSP_metadadoA$Group == 'Control',TRUE,FALSE), pValue = 0.01,
        trend = F, title = "Group A",DEsymbols) 

t_B <- differentiallyExpressed(object = t_B,
        levels = ifelse(PSP_metadadoB$Group == 'Control',TRUE,FALSE), pValue = 0.01,
        trend = F, title = "Group B",DEsymbols) 



###Cluster enrinchement
t_all_cluster <- clusterEnrichment(object = t_all, species = 'Homo sapiens',
                                   pValue = 0.005, nCores = 5)
saveRDS(t_all_cluster,'trans_results/all_cluster_PSP.rds')

t_A_cluster <- clusterEnrichment(object = t_A, species = 'Homo sapiens',
                                 pValue = 0.005, nCores = 5)
saveRDS(t_A_cluster,'trans_results/t_A_cluster_PSP.rds')

t_B_cluster <- clusterEnrichment(object = t_B, species = 'Homo sapiens',
                                 pValue = 0.005, nCores = 5)
saveRDS(t_B_cluster,'trans_results/t_B_cluster_PSP.rds')


#### Mouse

FAD5X <- read.csv('archives/FAD5X_all.csv', header = T, stringsAsFactors = F,
                  row.names = 1)

FAD5X$ENSG <- str_sub(FAD5X$ENSG,start = 1,end = 18)

FAD_metadado <- read.csv('metadata/FAD_metadata.csv',header = T, stringsAsFactors = F,
                          row.names = 1)

dictionary_FAD <- read.csv('results/dictionary_FAD.csv',header = T, stringsAsFactors = F,
                           row.names = 1)

FAD_CCX <- FAD5X[FAD5X$Area == 'Cortex',]
FAD_HIP <- FAD5X[FAD5X$Area == 'Hippocampus',]

vsd_FAD4_CCX <- vsd_FAD[,FAD_metadado$Region == 'CCX' & FAD_metadado$Month == '4']
FAD_4_CCX <- FAD5X[FAD5X$Area == 'Cortex' & FAD5X$Month == '4',]
vsd_FAD4_CCX <- vsd_FAD4_CCX[rownames(vsd_FAD4_CCX) %in% FAD_4_CCX$ENSG,]

vsd_FAD12_CCX <- vsd_FAD[,FAD_metadado$Region == 'CCX' & FAD_metadado$Month == '12']
FAD_12_CCX <- FAD5X[FAD5X$Area == 'Cortex' & FAD5X$Month == '12',]
vsd_FAD12_CCX <- vsd_FAD12_CCX[rownames(vsd_FAD12_CCX) %in% FAD_12_CCX$ENSG,]

vsd_FAD18_CCX <- vsd_FAD[,FAD_metadado$Region == 'CCX' & FAD_metadado$Month == '18']
FAD_18_CCX <- FAD5X[FAD5X$Area == 'Cortex' & FAD5X$Month == '18',]
vsd_FAD18_CCX <- vsd_FAD18_CCX[rownames(vsd_FAD18_CCX) %in% FAD_18_CCX$ENSG,]

vsd_FAD4_HIP <- vsd_FAD[,FAD_metadado$Region == 'HIP' & FAD_metadado$Month == '4']
FAD_4_HIP <- FAD5X[FAD5X$Area == 'Hippocampus' & FAD5X$Month == '4',]
vsd_FAD4_HIP <- vsd_FAD4_HIP[rownames(vsd_FAD4_HIP) %in% FAD_4_HIP$ENSG,]

vsd_FAD12_HIP <- vsd_FAD[,FAD_metadado$Region == 'HIP' & FAD_metadado$Month == '12']
FAD_12_HIP <- FAD5X[FAD5X$Area == 'Hippocampus' & FAD5X$Month == '12',]
vsd_FAD12_HIP <- vsd_FAD12_HIP[rownames(vsd_FAD12_HIP) %in% FAD_12_HIP$ENSG,]

vsd_FAD18_HIP <- vsd_FAD[,FAD_metadado$Region == 'HIP' & FAD_metadado$Month == '18']
FAD_18_HIP <- FAD5X[FAD5X$Area == 'Hippocampus' & FAD5X$Month == '18',]
vsd_FAD18_HIP <- vsd_FAD18_HIP[rownames(vsd_FAD18_HIP) %in% FAD_18_HIP$ENSG,]

association_mouse <- read.table('association_mouse.txt',header = F, stringsAsFactors = F)
colnames(association_mouse) <- c('protein1','protein2')
ordering_Mm900 <- read.table('ordering_Mm900.txt',header = F,stringsAsFactors = F )
ordering_Mm900 <- ordering_Mm900$V1

t_mouse <- transcriptogramPreprocess(association = association_mouse, ordering = ordering_Mm900,
                        radius = 80)

dictionary_FAD <- read.csv('results/dictionary_FAD.csv',header = T, stringsAsFactors = F,
                            row.names = 1)
dictionary_FAD <- dictionary_FAD[,1:2]

t_all_CCX <- transcriptogramStep1(object = t_mouse, 
                              expression = vsd_FAD[,FAD_metadado$Region == 'CCX'],
                              dictionary = dictionary_FAD,nCores = 4)

t_4_CCX <- transcriptogramStep1(object = t_mouse, 
                            expression = vsd_FAD4_CCX,
                            dictionary = dictionary_FAD,nCores = 4)

t_12_CCX <- transcriptogramStep1(object = t_mouse, 
                                 expression = vsd_FAD12_CCX,
                                 dictionary = dictionary_FAD,nCores = 4)

t_18_CCX <- transcriptogramStep1(object = t_mouse, 
                                 expression = vsd_FAD18_CCX,
                                 dictionary = dictionary_FAD,nCores = 4)


t_all_HIP <- transcriptogramStep1(object = t_mouse, 
                                  expression = vsd_FAD[,FAD_metadado$Region == 'HIP'],
                                  dictionary = dictionary_FAD,nCores = 4)

t_4_HIP <- transcriptogramStep1(object = t_mouse, 
                                expression = vsd_FAD4_HIP,
                                dictionary = dictionary_FAD,nCores = 4)

t_12_HIP <- transcriptogramStep1(object = t_mouse, 
                                 expression = vsd_FAD12_HIP,
                                 dictionary = dictionary_FAD,nCores = 4)

t_18_HIP <- transcriptogramStep1(object = t_mouse, 
                                 expression = vsd_FAD18_HIP,
                                 dictionary = dictionary_FAD,nCores = 4)
### Second Step

##CCX
t_all_CCX <- transcriptogramStep2(object = t_all_CCX, nCores = 4 )
saveRDS(t_all_CCX,'trans_results/t_all_CCX.rds')

t_4_CCX<- transcriptogramStep2(object = t_4_CCX, nCores = 4 )
saveRDS(t_4_CCX,'trans_results/t_4_CCX.rds')

t_12_CCX<- transcriptogramStep2(object = t_12_CCX, nCores = 4 )
saveRDS(t_12_CCX,'trans_results/t_12_CCX.rds')

t_18_CCX<- transcriptogramStep2(object = t_18_CCX, nCores = 4 )
saveRDS(t_18_CCX,'trans_results/t_18_CCX.rds')

##HIP
t_all_HIP <- transcriptogramStep2(object = t_all_HIP, nCores = 4 )
saveRDS(t_all_HIP,'trans_results/t_all_HIP.rds')

t_4_HIP<- transcriptogramStep2(object = t_4_HIP, nCores = 4 )
saveRDS(t_4_HIP,'trans_results/t_4_HIP.rds')

t_12_HIP<- transcriptogramStep2(object = t_12_HIP, nCores = 4 )
saveRDS(t_12_HIP,'trans_results/t_12_HIP.rds')

t_18_HIP<- transcriptogramStep2(object = t_18_HIP, nCores = 4 )
saveRDS(t_18_HIP,'trans_results/t_18_HIP.rds')

### Third Step
FAD_CCX <- FAD_metadado[FAD_metadado$Region == 'CCX',]

FAD_CCX_4 <- FAD_CCX[FAD_CCX$Month == '4',]
FAD_CCX_12 <- FAD_CCX[FAD_CCX$Month == '12',]
FAD_CCX_18 <- FAD_CCX[FAD_CCX$Month == '18',]

FAD_HIP <- FAD_metadado[FAD_metadado$Region == 'HIP',]
FAD_HIP_4 <- FAD_HIP[FAD_HIP$Month == '4',]
FAD_HIP_12 <- FAD_HIP[FAD_HIP$Month == '12',]
FAD_HIP_18 <- FAD_HIP[FAD_HIP$Month == '18',]

t_all_CCX <- differentiallyExpressed(object = t_all_CCX,
        levels = ifelse(FAD_CCX$Group == 'Control',TRUE,FALSE), 
        pValue = 0.01,
        trend = F, title = "All Groups",DEsymbols)

t_4_CCX <- differentiallyExpressed(object = t_4_CCX,
          levels = ifelse(FAD_CCX_4$Group == 'Control',TRUE,FALSE), pValue = 0.01,
          trend = F, title = "4 months",DEsymbols) #### No clusters differentially expressed


t_12_CCX <- differentiallyExpressed(object = t_12_CCX,
          levels = ifelse(FAD_CCX_12$Group == 'Control',TRUE,FALSE), pValue = 0.01,
        trend = F, title = "12 months",DEsymbols) 


t_18_CCX <- differentiallyExpressed(object = t_18_CCX,
            levels = ifelse(FAD_CCX_18$Group == 'Control',TRUE,FALSE), pValue = 0.01,
            trend = F, title = "18 months",DEsymbols) 
###Hippocampus
t_all_HIP <- differentiallyExpressed(object = t_all_HIP,
             levels = ifelse(FAD_HIP$Group == 'Control',TRUE,FALSE), 
             pValue = 0.01,
             trend = F, title = "All Groups",DEsymbols)

t_4_HIP <- differentiallyExpressed(object = t_4_HIP,
           levels = ifelse(FAD_HIP_4$Group == 'Control',TRUE,FALSE), pValue = 0.01,
           trend = F, title = "4 months",DEsymbols)


t_12_HIP <- differentiallyExpressed(object = t_12_HIP,
            levels = ifelse(FAD_HIP_12$Group == 'Control',TRUE,FALSE), pValue = 0.01,
            trend = F, title = "12 months",DEsymbols) 


t_18_HIP <- differentiallyExpressed(object = t_18_HIP,
            levels = ifelse(FAD_HIP_18$Group == 'Control',TRUE,FALSE), pValue = 0.01,
            trend = F, title = "18 months",DEsymbols) 

###Cluster enrinchement

## CCX
t_all_cluster_CCX <- clusterEnrichment(object = t_all_CCX, species = 'Mus musculus',
                                   pValue = 0.005, nCores = 5)
saveRDS(t_all_cluster_CCX,'trans_results/all_cluster_CCX.rds')

t_12_cluster_CCX <- clusterEnrichment(object = t_12_CCX, species = 'Mus musculus',
                                 pValue = 0.005, nCores = 5)
saveRDS(t_12_cluster_CCX,'trans_results/t_cluster_CCX12.rds')

t_18_cluster_CCX <- clusterEnrichment(object = t_18_CCX, species = 'Mus musculus',
                                      pValue = 0.005, nCores = 5)
saveRDS(t_18_cluster_CCX,'trans_results/t_cluster_CCX18.rds')

##HIP

t_all_cluster_HIP <- clusterEnrichment(object = t_all_HIP, species = 'Mus musculus',
                 pValue = 0.005, nCores = 5)
saveRDS(t_all_cluster_HIP,'trans_results/all_cluster_HIP.rds')

t_4_cluster_HIP <- clusterEnrichment(object = t_4_HIP, species = 'Mus musculus',
                  pValue = 0.005, nCores = 5)
saveRDS(t_4_cluster_HIP,'trans_results/t_cluster_HIP4.rds')

t_12_cluster_HIP <- clusterEnrichment(object = t_12_HIP, species = 'Mus musculus',
                  pValue = 0.005, nCores = 5)
saveRDS(t_12_cluster_HIP,'trans_results/t_cluster_HIP12.rds')

t_18_cluster_HIP <- clusterEnrichment(object = t_18_HIP, species = 'Mus musculus',
                    pValue = 0.005, nCores = 5)
saveRDS(t_18_cluster_HIP,'trans_results/t_cluster_HIP18.rds')


### Tau

association_mouse <- read.table('association_mouse.txt',header = F, stringsAsFactors = F)
colnames(association_mouse) <- c('protein1','protein2')
ordering_Mm900 <- read.table('ordering_Mm900.txt',header = F,stringsAsFactors = F )
ordering_Mm900 <- ordering_Mm900$V1

TAU_metadado <- read.csv('metadata/Tau_metadado.csv',header = T, stringsAsFactors = F,
                         row.names = 1)

t_mouse <- transcriptogramPreprocess(association = association_mouse, ordering = ordering_Mm900,
                                     radius = 80)

dictionary_TAU <- read.csv('results/dictionary_TAU.csv',header = T, stringsAsFactors = F,
                           row.names = 1)
TAU_4 <- TAU_metadado[TAU_metadado$Month == '4_months',]
TAU_17 <- TAU_metadado[TAU_metadado$Month == '17_months',]

## First Step

Tau <- read.csv('archives/Tau_all.csv',header = T, stringsAsFactors = F,
                row.names = 1)

Tau_4 <- Tau[Tau$Month == '4',]
Tau_17 <- Tau[Tau$Month == '17',]

dictionary_TAU_all <- dictionary_TAU[dictionary_TAU$external_gene_name %in% Tau$Gene_Name,]
dictionary_TAU_all <- dictionary_TAU_all[,1:2]

dictionary_TAU_4 <- dictionary_TAU[dictionary_TAU$external_gene_name %in% Tau_4$Gene_Name,]
dictionary_TAU_4 <- dictionary_TAU_4[,1:2]

dictionary_TAU_17 <- dictionary_TAU[dictionary_TAU$external_gene_name %in% Tau_17$Gene_Name,]
dictionary_TAU_17 <- dictionary_TAU_17[,1:2]

t_all <- transcriptogramStep1(object = t_mouse, 
                                  expression = vsd_TAU,
                                  dictionary = dictionary_TAU_all,nCores = 4)

t_4 <- transcriptogramStep1(object = t_mouse, 
                                expression = vsd_TAU[,TAU_4$id],
                                dictionary = dictionary_TAU_4,nCores = 4)

t_17 <- transcriptogramStep1(object = t_mouse, 
                            expression = vsd_TAU[,TAU_17$id],
                            dictionary = dictionary_TAU_17,nCores = 4)

## Second Step

t_all <- transcriptogramStep2(object = t_all, nCores = 4 )
saveRDS(t_all,'trans_results/t_all.rds')

t_4 <- transcriptogramStep2(object = t_4, nCores = 4 )
saveRDS(t_4,'trans_results/t_4.rds')

t_17 <- transcriptogramStep2(object = t_17, nCores = 4 )
saveRDS(t_17,'trans_results/t_17.rds')

## Third Step

t_all <- differentiallyExpressed(object = t_all,
        levels = ifelse(TAU_metadado$mutation == 'Control',TRUE,FALSE), 
        pValue = 0.05,
        trend = F, title = "All Groups",DEsymbols) ### #### No clusters differentially expressed

t_4 <- differentiallyExpressed(object = t_4,
      levels = ifelse(TAU_4$mutation == 'Control',TRUE,FALSE), pValue = 0.005,
      trend = F, title = "4 months",DEsymbols) #### No clusters differentially expressed


t_17 <- differentiallyExpressed(object = t_17,
        levels = ifelse(TAU_17$mutation == 'Control',TRUE,FALSE), pValue = 0.005,
        trend = F, title = "17 months",DEsymbols) 

###Cluster enrinchement

t_17 <- clusterEnrichment(object = t_17, species = 'Mus musculus',
                      pValue = 0.005, nCores = 5)
saveRDS(t_17,'trans_results/t_17_cluster.rds')
