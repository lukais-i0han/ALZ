#### Cluster Analysis
library('transcriptogramer')
library('ggplot2')
library('stringr')

### Human

ALZ_all <- readRDS('trans_results/MAYO/clusters/all_ALZ_cluster.rds')
ALZ_B <- readRDS('trans_results/MAYO/clusters/t_B_cluster_ALZ.rds')
ALZ_C <- readRDS('trans_results/MAYO/clusters/t_C_cluster_ALZ.rds')

PSP_all <- readRDS('trans_results/PSP/clusters/all_cluster_PSP.rds')
PSP_A <- readRDS('trans_results/PSP/clusters/t_A_cluster_PSP.rds')
PSP_B <- readRDS('trans_results/PSP/clusters/t_B_cluster_PSP.rds')

inters_ALZ_PSP <- readRDS('results/tables/inters_ALZ_PSP.rds')

### ALZ

ALZ_GO_ENSP <- ALZ_all@Protein2GO
ALZ_GO_ENSP$ensembl_peptide_id <- gsub('^ENSP','9606.ENSP',
                                       ALZ_GO_ENSP$ensembl_peptide_id)
colnames(ALZ_GO_ENSP) <- c('Protein','GO:ID')

## all 

inters_ALZ <- inters_ALZ_PSP[,!(str_detect(colnames(inters_ALZ_PSP),'PSP'))]
inters_ALZ <- c(inters_ALZ$ALZ_C,inters_ALZ$ALZ_B,inters_ALZ$`ALZ_B:ALZ_C`,
                inters_ALZ$ALZ_A,inters_ALZ$`ALZ_A:ALZ_C`,
                inters_ALZ$`ALZ_A:ALZ_B`,inters_ALZ$`ALZ_A:ALZ_B:ALZ_C`)

ALZ_all_DEG <- ALZ_all@DE
ALZ_all_DEG$'Age_Group' <- 'A_B_C'
ALZ_all_DEG <- ALZ_all_DEG[ALZ_all_DEG$Symbol %in% inters_ALZ,]

ALZ_GO_all <- ALZ_GO_ENSP[ALZ_GO_ENSP$Protein %in% ALZ_all_DEG$Protein,]

ALZ_all_Terms <- ALZ_all@Terms
ALZ_all_Terms$'Age_Group' <- 'A_B_C'
ALZ_all_Terms <- ALZ_all_Terms[ALZ_all_Terms$GO.ID %in% ALZ_GO_all$`GO:ID`,]


##B

inters_ALZ_B <- inters_ALZ_PSP[,!(str_detect(colnames(inters_ALZ_PSP),'PSP'))]
inters_ALZ_B <- c(inters_ALZ_B$ALZ_B,inters_ALZ_B$`ALZ_B:ALZ_C`,
                  inters_ALZ_B$`ALZ_A:ALZ_B`,inters_ALZ_B$`ALZ_A:ALZ_B:ALZ_C`)

ALZ_B_DEG <- ALZ_B@DE
ALZ_B_DEG$'Age_Group' <- 'B'
ALZ_B_DEG <- ALZ_B_DEG[ALZ_B_DEG$Symbol %in% inters_ALZ_B,]

ALZ_GO_B <- ALZ_GO_ENSP[ALZ_GO_ENSP$Protein %in% ALZ_B_DEG$Protein,]

ALZ_B_Terms <- ALZ_B@Terms
ALZ_B_Terms$'Age_Group' <- 'B'
ALZ_B_Terms <- ALZ_B_Terms[ALZ_B_Terms$GO.ID %in% ALZ_GO_B$`GO:ID`,]

##C

inters_ALZ_C <- inters_ALZ_PSP[,!(str_detect(colnames(inters_ALZ_PSP),'PSP'))]
inters_ALZ_C <- c(inters_ALZ_C$ALZ_C,inters_ALZ_C$`ALZ_B:ALZ_C`,
                 inters_ALZ_C$`ALZ_A:ALZ_B:ALZ_C`,inters_ALZ_C$`ALZ_A:ALZ_C`)

ALZ_C_DEG <- ALZ_C@DE
ALZ_C_DEG$'Age_Group' <- 'C'
ALZ_C_DEG <- ALZ_C_DEG[ALZ_C_DEG$Symbol %in% inters_ALZ_C,]

ALZ_GO_C <- ALZ_GO_ENSP[ALZ_GO_ENSP$Protein %in% ALZ_C_DEG$Protein,]

ALZ_C_Terms <- ALZ_C@Terms
ALZ_C_Terms$'Age_Group' <- 'C'
ALZ_C_Terms <- ALZ_C_Terms[ALZ_C_Terms$GO.ID %in% ALZ_GO_C$`GO:ID`,]

## 
ALZ_DEG <- rbind(ALZ_all_DEG,ALZ_B_DEG,ALZ_C_DEG)
ALZ_DEG$'Group' <- 'ALZ'

ALZ_Terms <- rbind(ALZ_all_Terms,ALZ_B_Terms,ALZ_C_Terms)
ALZ_Terms$'Group' <- 'ALZ'

### PSP

PSP_GO_ENSP <- PSP_all@Protein2GO
PSP_GO_ENSP$ensembl_peptide_id <- gsub('^ENSP','9606.ENSP',
                                       PSP_GO_ENSP$ensembl_peptide_id)
colnames(PSP_GO_ENSP) <- c('Protein','GO:ID')

## All

inters_PSP <- inters_ALZ_PSP[,!(str_detect(colnames(inters_ALZ_PSP),'ALZ'))]
inters_PSP <- c(inters_PSP$PSP_A,inters_PSP$PSP_B,inters_PSP$`PSP_A:PSP_B`)

PSP_all_DEG <- PSP_all@DE
PSP_all_DEG$'Age_Group' <- 'A_B'
PSP_all_DEG <- PSP_all_DEG[PSP_all_DEG$Symbol %in% inters_PSP,]

PSP_GO_all <- PSP_GO_ENSP[PSP_GO_ENSP$Protein %in% PSP_all_DEG$Protein,]

PSP_all_Terms <- PSP_all@Terms
PSP_all_Terms$'Age_Group' <- 'A_B'
PSP_all_Terms <- PSP_all_Terms[PSP_all_Terms$GO.ID %in% PSP_GO_all$`GO:ID`,]

## A

inters_PSP <- inters_ALZ_PSP[,!(str_detect(colnames(inters_ALZ_PSP),'ALZ'))]
inters_PSP_A <- c(inters_PSP$PSP_A,inters_PSP$`PSP_A:PSP_B`)

PSP_A_DEG <- PSP_A@DE
PSP_A_DEG$'Age_Group' <- 'A'
PSP_A_DEG <- PSP_A_DEG[PSP_A_DEG$Symbol %in% inters_PSP_A,]

PSP_GO_A <- PSP_GO_ENSP[PSP_GO_ENSP$Protein %in% PSP_A_DEG$Protein,]

PSP_A_Terms <- PSP_A@Terms
PSP_A_Terms$'Age_Group' <- 'A'
PSP_A_Terms <- PSP_A_Terms[PSP_A_Terms$GO.ID %in% PSP_GO_A$`GO:ID`,]

## B

inters_PSP <- inters_ALZ_PSP[,!(str_detect(colnames(inters_ALZ_PSP),'ALZ'))]
inters_PSP_B <- c(inters_PSP$PSP_B,inters_PSP$`PSP_A:PSP_B`)

PSP_B_DEG <- PSP_B@DE
PSP_B_DEG$'Age_Group' <- 'B'
PSP_B_DEG <- PSP_B_DEG[PSP_B_DEG$Symbol %in% inters_PSP_B,]

PSP_GO_B <- PSP_GO_ENSP[PSP_GO_ENSP$Protein %in% PSP_B_DEG$Protein,]

PSP_B_Terms <- PSP_B@Terms
PSP_B_Terms$'Age_Group' <- 'B'
PSP_B_Terms <- PSP_B_Terms[PSP_B_Terms$GO.ID %in% PSP_GO_B$`GO:ID`,]

## 

PSP_DEG <- rbind(PSP_all_DEG,PSP_A_DEG,PSP_B_DEG)
PSP_DEG$'Group' <- 'PSP'

PSP_Terms <- rbind(PSP_all_Terms,PSP_A_Terms,PSP_B_Terms)
PSP_Terms$'Group' <- 'PSP'

##

ALZ_PSP_DEG <- rbind(ALZ_DEG,PSP_DEG)
ALZ_PSP_Terms <- rbind(ALZ_Terms,PSP_Terms)

saveRDS(ALZ_PSP_DEG,'trans_results/ALZ_PSP_DEG.rds')
saveRDS(ALZ_PSP_Terms,'trans_results/ALZ_PSP_Terms.rds')

#### Mouse

### CCX

dictionary_FAD <- read.csv('results/tables/dictionary_FAD.csv',
                  header = T, stringsAsFactors = F, row.names = 1)

inters_FAD <- readRDS('results/tables/inters_FAD_TAU_DGE.rds')

inters_FAD <- inters_FAD[,!(str_detect(colnames(inters_FAD),'TAU'))]
inters_FAD <- inters_FAD[,!(str_detect(colnames(inters_FAD),'HIP'))]
inters_FAD <- c(inters_FAD$FAD_4_CCX,inters_FAD$FAD_12_CCX,
                inters_FAD$FAD_18_CCX,inters_FAD$`FAD_12_CCX:FAD_18_CCX`,
                inters_FAD$`FAD_4_CCX:FAD_18_CCX`)

inters_FAD <- inters_FAD[inters_FAD != '']

CCX_all <- readRDS('trans_results/FAD/clusters/all_cluster_CCX.rds')

CCX_GO <- CCX_all@Protein2GO

####

CCX_all <- list(CCX_all_DEG = CCX_all@DE,
                CCX_all_Terms = CCX_all@Terms)
CCX_all$CCX_all_DEG$'Month_Group' <- '4_12_18'
CCX_all$CCX_all_Terms$'Month_Group' <- '4_12_18'


CCX_12 <- readRDS('trans_results/FAD/clusters/t_cluster_CCX12.rds')
CCX_12 <- list(CCX_12_DEG = CCX_12@DE,
               CCX_12_Terms = CCX_12@Terms)
CCX_12$CCX_12_DEG$'Month_Group' <- '12'
CCX_12$CCX_12_Terms$'Month_Group' <- '12'

CCX_18 <- readRDS('trans_results/FAD/clusters/t_cluster_CCX18.rds')
CCX_18 <- list(CCX_18_DEG = CCX_18@DE,
               CCX_18_Terms = CCX_18@Terms)
CCX_18$CCX_18_DEG$'Month_Group' <- '18'
CCX_18$CCX_18_Terms$'Month_Group' <- '18'

CCX_DEG <- rbind(CCX_all$CCX_all_DEG,CCX_12$CCX_12_DEG,
                 CCX_18$CCX_18_DEG)
CCX_DEG$'Region' <- 'Cortex'

CCX_Terms <- rbind(CCX_all$CCX_all_Terms,CCX_12$CCX_12_Terms,
                   CCX_18$CCX_18_Terms)
CCX_Terms$'Region' <- 'Cortex'

####
colnames(CCX_GO) <- c('Protein','GO:ID')
CCX_GO$Protein <- gsub('^ENSMUSP','10090.ENSMUSP',CCX_GO$Protein)

dictionary_FAD <- dictionary_FAD[dictionary_FAD$external_gene_name %in% 
                                   inters_FAD,]

CCX_DEG <- CCX_DEG[CCX_DEG$Symbol %in% dictionary_FAD$ensembl_peptide_id,]

CCX_GO <- CCX_GO[CCX_GO$Protein %in% CCX_DEG$Protein, ]

CCX_Terms <- CCX_Terms[CCX_Terms$GO.ID %in% CCX_GO$`GO:ID`,]

### HIP


dictionary_FAD <- read.csv('results/tables/dictionary_FAD.csv',
                           header = T, stringsAsFactors = F, row.names = 1)

inters_FAD <- readRDS('results/tables/inters_FAD_TAU_DGE.rds')

inters_FAD <- inters_FAD[,!(str_detect(colnames(inters_FAD),'TAU'))]
inters_FAD <- inters_FAD[,!(str_detect(colnames(inters_FAD),'CCX'))]
inters_FAD <- c(inters_FAD$FAD_4_HIP,inters_FAD$FAD_12_HIP,
              inters_FAD$FAD_18_HIP, inters_FAD$`FAD_4_HIP:FAD_12_HIP`,
              inters_FAD$`FAD_12_HIP:FAD_18_HIP`,
              inters_FAD$`FAD_4_HIP:FAD_12_HIP`)

inters_FAD <- inters_FAD[inters_FAD != '']

HIP_all <- readRDS('trans_results/FAD/clusters/all_cluster_HIP.rds')

HIP_GO <- HIP_all@Protein2GO

####

HIP_all <- list(HIP_all_DEG = HIP_all@DE,
                HIP_all_Terms = HIP_all@Terms)
HIP_all$HIP_all_DEG$'Month_Group' <- '4_12_18'
HIP_all$HIP_all_Terms$'Month_Group' <- '4_12_18'

HIP_4 <- readRDS('trans_results/FAD/clusters/t_cluster_HIP4.rds')
HIP_4 <- list(HIP_4_DEG = HIP_4@DE,
              HIP_4_Terms = HIP_4@Terms)
HIP_4$HIP_4_DEG$'Month_Group' <- '4'
HIP_4$HIP_4_Terms$'Month_Group' <- '4'

HIP_12 <- readRDS('trans_results/FAD/clusters/t_cluster_HIP12.rds')
HIP_12 <- list(HIP_12_DEG = HIP_12@DE,
               HIP_12_Terms = HIP_12@Terms)
HIP_12$HIP_12_DEG$'Month_Group' <- '12'
HIP_12$HIP_12_Terms$'Month_Group' <- '12'

HIP_18 <- readRDS('trans_results/FAD/clusters/t_cluster_HIP18.rds')
HIP_18 <- list(HIP_18_DEG = HIP_18@DE,
               HIP_18_Terms = HIP_18@Terms)
HIP_18$HIP_18_DEG$'Month_Group' <- '18'
HIP_18$HIP_18_Terms$'Month_Group' <- '18'

HIP_DEG <- rbind(HIP_all$HIP_all_DEG,HIP_4$HIP_4_DEG,HIP_12$HIP_12_DEG,
                 HIP_18$HIP_18_DEG)
HIP_DEG$'Region' <- 'Hippocampus'

HIP_Terms <- rbind(HIP_all$HIP_all_Terms,HIP_4$HIP_4_Terms,HIP_12$HIP_12_Terms,
                   HIP_18$HIP_18_Terms)
HIP_Terms$'Region' <- 'Hippocampus'

####
colnames(HIP_GO) <- c('Protein','GO:ID')
HIP_GO$Protein <- gsub('^ENSMUSP','10090.ENSMUSP',HIP_GO$Protein)

dictionary_FAD <- dictionary_FAD[dictionary_FAD$external_gene_name %in% 
                                   inters_FAD,]

HIP_DEG <- HIP_DEG[HIP_DEG$Symbol %in% dictionary_FAD$ensembl_peptide_id,]

HIP_GO <- HIP_GO[HIP_GO$Protein %in% HIP_DEG$Protein, ]

HIP_Terms <- HIP_Terms[HIP_Terms$GO.ID %in% HIP_GO$`GO:ID`,]



###
CCX_HIP_DEG <- rbind(CCX_DEG,HIP_DEG)

CCX_HIP_Terms <- rbind(CCX_Terms,HIP_Terms)

saveRDS(CCX_HIP_DEG,'trans_results/CCX_HIP_DEG.rds')
saveRDS(CCX_HIP_Terms,'trans_results/CCX_HIP_Terms.rds')
