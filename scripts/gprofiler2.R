#### gprofiler
library('gprofiler2')
library('dplyr')

ALZ_sig <- read.csv('archives/ALZ_sig.csv',header = T, stringsAsFactors = F, 
                    row.names = 1)
PSP_sig <- read.csv('archives/PSP_sig.csv',header = T, stringsAsFactors = F,
                    row.names = 1)

Mayo_DTU <- readRDS('results/isoform_results/mayo_ISW_sig.rds')

ALZ_DTU <- Mayo_DTU[Mayo_DTU$condition_1 == 'AD',]
PSP_DTU <- Mayo_DTU[Mayo_DTU$condition_1 == 'Control',]

### DEG

DEG_human <- gost(query = list(
  'ALZ_A' = ALZ_sig[ALZ_sig$Age_Group == 'A',9],
  'ALZ_B' = ALZ_sig[ALZ_sig$Age_Group == 'B',9],
  'ALZ_C' = ALZ_sig[ALZ_sig$Age_Group == 'C',9],
  'PSP_A' = PSP_sig[PSP_sig$Age_Group == 'A',9],
  'PSP_B' = PSP_sig[PSP_sig$Age_Group == 'B',9]),
  organism = 'hsapiens',
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')

DEG_human <- DEG_human$result
DEG_human$'Query' <- 'DEG'

DTU_human <- gost(query = list(
  'ALZ_A' = ALZ_DTU[ALZ_DTU$Age_Group == 'A',5],
  'ALZ_B' = ALZ_DTU[ALZ_DTU$Age_Group == 'B',5],
  'ALZ_C' = ALZ_DTU[ALZ_DTU$Age_Group == 'C',5],
  'PSP_A' = PSP_DTU[PSP_DTU$Age_Group == 'A',5],
  'PSP_B' = PSP_DTU[PSP_DTU$Age_Group == 'B',5]),
  organism = 'hsapiens',
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')

DTU_human <- DTU_human$result
DTU_human$'Query' <- 'DTU'

Human_ontology <- rbind(DEG_human,DTU_human)

saveRDS(Human_ontology,'results/Human_gp2.rds')

### Mouse - FAD

FAD5X <- read.csv('archives/FAD5X_sig.csv',header = T,stringsAsFactors = F,
                  row.names = 1)

FAD5X_ISW <- readRDS('results/isoform_results/FAD5X_ISW_sig.rds')
FAD5X_ISW_CCX <- FAD5X_ISW[FAD5X_ISW$Region == 'Cortex',]
FAD5X_ISW_HIP <- FAD5X_ISW[FAD5X_ISW$Region == 'Hippocampus',]

FAD5X_CCX <- FAD5X[FAD5X$Area == 'Cortex',]
FAD5X_HIP <- FAD5X[FAD5X$Area == 'Hippocampus',]


### DEG
DEG_FAD <- gost(query = list(
  'CCX-4' = FAD5X_CCX[FAD5X_CCX$Month == '4',4],
  'CCX-12' = FAD5X_CCX[FAD5X_CCX$Month == '12',4],
  'CCX-18' = FAD5X_CCX[FAD5X_CCX$Month == '18',4],
  'HIP-4' = FAD5X_HIP[FAD5X_HIP$Month == '4',4],
  'HIP-12' = FAD5X_HIP[FAD5X_HIP$Month == '12',4],
  'HIP-18' = FAD5X_HIP[FAD5X_HIP$Month == '18',4]),
  organism = 'mmusculus',
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')

DEG_FAD <- DEG_FAD$result
DEG_FAD$'Query' <- 'DEG'

### DTU

DTU_FAD <- gost(query = list(
  'CCX-4' = FAD5X_ISW_CCX[FAD5X_ISW_CCX$Age_Group == '4',5],
  'CCX-12' = FAD5X_ISW_CCX[FAD5X_ISW_CCX$Age_Group == '12',5],
  'CCX-18' = FAD5X_ISW_CCX[FAD5X_ISW_CCX$Age_Group == '18',5],
  'HIP-4' = FAD5X_ISW_HIP[FAD5X_ISW_HIP$Age_Group == '4',5],
  'HIP-12' = FAD5X_ISW_HIP[FAD5X_ISW_HIP$Age_Group == '12',5],
  'HIP-18' = FAD5X_ISW_HIP[FAD5X_ISW_HIP$Age_Group == '18',5]),
  organism = 'mmusculus',
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')

DTU_FAD <- DTU_FAD$result
DTU_FAD$'Query' <- 'DTU'

FAD_ontology <- rbind(DEG_FAD,DTU_FAD)

#### Tau

TAU_DEG <- read.csv('archives/Tau_sig.csv',header = T, stringsAsFactors = F,
                    row.names = 1)
TAU_DTU <- readRDS('results/isoform_results/Tau_ISW_sig.rds')

## DEG
DEG_TAU <- gost(query = list(
  'HIP-4' = TAU_DEG[TAU_DEG$Month == '4',6],
  'HIP-17' = TAU_DEG[TAU_DEG$Month == '17',6]),
  organism = 'mmusculus',
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')

DEG_TAU <- DEG_TAU$result
DEG_TAU$'Query' <- 'DEG'

## DTU

DTU_TAU <- gost(query = list(
  'HIP-4' = TAU_DTU[TAU_DTU$Age_Group == '4',5],
  'HIP-17' = TAU_DTU[TAU_DTU$Age_Group == '17',5]),
  organism = 'mmusculus',
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')

DTU_TAU <- DTU_TAU$result
DTU_TAU$'Query' <- 'DTU'

Tau_ontology <- rbind(DEG_TAU,DTU_TAU)

Mouse_ontology <- list(FAD5X = FAD_ontology,
                       Taud35 = Tau_ontology)

saveRDS(Mouse_ontology,'results/Mouse_gp2.rds')
