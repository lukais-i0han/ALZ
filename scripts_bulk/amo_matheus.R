#### Integrated analysis of Cemitool and DEGs/DTUs

library('CEMiTool')
library('readr')
library('dplyr')

## DEGs


ALZ_DEG <- read.csv('archives/ALZ_sig.csv',header = T, stringsAsFactors = F,row.names = 1)
PSP_DEG <- read.csv('archives/PSP_sig.csv',header = T, stringsAsFactors = F, row.names = 1)
FAD5X_DEG <- read.csv('archives/FAD5X_sig.csv',header = T, stringsAsFactors = F, row.names = 1)
TAU_DEG <- read.csv('archives/Tau_sig.csv',header = T, stringsAsFactors = F, row.names = 1)

## DTUs

MAYO_DTU <- readRDS('results/isoform_results/mayo_ISW_sig.rds')
FAD5X_DTU <- readRDS('results/isoform_results/FAD5X_ISW_sig.rds')
TAU_DTU <- readRDS('results/isoform_results/Tau_ISW_sig.rds')


## Organizing  human DEGs and DTU

ALZ_DEG <- table(ALZ_DEG$Gene_Symbol,ALZ_DEG$Age_Group)
ALZ_DEG <- as.data.frame.matrix(ALZ_DEG,stringsAsFactors = F)
ALZ_DEG$'genes' <- rownames(ALZ_DEG)
ALZ_DEG$'Group' <- 'ALZ'

ALZ_DEG <- ALZ_DEG %>% mutate(Age_Group = case_when(
  A == 1 & B == 0 & C == 0 ~ "[A]",
  B == 1 & A == 0 & C == 0 ~ "[B]",
  C == 1 & A == 0 & B == 0 ~ "[C]",
  A == 1 & B == 1 & C == 0 ~ "[A-B]",
  A == 1 & C == 1 & B == 0 ~ "[A-C]",
  B == 1 & C == 1 & A == 0 ~ "[B-C]",
  TRUE ~ "[A-B-C]"
))

PSP_DEG <- table(PSP_DEG$Gene_Symbol,PSP_DEG$Age_Group)
PSP_DEG <- as.data.frame.matrix(PSP_DEG,stringsAsFactors = F)
PSP_DEG$'genes' <- rownames(PSP_DEG)

PSP_DEG <- PSP_DEG %>% mutate(Age_Group = case_when(
  A == 1 & B == 0  ~ "[A]",
  B == 1 & A == 0 ~ "[B]",
  TRUE ~ '[A-B]'
))
PSP_DEG$'Group' <- 'PSP'

Human_DEG <- data.frame(genes = c(ALZ_DEG$genes,PSP_DEG$genes),
                        Age_Group = c(ALZ_DEG$Age_Group,PSP_DEG$Age_Group),
                        Group = c(ALZ_DEG$Group,PSP_DEG$Group),stringsAsFactors = F)

Human_DEG$'Class' <- 'DEG'

ALZ_DTU <- MAYO_DTU[MAYO_DTU$condition_1 == 'AD',]
ALZ_DTU <- ALZ_DTU %>% distinct(gene_name,Age_Group,.keep_all = T)
ALZ_DTU <- table(ALZ_DTU$gene_name,ALZ_DTU$Age_Group)
ALZ_DTU <- as.data.frame.matrix(ALZ_DTU,stringsAsFactors = F)
ALZ_DTU$'genes' <- rownames(ALZ_DTU)

ALZ_DTU <- ALZ_DTU %>% mutate(Age_Group = case_when(
  A == 1 & B == 0 & C == 0 ~ "[A]",
  B == 1 & A == 0 & C == 0 ~ "[B]",
  C == 1 & A == 0 & B == 0 ~ "[C]",
  A == 1 & B == 1 & C == 0 ~ "[A-B]",
  A == 1 & C == 1 & B == 0 ~ "[A-C]",
  B == 1 & C == 1 & A == 0 ~ "[B-C]",
  TRUE ~ "[A-B-C]"
))

ALZ_DTU$'Group' <- 'ALZ'

PSP_DTU <- MAYO_DTU[MAYO_DTU$condition_1 != 'AD',]
PSP_DTU <- PSP_DTU %>% distinct(gene_name,Age_Group)
PSP_DTU <- table(PSP_DTU$gene_name,PSP_DTU$Age_Group)
PSP_DTU <- as.data.frame.matrix(PSP_DTU,stringsAsFactors = F)
PSP_DTU$'genes' <- rownames(PSP_DTU)

PSP_DTU <- PSP_DTU %>% mutate(Age_Group = case_when(
  A == 1 & B == 0  ~ "[A]",
  B == 1 & A == 0 ~ "[B]",
  TRUE ~ '[A-B]'
))
PSP_DTU$'Group' <- 'PSP'

Human_DTU <- data.frame(genes = c(ALZ_DTU$genes,PSP_DTU$genes),
                        Age_Group = c(ALZ_DTU$Age_Group,PSP_DTU$Age_Group),
                        Group = c(ALZ_DTU$Group,PSP_DTU$Group),stringsAsFactors = F)
Human_DTU$'Class' <- 'DTU'

Human_DEG_DTU <- rbind(Human_DEG,Human_DTU)



## Organizing mouse DEGs and DTU

FAD5X_DEG_CCX <- FAD5X_DEG[FAD5X_DEG$Area == 'Cortex',]
FAD5X_DEG_CCX <- table(FAD5X_DEG_CCX$Gene_Symbol,FAD5X_DEG_CCX$Month)
FAD5X_DEG_CCX <- as.data.frame.matrix(FAD5X_DEG_CCX,stringsAsFactors = F)
FAD5X_DEG_CCX$'genes' <- rownames(FAD5X_DEG_CCX)
FAD5X_DEG_CCX$'Area' <- 'Cortex'


FAD5X_DEG_HIP <- FAD5X_DEG[FAD5X_DEG$Area == 'Hippocampus',]
FAD5X_DEG_HIP <- table(FAD5X_DEG_HIP$Gene_Symbol,FAD5X_DEG_HIP$Month)
FAD5X_DEG_HIP <- as.data.frame.matrix(FAD5X_DEG_HIP,stringsAsFactors = F)
FAD5X_DEG_HIP$'genes' <- rownames(FAD5X_DEG_HIP)
FAD5X_DEG_HIP$'Area' <- 'Hippocampus'


FAD5X_DEG <- rbind(FAD5X_DEG_CCX,FAD5X_DEG_HIP)

FAD5X_DEG <- FAD5X_DEG %>% mutate(Age_Group = case_when(
  .$`4` == 1 & .$`12` == 0 & .$`18` == 0 ~ "[4]",
  .$`12` == 1 & .$`4` == 0 & .$`18` == 0 ~ "[12]",
  .$`18` == 1 & .$`4` == 0 & .$`12` == 0 ~ "[18]",
  .$`4` == 1 & .$`12` == 1 & .$`18` == 0 ~ "[4-12]",
  .$`4` == 1 & .$`18` == 1 & .$`12` == 0 ~ "[4-18]",
  .$`12` == 1 & .$`18` == 1 & .$`4` == 0 ~ "[12-18]",
  TRUE ~ "[4-12-18]"
))

FAD5X_DEG$'Class' <- 'DEG'


FAD5X_DTU_CCX <- FAD5X_DTU[FAD5X_DTU$Region == 'Cortex',]
FAD5X_DTU_CCX <- FAD5X_DTU_CCX %>% distinct(gene_name,Age_Group)
FAD5X_DTU_CCX <- table(FAD5X_DTU_CCX$gene_name,FAD5X_DTU_CCX$Age_Group)
FAD5X_DTU_CCX <- as.data.frame.matrix(FAD5X_DTU_CCX,stringsAsFactors = F)
FAD5X_DTU_CCX$'genes' <- rownames(FAD5X_DTU_CCX)
FAD5X_DTU_CCX$'Area' <- 'Cortex'


FAD5X_DTU_HIP <- FAD5X_DTU[FAD5X_DTU$Region == 'Hippocampus',]
FAD5X_DTU_HIP <- FAD5X_DTU_HIP %>% distinct(gene_name,Age_Group)
FAD5X_DTU_HIP <- table(FAD5X_DTU_HIP$gene_name,FAD5X_DTU_HIP$Age_Group)
FAD5X_DTU_HIP <- as.data.frame.matrix(FAD5X_DTU_HIP,stringsAsFactors = F)
FAD5X_DTU_HIP$'genes' <- rownames(FAD5X_DTU_HIP)
FAD5X_DTU_HIP$'Area' <- 'Hippocampus'


FAD5X_DTU <- rbind(FAD5X_DTU_CCX,FAD5X_DTU_HIP)

FAD5X_DTU <- FAD5X_DTU %>% mutate(Age_Group = case_when(
  .$`4` == 1 & .$`12` == 0 & .$`18` == 0 ~ "[4]",
  .$`12` == 1 & .$`4` == 0 & .$`18` == 0 ~ "[12]",
  .$`18` == 1 & .$`4` == 0 & .$`12` == 0 ~ "[18]",
  .$`4` == 1 & .$`12` == 1 & .$`18` == 0 ~ "[4-12]",
  .$`4` == 1 & .$`18` == 1 & .$`12` == 0 ~ "[4-18]",
  .$`12` == 1 & .$`18` == 1 & .$`4` == 0 ~ "[12-18]",
  TRUE ~ "[4-12-18]"
))

FAD5X_DTU$'Class' <- 'DTU'

FAD5X_DEG_DTU <- rbind(FAD5X_DEG[,4:7],FAD5X_DTU[,4:7])


##TAU
TAU_DEG <- table(TAU_DEG$Gene_Name,TAU_DEG$Month)
TAU_DEG <- as.data.frame.matrix(TAU_DEG)
TAU_DEG$'genes' <- rownames(TAU_DEG)

TAU_DEG <- TAU_DEG %>% mutate(Age_Group = case_when(
  .$`4` == 1 & .$`17` == 0 ~ '[4]',
  .$`4` == 0 & .$`17` == 1 ~ '[17]',
  TRUE ~ '[4-17]'
))

TAU_DEG$'Area' <- 'Hippocampus'
TAU_DEG$'Class' <- 'DEG'

TAU_DTU <- TAU_DTU %>% distinct(gene_name,Age_Group,.keep_all = T) 
TAU_DTU <- table(TAU_DTU$gene_name,TAU_DTU$Age_Group)
TAU_DTU <- as.data.frame.matrix(TAU_DTU)
TAU_DTU$'genes' <- rownames(TAU_DTU)

TAU_DTU <- TAU_DTU %>% mutate(Age_Group = case_when(
  .$`4` == 1 & .$`17` == 0 ~ '[4]',
  .$`4` == 0 & .$`17` == 1 ~ '[17]',
  TRUE ~ '[4-17]'
))

TAU_DTU$'Area' <- 'Hippocampus'
TAU_DTU$'Class' <- 'DTU'

TAU_DEG_DTU <- rbind(TAU_DEG[,3:6],TAU_DTU[,3:6])


saveRDS(Human_DEG_DTU,'CEMtool/results/human/H_DEG_DTU.rds')
saveRDS(FAD5X_DEG_DTU,'CEMtool/results/mouse_fad/F_DEG_DTU.rds')
saveRDS(TAU_DEG_DTU,'CEMtool/results/mouse_tau/T_DEG_DTU.rds')
### Import of CEMtool results

## Human

human_cem <- list(modules= read_tsv('CEMtool/results/Tables/Human/module.tsv'),
                  ora = read_tsv('CEMtool/results/Tables/Human/ora.tsv'),
                  ppi = read_tsv('CEMtool/results/Tables/Human/interactions.tsv'))

## Mouse

FAD_cem <- list(modules= read_tsv('CEMtool/results/Tables/FAD//module.tsv'),
                ora = read_tsv('CEMtool/results/Tables/FAD/ora.tsv'),
                ppi = read_tsv('CEMtool/results/Tables/FAD/interactions.tsv'))

TAU_cem <- list(modules= read_tsv('CEMtool/results/Tables/TAU/module.tsv'),
                ora = read_tsv('CEMtool/results/Tables/TAU/ora.tsv'),
                ppi = read_tsv('CEMtool/results/Tables/TAU/interactions.tsv'))

