#### Construct of transcript abundance for stageR
library('dplyr')
library('tximport')


##Metadado

mayo_metadado <- read.csv("all_data/Mayo_TCX/refs/MayoRNAseq_RNAseq_TCX_covariates.csv",
                          stringsAsFactors = F,header = T)
mayo_metadado <- mayo_metadado %>% dplyr::select(ID,Diagnosis,AgeAtDeath) %>% dplyr::filter(Diagnosis %in% 
                       c('Control','AD','PSP')) %>% arrange(ID)
rownames(mayo_metadado) <- mayo_metadado$ID

colnames(mayo_metadado) <- c('ID','condition','AgeAtDeath')

mayo_metadado <- mayo_metadado %>%  
  mutate(AgeAtDeath = case_when(
    AgeAtDeath == '90_or_above' ~ '90',
    AgeAtDeath != '90_or_above' ~ mayo_metadado$AgeAtDeath
  )) 

mayo_metadado$AgeAtDeath <- sapply(mayo_metadado$AgeAtDeath,as.integer)

mayo_metadado <- mayo_metadado %>%  
  mutate(AgeAtDeath = case_when(
    AgeAtDeath >= 70 & AgeAtDeath <= 80 ~ 'A',
    AgeAtDeath >= 81 & AgeAtDeath <= 89 ~ 'B',
    AgeAtDeath == 90 ~ 'C',
    TRUE ~ 'Z'
  )) 

mayo_metadado$'Diag_Age' <- paste(mayo_metadado$condition,mayo_metadado$AgeAtDeath,sep='_')


mayo_metadado <- mayo_metadado[mayo_metadado$AgeAtDeath != 'Z',]

### Import Kallisto files

### ALZ

mayo_metadado_ALZ<-mayo_metadado[mayo_metadado$condition != 'PSP', ]

files_ALZ = paste0(list.files("all_data/Mayo_TCX/kallisto_TCX", full.names=T), "/abundance.tsv")
files_ALZ = grep(paste0(paste0("/",mayo_metadado_ALZ$ID), collapse="|"), files_ALZ, value=T)

kallistoQuant_ALZ <- tximport(files_ALZ, type = 'kallisto',txOut = T)

counts_ALZ <- kallistoQuant_ALZ$counts
counts_ALZ <- counts_ALZ[rowSums(counts_ALZ)>=10,]

### PSP

mayo_metadado_PSP<-mayo_metadado[mayo_metadado$condition != 'AD', ]

files_PSP = paste0(list.files("all_data/Mayo_TCX/kallisto_TCX", full.names=T), "/abundance.tsv")
files_PSP = grep(paste0(paste0("/",mayo_metadado_PSP$ID), collapse="|"), files_PSP, value=T)

kallistoQuant_PSP <- tximport(files_PSP, type = 'kallisto',txOut = T)

counts_PSP <- kallistoQuant_PSP$counts
counts_PSP <- counts_PSP[rowSums(counts_PSP)>=10,]


###

write.csv(mayo_metadado_ALZ,'stage_archives/metadado_ALZ.csv')
write.csv(mayo_metadado_PSP,'stage_archives/metadado_PSP.csv')

write.csv(counts_ALZ,'stage_archives/counts_ALZ.csv')
write.csv(counts_PSP,'stage_archives/counts_PSP.csv')



#### FAD5X

FAD5X_metadado <- read.csv('metadado/FAD_metadado.csv',header = T, 
                           stringsAsFactors = F, row.names = 1)

files_FAD = paste0(list.files("all_data/5XFAD/kallisto", full.names=T), "/abundance.tsv")
files_FAD = grep(paste0(paste0("/",FAD5X_metadado$individualID), collapse="|"), files_FAD, value=T)

kallistoQuant_FAD <- tximport(files_FAD, type = 'kallisto',txOut = T)

counts_FAD <- kallistoQuant_FAD$counts

colnames(counts_FAD) <- FAD5X_metadado$SpecimenID

counts_FAD <- counts_FAD[rowSums(counts_FAD)>=10,]

write.csv(counts_FAD,'stage_archives/counts_FAD.csv')



### Tau 

TAU_metadado <- read.csv('metadado/Tau_metadado.csv',header = T, 
                           stringsAsFactors = F, row.names = 1)

files_TAU = paste0(list.files("all_data/Tau35/deseq2/kallisto", full.names=T), "/abundance.tsv")


kallistoQuant_TAU <- tximport(files_TAU, type = 'kallisto',txOut = T)

counts_TAU <- kallistoQuant_TAU$counts

colnames(counts_TAU) <- TAU_metadado$BamFileName

counts_TAU <- counts_TAU[rowSums(counts_TAU)>=10,]

write.csv(counts_TAU,'stage_archives/refs/Mouse/counts_TAU.csv')
