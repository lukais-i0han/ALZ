#### Import of deseq rds data
library('dplyr')
library('stringr')
library('transcriptogramer')


#### Mayo_metadata
mayo_metadado <- read.csv('metadado/mayo_metadado.csv', 
                          header = T, 
                          stringsAsFactors = F, 
                          row.names = 1)

FAD_metadado <- read.csv('metadata/FAD_metadata.csv',
                         header = T,
                         stringsAsFactors = F,
                         row.names = 1)

TAU_metadado <- read.csv('metadata/Tau_metadado.csv',
                         header = T,
                         stringsAsFactors = F,
                         row.names = 1)

vsd_ALZ <- readRDS('results/vsd_ALZ.rds')
vsd_PSP <- readRDS('results/vsd_PSP.rds')
vsd_FAD <- readRDS('archives/vsd_FAD.rds')
vsd_TAU <- readRDS('archives/vsd_TAU.rds')

library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mart_mouse <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
### Mayo

ALZ_names <- str_sub(rownames(vsd_ALZ),start = 1,end = 15)
rownames(vsd_ALZ) <- ALZ_names


results_ALZ <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id",'entrezgene_id',
                                'hgnc_symbol'),
                 filters = "ensembl_gene_id", values = rownames(vsd_ALZ),
                 mart = mart)


results_ALZ$ensembl_peptide_id <- gsub("^ENSP",replacement = '9606.ENSP',
                                   results_ALZ$ensembl_peptide_id)

results_ALZ <- results_ALZ[results_ALZ$ensembl_peptide_id != '',]
results_ALZ <- results_ALZ[results_ALZ$ensembl_peptide_id %in% Hs900,]

vsd_ALZ  <- vsd_ALZ [rownames(vsd_ALZ ) %in% results_ALZ$ensembl_gene_id,]

ALZ_dictionary <- results_ALZ[,c(2,1)]

write.csv(ALZ_dictionary,'results/ALZ_dictionary.csv')
saveRDS(vsd_ALZ,'results/ALZ_filtered.rds')


### PSP

# query biomart

PSP_names <- str_sub(rownames(vsd_PSP),start = 1,end = 15)
rownames(vsd_PSP) <- PSP_names

results_PSP <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id",'entrezgene_id',
                                'hgnc_symbol'),
                 filters = "ensembl_gene_id", values = rownames(vsd_PSP),
                 mart = mart)


results_PSP$ensembl_peptide_id <- gsub("^ENSP",replacement = '9606.ENSP',
                                   results_PSP$ensembl_peptide_id)

results_PSP <- results_PSP[results_PSP$ensembl_peptide_id != '',]
results_PSP <- results_PSP[results_PSP$ensembl_peptide_id %in% Hs900,]
dictionary_PSP <- results_PSP[,c(2,1)]

vsd_PSP <- vsd_PSP[rownames(vsd_PSP) %in% results_PSP$ensembl_gene_id,]

write.csv(dictionary_PSP,'results/PSP_dictionary.csv')
saveRDS(vsd_PSP,'results/PSP_filtered.rds')


### FAD
association_mouse <- read.table('association_mouse.txt',header = F, stringsAsFactors = F)
colnames(association_mouse) <- c('protein1','protein2')
ordering_Mm900 <- read.table('ordering_Mm900.txt',header = F,stringsAsFactors = F )
ordering_Mm900 <- ordering_Mm900$V1

FAD_names <- str_sub(rownames(vsd_FAD),start = 1,end = 18)
rownames(vsd_FAD) <- FAD_names

results_FAD <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id",'entrezgene_id',
                                     'external_gene_name'),
                      filters = "ensembl_gene_id", values = rownames(vsd_FAD),
                      mart = mart_mouse)


results_FAD$ensembl_peptide_id <- gsub("^ENSMUSP",replacement = '10090.ENSMUSP',
                                        results_FAD$ensembl_peptide_id)

results_FAD <- results_FAD[results_FAD$ensembl_peptide_id != '',]
results_FAD <- results_FAD[results_FAD$ensembl_peptide_id %in% ordering_Mm900,]
dictionary_FAD <- results_FAD[,c(2,1,4)]

vsd_FAD <- vsd_FAD[rownames(vsd_FAD) %in% results_FAD$ensembl_gene_id,]

write.csv(dictionary_FAD,'results/dictionary_FAD.csv')
saveRDS(vsd_FAD,'results/FAD_filtered.rds')


### Tau

TAU_names <- str_sub(rownames(vsd_TAU),start = 1,end = 18)
rownames(vsd_TAU) <- TAU_names

results_TAU <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id",'entrezgene_id',
                                    'external_gene_name'),
                     filters = "ensembl_gene_id", values = rownames(vsd_TAU),
                     mart = mart_mouse)


results_TAU$ensembl_peptide_id <- gsub("^ENSMUSP",replacement = '10090.ENSMUSP',
                                       results_TAU$ensembl_peptide_id)

results_TAU <- results_TAU[results_TAU$ensembl_peptide_id != '',]
results_TAU <- results_TAU[results_TAU$ensembl_peptide_id %in% ordering_Mm900,]
dictionary_TAU <- results_TAU[,c(2,1,4)]

vsd_TAU <- vsd_TAU[rownames(vsd_TAU) %in% results_TAU$ensembl_gene_id,]

write.csv(dictionary_TAU,'results/dictionary_TAU.csv')
saveRDS(vsd_TAU,'results/TAU_filtered.rds')
