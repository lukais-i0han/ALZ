#### Import of deseq rds data
library('dplyr')
library('DESeq2')
library('stringr')
library('transcriptogramer')


#### Mayo_metadata
mayo_metadado <- read.csv('deseq/MAYO/mayo_metadado.csv', 
                  header = T, 
                  stringsAsFactors = F, 
                  row.names = 1)


mayo_names <- str_sub(rownames(vsd_mayo),start = 1,end = 15)
rownames(vsd_mayo) <- mayo_names

library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# query biomart
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id",'entrezgene_id',
                                'hgnc_symbol'),
                 filters = "hgnc_symbol", values = rownames(vsd_mayo),
                 mart = mart)


results$ensembl_peptide_id <- gsub("^ENSP",replacement = '9606.ENSP',results$ensembl_peptide_id)

results <- results[results$ensembl_peptide_id != '',]
results <- results[results$ensembl_peptide_id %in% Hs900,]

dictionary_mayo <- results[,c(2,4)]
write.csv(dictionary_mayo,'deseq/MAYO/dictionary_mayo.csv')

vsd_mayo <- vsd_mayo[rownames(vsd_mayo) %in% results$hgnc_symbol,]
saveRDS(vsd_mayo,'deseq/MAYO/vsd_mayo.rds')

##### Transcriptgrammer
vsd_mayo <- readRDS('deseq/MAYO/vsd_mayo.rds')

t <- transcriptogramPreprocess(association = association, ordering = Hs900,
                               radius = 80)

vsd_mayoA <- str_detect(mayo_metadado$Diag_Age,'_A')
vsd_mayoB <- str_detect(mayo_metadado$Diag_Age,'_B')
vsd_mayoC <- str_detect(mayo_metadado$Diag_Age,'_C')

### First Step
t_all <- transcriptogramStep1(object = t, expression = vsd_mayo,
                          dictionary = dictionary_mayo,nCores = 4)

t_A <- transcriptogramStep1(object = t, expression = vsd_mayo[,vsd_mayoA],
                              dictionary = dictionary_mayo,nCores = 4)

t_B <- transcriptogramStep1(object = t, expression = vsd_mayo[,vsd_mayoB],
                            dictionary = dictionary_mayo,nCores = 4)

t_C <- transcriptogramStep1(object = t, expression = vsd_mayo[,vsd_mayoC],
                            dictionary = dictionary_mayo,nCores = 4)

### Second Step
        
t_all <- transcriptogramStep2(object = t_all, nCores = 4 )
saveRDS(t_all,'deseq/MAYO/t_all.rds')

t_A <- transcriptogramStep2(object = t_A, nCores = 4 )
saveRDS(t_A,'deseq/MAYO/t_A.rds')

t_B <- transcriptogramStep2(object = t_B, nCores = 4 )
saveRDS(t_B,'deseq/MAYO/t_B.rds')

t_C <- transcriptogramStep2(object = t_C, nCores = 4 )
saveRDS(t_C,'deseq/MAYO/t_C.rds')


mayo_metadadoA <- mayo_metadado[mayo_metadado$AgeAtDeath == 'A',]
mayo_metadadoB <- mayo_metadado[mayo_metadado$AgeAtDeath == 'B',]
mayo_metadadoC <- mayo_metadado[mayo_metadado$AgeAtDeath == 'C',]


### Third Step

t_all <- differentiallyExpressed(object = t_all,
         levels = ifelse(mayo_metadado$condition == 'Control',TRUE,FALSE)
          , pValue = 0.001,
          trend = F, title = "All Groups",DEsymbols)

t_A <- differentiallyExpressed(object = t_A,
    levels = ifelse(mayo_metadadoA$condition == 'Control',TRUE,FALSE), pValue = 0.001,
    trend = F, title = "All Groups",DEsymbols) #### No clusters differentially expressed


t_B <- differentiallyExpressed(object = t_B,
    levels = ifelse(mayo_metadadoB$condition == 'Control',TRUE,FALSE), pValue = 0.001,
    trend = F, title = "Group B",DEsymbols) 

t_C <- differentiallyExpressed(object = t_C,
      levels = ifelse(mayo_metadadoC$condition == 'Control',TRUE,FALSE), pValue = 0.001,
      trend = F, title = "Group C",DEsymbols)




###Cluster enrinchement
t_all_cluster <- clusterEnrichment(object = t_all, species = 'Homo sapiens',
                       pValue = 0.001, nCores = 6)
saveRDS(t_all_cluster,'deseq/MAYO/all_cluster.rds')

t_B_cluster <- clusterEnrichment(object = t_B, species = 'Homo sapiens',
                pValue = 0.001, nCores = 6)
saveRDS(t_B_cluster,'deseq/MAYO/t_B_cluster.rds')

t_C_cluster <- clusterEnrichment(object = t_C, species = 'Homo sapiens',
                                 pValue = 0.001, nCores = 6)
saveRDS(t_C_cluster,'deseq/MAYO/t_C_cluster.rds')

#### Clusters and Graph

all_clusters <- t_all_cluster@Terms

all_properties <- connectivityProperties(t_all)

B_clusters <- t_B_cluster@Terms

C_clusters <- t_C_cluster@Terms

rdp <- clusterVisualization(object = t_B,clusters = c(2,9:13),
                            onlyGenesInDE = T,maincomp = T,connected = T )


#### DE genes
DE_all <- DE(object = t_all)

DE_B <- DE(object = t_B)

DE_C <- DE(object = t_C)


### Graph

library('igraph')
g <- igraph::graph.data.frame(d = t_all@association)


gt3  <- subg(g=g, dat=DE_all[DE_all$ClusterNumber == 12,], refcol=1)
gt3  <- att.setv(g=gt3, from="Symbol", to="nodeAlias")
rdp <- RedPort() 
calld(rdp)
addGraph(rdp, gt3, gcoord=c(10,25), gscale=20, isNest=TRUE, theme='tm1', zoom=30)
nestNodes(rdp, nodes=V(gt3)$name, parent="N1", theme='tm2')
resetd(rdp)
