### Import of DATA
library('igraph')
library('RedeR')
library('transcriptogramer')
library('RColorBrewer')
library('dplyr')

#### This is script is to construct the network in Reder with
### the informations of transcriptogramer analysis.

ALZ_PSP_DEG <- readRDS('trans_results/ALZ_PSP_DEG.rds')
CCX_HIP_DEG <- readRDS('trans_results/CCX_HIP_DEG.rds')
CCX_all <- readRDS('trans_results/FAD/clusters/all_cluster_CCX.rds')
dictionary_FAD <- read.csv('results/tables/dictionary_FAD.csv',
                  header = T, stringsAsFactors = F,row.names = 1)

ALZ_PSP_DEG$Age_Group <- paste(ALZ_PSP_DEG$Group,ALZ_PSP_DEG$Age_Group,
                          sep='-')
CCX_HIP_DEG$Month_Group <- paste(CCX_HIP_DEG$Region,
                            CCX_HIP_DEG$Month_Group,sep='-')
dictionary_FAD <- dictionary_FAD[dictionary_FAD$ensembl_peptide_id %in%
                  CCX_HIP_DEG$Protein,]
dictionary_FAD <- dictionary_FAD[,-2]
colnames(dictionary_FAD) <- c('Protein','Symbol')

CCX_HIP_DEG <- CCX_HIP_DEG[,-8] ### exclude the symbol column with EMSP

CCX_HIP_DEG <- merge(CCX_HIP_DEG,dictionary_FAD,by='Protein')

## Human


association_ABC <- association

association_ABC <- association_ABC[association_ABC$protein1 %in%
                                     ALZ_PSP_DEG$Protein |
                                     association_ABC$protein2 %in% 
                                     ALZ_PSP_DEG$Protein,]             

# graph
g <- igraph::graph.data.frame(d=association_ABC,
                 directed = F)

g_ABC <- RedeR::subg(g=g, dat = ALZ_PSP_DEG,refcol = 1,
              maincomp = T,connected = T)
                    
g_ABC <- att.setv(g=g_ABC,from='Symbol',to='nodeAlias')
g_ABC <-att.setv(g=g_ABC,from='Age_Group', to='nodeColor',
                 cols = brewer.pal(n = 6, name = "Accent"))

addGraph(rdp, g_ABC,gscale=30,zoom=30)

relax(rdp,50,400)
scl <- g_ABC$legNodeColor$scale 
leg <- g_ABC$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")



resetd(rdp)
graph_propriets <- as_data_frame(g_ABC,what = 'both')

write.csv(graph_propriets$vertices,
          'results/transcriptogramer/ALZ_PSP_V.csv')
write.csv(graph_propriets$edges,
          'results/transcriptogramer/ALZ_PSP_E.csv')
saveRDS(g_ABC,'results/transcriptogramer/ALZ_PSP_G.rds')

#### Mouse 

association_mouse <- CCX_all@association

association_mouse <- association_mouse[association_mouse$p1 %in%
                                      CCX_HIP_DEG$Protein |
                                      association_mouse$p2 %in%
                                      CCX_HIP_DEG$Protein,] 

g <- igraph::graph.data.frame(d=association_mouse,
                              directed = F)

g_mouse <- RedeR::subg(g=g, dat = CCX_HIP_DEG,refcol = 1,
                     maincomp = T,connected = T)

g_mouse <- att.setv(g=g_mouse,from='Symbol',to='nodeAlias')
g_mouse <-att.setv(g=g_mouse,from='Month_Group', to='nodeColor',
                 cols = brewer.pal(n = 7, name = "Accent"))


### RedeR

rdp <- RedPort() 
calld(rdp)


addGraph(rdp, g_mouse,gscale=30,zoom=30)

relax(rdp,50,400)
scl <- g_mouse$legNodeColor$scale 
leg <- g_mouse$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")
resetd(rdp)

graph_propriets_mouse <- as_data_frame(g_mouse,what = 'both')

write.csv(graph_propriets_mouse$vertices,
          'results/transcriptogramer/CCX_HIP_V.csv')
write.csv(graph_propriets_mouse$edges,
          'results/transcriptogramer/CCX_HIP_E.csv')
saveRDS(g_mouse,'results/transcriptogramer/CCX_HIP_G.rds')



#### Human

ALZ_PSP_E <- read.csv('results/transcriptogramer/ALZ_PSP_E.csv',
              header = T,stringsAsFactors = F)
ALZ_PSP_V <- read.csv('results/transcriptogramer/ALZ_PSP_V.csv',
                    header = T,stringsAsFactors = F)

ALZ_PSP_from <- data.frame(Protein=ALZ_PSP_E$Protein,stringsAsFactors = F)
ALZ_PSP_from <- merge(ALZ_PSP_from,ALZ_PSP_V,by='Protein')

ALZ_PSP_to <- data.frame(Protein=ALZ_PSP_E$Protein.1,stringsAsFactors = F)
ALZ_PSP_to <- merge(ALZ_PSP_to,ALZ_PSP_V,by='Protein')

ALZ_PSP_from <- ALZ_PSP_from[match(ALZ_PSP_E$Protein,ALZ_PSP_from$Protein),]
ALZ_PSP_to <- ALZ_PSP_to[match(ALZ_PSP_E$Protein.1,ALZ_PSP_to$Protein),]

ALZ_PSP_E$Protein <- ALZ_PSP_from$Symbol
ALZ_PSP_E$Protein.1 <- ALZ_PSP_to$Symbol
colnames(ALZ_PSP_E) <- c('from','to')
write.csv(ALZ_PSP_E,'results/transcriptogramer/ALZ_PSP_E.csv')

### Mouse
CCX_HIP_E <- read.csv('results/transcriptogramer/CCX_HIP_E.csv',
                      header = T,stringsAsFactors = F)
CCX_HIP_V <- read.csv('results/transcriptogramer/CCX_HIP_V.csv',
                      header = T,stringsAsFactors = F)

CCX_HIP_E_from <- data.frame(Protein=CCX_HIP_E$from,stringsAsFactors = F)
colnames(CCX_HIP_E_from) <- 'Protein'
CCX_HIP_E_from <- merge(CCX_HIP_E_from,CCX_HIP_V,by='Protein')

CCX_HIP_E_to <- data.frame(Protein=CCX_HIP_E$to,stringsAsFactors = F)
colnames(CCX_HIP_E_to) <- 'Protein'
CCX_HIP_E_to <- merge(CCX_HIP_E_to,CCX_HIP_V,by='Protein')

CCX_HIP_E_from <- CCX_HIP_E_from[match(CCX_HIP_E$from,
                  CCX_HIP_E_from$Protein),]
CCX_HIP_E_to <- CCX_HIP_E_to[match(CCX_HIP_E$to,
                  CCX_HIP_E_to$Protein),]

CCX_HIP_E$from <- CCX_HIP_E_from$Symbol
CCX_HIP_E$to <- CCX_HIP_E_to$Symbol

write.csv(CCX_HIP_E,'results/transcriptogramer/CCX_HIP_E.csv')
