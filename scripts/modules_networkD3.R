library('RColorBrewer')
library('stringr')
library('ggplot2')
library('readr')
library('dplyr')
library('networkD3')



useRtreeList <- ToListExplicit(population,unname = TRUE)
radialNetwork( useRtreeList,fontSize = 10)


### Analysis of DEG-DTU and CEMitool results (Ora plots)

H_DEG_DTU <- readRDS('CEMtool/results/human/H_DEG_DTU.rds')

F_DEG_DTU <- readRDS('CEMtool/results/mouse_fad/F_DEG_DTU.rds')

T_DEG_DTU <- readRDS('CEMtool/results/mouse_tau/T_DEG_DTU.rds')


# Human


human_cem <- list(modules= read_tsv('CEMtool/results/Tables/Human/module.tsv'),
                  ora = read_tsv('CEMtool/results/Tables/Human/ora.tsv'),
                  ppi = read_tsv('CEMtool/results/Tables/Human/interactions.tsv'))


human_modules <- human_cem$modules

H_DEG_DTU <- H_DEG_DTU[H_DEG_DTU$genes %in% human_modules$genes,]
human_modules <- merge(human_modules,H_DEG_DTU,by='genes')

human_modules <- human_modules[human_modules$modules != 'Not.Correlated',]

human_modules_DEG <- human_modules[human_modules$Class == 'DEG',]
human_modules_DEG <- human_modules_DEG[order(human_modules_DEG$modules),]
human_modules_DEG <- table(human_modules_DEG$modules,human_modules_DEG$Age_Group,human_modules_DEG$Group)
human_modules_DEG <- data.frame(human_modules_DEG,stringsAsFactors = F)
human_modules_DEG <- human_modules_DEG[human_modules_DEG$Freq != 0,]

colnames(human_modules_DEG) <- c('modules','Age_Group','Group','Number.of.Genes')

human_modules_DEG$modules <- factor(human_modules_DEG$modules,
                                    levels = c('M1','M2','M3','M4','M5'))

human_modules_DEG$Age_Group <- factor(human_modules_DEG$Age_Group,
                                      levels = c('[A]','[A-B]','[A-C]','[A-B-C]','[B]','[B-C]','[C]'))

human_modules_DTU <- human_modules[human_modules$Class == 'DTU',]
human_modules_DTU <- human_modules_DTU[order(human_modules_DTU$modules),]
human_modules_DTU <- table(human_modules_DTU$modules,human_modules_DTU$Age_Group,human_modules_DTU$Group)
human_modules_DTU <- data.frame(human_modules_DTU,stringsAsFactors = F)
human_modules_DTU <- human_modules_DTU[human_modules_DTU$Freq != 0,]

colnames(human_modules_DTU) <- c('modules','Age_Group','Group','Number.of.Genes')

human_modules_DTU$modules <- factor(human_modules_DTU$modules,
                                    levels = c('M1','M2','M3','M4','M5'))

human_modules_DTU$Age_Group <- factor(human_modules_DTU$Age_Group,
                levels = c('[A]','[A-B]','[A-C]','[A-B-C]','[B]','[B-C]','[C]'))

human_modules_DEG$'Class' <- 'DEG'
human_modules_DTU$'Class' <- 'DTU'

human_modules_DEG_DTU <- rbind(human_modules_DEG,human_modules_DTU)

human_modules_DEG_DTU <- human_modules_DEG_DTU[
human_modules_DEG_DTU$Group == 'PSP',]

human_modules_DEG_DTU$pathString <- paste('Modules',
                              human_modules_DEG_DTU$modules,
                              human_modules_DEG_DTU$Group,
                              human_modules_DEG_DTU$Class,
                              human_modules_DEG_DTU$Age_Group,
                              human_modules_DEG_DTU$Number.of.Genes,
                              sep = '/')

nodes_human_modules <- as.Node(human_modules_DEG_DTU)


useRtreeList <- ToListExplicit(nodes_human_modules,unname = TRUE)
radialNetwork( useRtreeList,fontSize = 10,nodeStroke = 'blue',opacity = 5)
diagonalNetwork(List = useRtreeList, fontSize = 10, opacity = 0.9)
