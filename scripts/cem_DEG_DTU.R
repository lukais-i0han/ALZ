library('RColorBrewer')
library('stringr')
library('ggplot2')
library('readr')
library('dplyr')

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

colnames(human_modules_DEG) <- c('Modules','Age_Group','Group','Number of Genes')

human_modules_DEG$Modules <- factor(human_modules_DEG$Modules,
              levels = c('M1','M2','M3','M4','M5'))

human_modules_DEG$Age_Group <- factor(human_modules_DEG$Age_Group,
              levels = c('[A]','[A-B]','[A-C]','[A-B-C]','[B]','[B-C]','[C]'))

human_DEG_plot <-ggplot(human_modules_DEG,aes(x=Modules, y=Value, fill=Age_Group))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Group) 

human_modules_DTU <- human_modules[human_modules$Class == 'DTU',]
human_modules_DTU <- human_modules_DTU[order(human_modules_DTU$modules),]
human_modules_DTU <- table(human_modules_DTU$modules,human_modules_DTU$Age_Group,human_modules_DTU$Group)
human_modules_DTU <- data.frame(human_modules_DTU,stringsAsFactors = F)
human_modules_DTU <- human_modules_DTU[human_modules_DTU$Freq != 0,]

colnames(human_modules_DTU) <- c('Modules','Age_Group','Group','Number of Genes')

human_modules_DTU$Modules <- factor(human_modules_DTU$Modules,
                                    levels = c('M1','M2','M3','M4','M5'))

human_modules_DTU$Age_Group <- factor(human_modules_DTU$Age_Group,
          levels = c('[A]','[A-B]','[A-C]','[A-B-C]','[B]','[B-C]','[C]'))

human_DTU_plot <- ggplot(human_modules_DTU,aes(x=Modules, y=Value, fill=Age_Group))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Group) 



human_ora <- human_cem$ora
human_ora <- human_ora[str_detect(human_ora$ID,'%GO:'),]
ID_human_ora <- str_split_fixed(human_ora$ID,'%',3)
ID_human_ora <- data.frame(ID_human_ora)

human_ora$ID <- ID_human_ora$X2
human_ora$Description <- ID_human_ora$X1
human_ora <- human_ora[human_ora$Module != 'Not.Correlated',]
human_ora <- human_ora[human_ora$p.adjust < 0.01,]

human_ora <- human_ora %>%
  group_by(Module) %>%
  top_n(n=10,wt = Count)

human_ora <- human_ora[order(human_ora$Description,decreasing = T),]
human_ora$Description <- factor(human_ora$Description,levels = unique(human_ora$Description))


ora_human_plot <- ggplot(human_ora,aes(x=Description, y=Count))+
  geom_bar(stat="identity")  +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Module) 

## Mouse

FAD_cem <- list(modules= read_tsv('CEMtool/results/Tables/FAD//module.tsv'),
                ora = read_tsv('CEMtool/results/Tables/FAD/ora.tsv'),
                ppi = read_tsv('CEMtool/results/Tables/FAD/interactions.tsv'))


FAD_modules <- FAD_cem$modules

F_DEG_DTU <- F_DEG_DTU[F_DEG_DTU$genes %in% FAD_modules$genes,]
FAD_modules <- merge(FAD_modules,F_DEG_DTU,by='genes')

FAD_modules <- FAD_modules[FAD_modules$modules != 'Not.Correlated',]

FAD_modules_CCX <- FAD_modules[FAD_modules$Area == 'Cortex',]
FAD_modules_CCX <- FAD_modules_CCX[order(FAD_modules_CCX$modules),]
FAD_modules_CCX <- table(FAD_modules_CCX$modules,FAD_modules_CCX$Age_Group,
                         FAD_modules_CCX$Class)
FAD_modules_CCX<- data.frame(FAD_modules_CCX,stringsAsFactors = F)
FAD_modules_CCX <- FAD_modules_CCX[FAD_modules_CCX$Freq != 0,]

colnames(FAD_modules_CCX) <- c('Modules','Age_Group','Group','No.Genes')

FAD_modules_CCX$Modules <- factor(FAD_modules_CCX$Modules,
                                    levels = c('M1','M2','M3','M4','M5'))

FAD_modules_CCX$Age_Group <- factor(FAD_modules_CCX$Age_Group,
          levels = c('[4]','[4-12]','[4-18]','[4-12-18]','[12]','[12-18]','[18]'))

CCX_FAD <-ggplot(FAD_modules_CCX,aes(x=Modules, y=No.Genes, fill=Age_Group))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Group) 

#HIPPOCAMPUS

FAD_modules_HIP <- FAD_modules[FAD_modules$Area == 'Hippocampus',]
FAD_modules_HIP <- FAD_modules_HIP[order(FAD_modules_HIP$modules),]
FAD_modules_HIP <- table(FAD_modules_HIP$modules,FAD_modules_HIP$Age_Group,
                         FAD_modules_HIP$Class)
FAD_modules_HIP<- data.frame(FAD_modules_HIP,stringsAsFactors = F)
FAD_modules_HIP <- FAD_modules_HIP[FAD_modules_HIP$Freq != 0,]

colnames(FAD_modules_HIP) <- c('Modules','Age_Group','Group','No.Genes')

FAD_modules_HIP$Modules <- factor(FAD_modules_HIP$Modules,
         levels = c('M1','M2','M3','M4','M5','M6'))

FAD_modules_HIP$Age_Group <- factor(FAD_modules_HIP$Age_Group,
              levels = c('[4]','[4-12]','[4-18]','[4-12-18]','[12]','[12-18]','[18]'))

HIP_FAD<-ggplot(FAD_modules_HIP,aes(x=Modules, y=No.Genes, fill=Age_Group))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Group)  



FAD_ora <- FAD_cem$ora
FAD_ora <- FAD_ora[str_detect(FAD_ora$ID,'%GO:'),]
ID_FAD_ora <- str_split_fixed(FAD_ora$ID,'%',3)
ID_FAD_ora <- data.frame(ID_FAD_ora)

FAD_ora$ID <- ID_FAD_ora$X2
FAD_ora$Description <- ID_FAD_ora$X1
FAD_ora <- FAD_ora[FAD_ora$Module != 'Not.Correlated',]
FAD_ora <- FAD_ora[FAD_ora$p.adjust < 0.01,]

FAD_ora <- FAD_ora %>%
  group_by(Module) %>%
  top_n(n=15,wt = Count)

FAD_ora <- FAD_ora[order(FAD_ora$Description,decreasing = T),]
FAD_ora$Description <- factor(FAD_ora$Description,levels = unique(FAD_ora$Description))


ora_FAD_plot<-ggplot(FAD_ora,aes(x=Description, y=Count))+
  geom_bar(stat="identity")  +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Module) 








####


TAU_cem <- list(modules= read_tsv('CEMtool/results/Tables/TAU/module.tsv'),
                ora = read_tsv('CEMtool/results/Tables/TAU/ora.tsv'),
                ppi = read_tsv('CEMtool/results/Tables/TAU/interactions.tsv'))

Tau_modules <- TAU_cem$modules

T_DEG_DTU <- T_DEG_DTU[T_DEG_DTU$genes %in% Tau_modules$genes,]
Tau_modules <- merge(Tau_modules,T_DEG_DTU,by='genes')

Tau_modules <- Tau_modules[Tau_modules$modules != 'Not.Correlated',]

Tau_modules <- Tau_modules[order(Tau_modules$modules),]
Tau_modules <- table(Tau_modules$modules,Tau_modules$Age_Group,
                     Tau_modules$Class)

Tau_modules <- data.frame(Tau_modules,stringsAsFactors = F)
Tau_modules <- Tau_modules[Tau_modules$Freq != 0,]

colnames(Tau_modules) <- c('Modules','Age_Group','Group','No.Genes')

Tau_modules$Modules <- factor(Tau_modules$Modules,
                                  levels = c('M1','M2','M3','M4','M7'))


Tau_modules$Age_Group <- factor(Tau_modules$Age_Group,
      levels = c('[4]','[4-17]','[17]'))

HIP_TAU<-ggplot(Tau_modules,aes(x=Modules, y=No.Genes, fill=Age_Group))+
  geom_bar(stat="identity") + scale_fill_brewer(palette = 'Accent') +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Group)  


Tau_ora <- TAU_cem$ora
Tau_ora <- Tau_ora[str_detect(Tau_ora$ID,'%GO:'),]
ID_Tau_ora <- str_split_fixed(Tau_ora$ID,'%',3)
ID_Tau_ora <- data.frame(ID_Tau_ora)

Tau_ora$ID <- ID_Tau_ora$X2
Tau_ora$Description <- ID_Tau_ora$X1
Tau_ora <- Tau_ora[Tau_ora$Module != 'Not.Correlated',]
Tau_ora <- Tau_ora[Tau_ora$p.adjust < 0.01,]

Tau_ora <- Tau_ora %>%
  group_by(Module) %>%
  top_n(n=10,wt = Count)

Tau_ora <- Tau_ora[order(Tau_ora$Description,decreasing = T),]
Tau_ora$Description <- factor(Tau_ora$Description,levels = unique(Tau_ora$Description))


ora_Tau_plot<-ggplot(Tau_ora,aes(x=Description, y=Count))+
  geom_bar(stat="identity")  +
  coord_flip() +theme_minimal()+ 
  theme(
    axis.text.y = element_text(face = 'bold'),
    axis.text.x = element_text(face='bold'),legend.position = 'top')+
  facet_grid(~Module) 

