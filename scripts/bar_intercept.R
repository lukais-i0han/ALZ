### Intercept - Bar graphics - Up down
library('Vennerable')
library("UpSetR")
library('dplyr')
library('stringr')
library('tidyr')
library('ggplot2')
## Human

ALZ_sig <- read.csv('archives/ALZ_sig.csv', header = T, stringsAsFactors = F,
           row.names = 1)

names(ALZ_sig)[names(ALZ_sig) == 'AgeAtDeath'] <- 'Age_Group'

ALZ_sig <- ALZ_sig %>% mutate(DEG = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))


PSP_sig <- read.csv('archives/PSP_sig.csv',header = T, stringsAsFactors = F,
            row.names = 1)

names(PSP_sig)[names(PSP_sig) == 'AgeAtDeath'] <- 'Age_Group'

PSP_sig <- PSP_sig %>% mutate(DEG = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))


#####


a<-ggplot(ALZ_sig,aes(x=Age_Group, y=factor(Gene_Symbol, levels = unique((Gene_Symbol))), fill=DEG)) +
  geom_tile(size=0.3,width = 1)  +# coord_equal() +
  labs(y='Genes', x=NULL,fill = 'DEG') + scale_fill_manual(values = c('firebrick','blue'),labels = c('Up','Down'))+
  theme_classic() + theme(strip.background = element_blank(),
                          panel.grid.major.x = element_line(color='black',size = 0.5),
                          axis.text.y = element_blank(),
                          legend.text = element_text(size = 12,face = 'bold'),
                          legend.justification = 'top',
                          legend.box.just = 'right')
                          
a  

b <- ggplot(PSP_sig,aes(x=Age_Group, y=factor(Gene_Symbol, levels = unique((Gene_Symbol))), fill=DEG)) +
  geom_tile(size=0.3,width = 1)  +# coord_equal() +
  labs(y='Genes', x=NULL,fill = 'DEG') + scale_fill_manual(values = c('firebrick','blue'),labels = c('Up','Down')) +
  theme_classic() + theme(strip.background = element_blank(),
                          panel.grid.major.x = element_line(color='black',size = 0.5),
                          strip.text.x = element_text(size = 10,face = 'bold'), plot.title = element_text(size = 10),
                          axis.text.y = element_blank(),
                          legend.text = element_text(size = 12,face = 'bold'),
                          legend.justification = 'top',
                          legend.box.just = 'right')
b


png(file= 'results/figures/up_donw_ALZ.png',height = 1200,width = 1000,res = 200)
a
dev.off()

png(file= 'results/figures/up_donw_psp.png',height = 1200,width = 900,res = 200)
b
dev.off()

### Mouse

FAD5X <- read.csv('archives/FAD5X_sig.csv',header = T, stringsAsFactors = F,
                  row.names = 1)

FAD5X$Month <- factor(FAD5X$Month,levels = c('4','12','18'))

FAD5X <- FAD5X %>% mutate(DEG = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))

Tau35 <- read.csv('archives/Tau_sig.csv',header = T, stringsAsFactors = F,
        row.names = 1)

Tau35$Month <- factor(Tau35$Month,levels = c('4','17'))

Tau35 <- Tau35 %>% mutate(DEG = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))


a<-ggplot(FAD5X,aes(x=Month, y=factor(Gene_Symbol, levels = unique((Gene_Symbol))), fill=DEG)) +
  geom_tile(size=0.3,width = 1)  +
  facet_grid(~Area, scales = "free", space = "free")+# coord_equal() +
  labs(y='Genes', x=NULL,fill = 'DEG') + scale_fill_manual(values = c('firebrick','blue'),labels = c('Up','Down'))+
  theme_classic() + theme(strip.background = element_blank(),
                          panel.grid.major.x = element_line(color='black',size = 0.5),
                          axis.text.y = element_blank(),
                          legend.text = element_text(size = 12,face = 'bold'),
                          legend.justification = 'top',
                          legend.box.just = 'right')

a

png(file= 'results/figures/up_donw_FAD5X.png',height = 1200,width = 1000,res = 200)
a
dev.off()

b <- ggplot(Tau35,aes(x=Month, y=factor(Gene_Name, levels = unique((Gene_Name))), fill=DEG)) +
  geom_tile(size=0.3,width = 1)  +# coord_equal() +
  labs(y='Genes', x=NULL,fill = 'DEG') + scale_fill_manual(values = c('firebrick','blue'),labels = c('Up','Down')) +
  theme_classic() + theme(strip.background = element_blank(),
                          panel.grid.major.x = element_line(color='black',size = 0.5),
                          strip.text.x = element_text(size = 10,face = 'bold'), plot.title = element_text(size = 10),
                          axis.text.y = element_blank(),
                          legend.text = element_text(size = 12,face = 'bold'),
                          legend.justification = 'top',
                          legend.box.just = 'right')
b


png(file= 'results/figures/up_donw_tau.png',height = 1200,width = 1000,res = 200)
b
dev.off()


### Intercept

mayo_A <- mayo_sig[mayo_sig$Age_Group == 'A',]
mayo_B <- mayo_sig[mayo_sig$Age_Group == 'B',]
mayo_C <- mayo_sig[mayo_sig$Age_Group == 'C',]

PSP_A <- PSP_sig[PSP_sig$Age_Group == 'A',]
PSP_B <- PSP_sig[PSP_sig$Age_Group == 'B',]


mayo_list <- list(A = unique(mayo_A$Gene_Symbol),
                  B = unique(mayo_B$Gene_Symbol),
                  C = unique(mayo_C$Gene_Symbol))

intercep_mayo <- venn::venn(mayo_list,
    zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)
inters_AD <-attr(intercep_mayo,"intersections")

## 
PSP_list <- list(A=unique(PSP_A$Gene_Symbol),
                  B = unique(PSP_B$Gene_Symbol))

intercept_PSP <-venn::venn(PSP_list,
    zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)
inters_PSP <-attr(intercept_PSP,'intersections') 

inters <- plyr::ldply(inters_AD, rbind)
inters <- t(inters)
inters <- data.frame(inters,stringsAsFactors = F)
colnames(inters) <- inters[1,]
inters <- inters[-1,]
inters[is.na(inters)] <-""

saveRDS(inters,'results/mayo_intersection.rds')

inters_PSP <- plyr::ldply(inters_PSP, rbind)
inters_PSP <- t(inters_PSP)
inters_PSP <- data.frame(inters_PSP,stringsAsFactors = F)
colnames(inters_PSP) <- inters_PSP[1,]
inters_PSP <- inters_PSP[-1,]
inters_PSP[is.na(inters_PSP)] <-""

saveRDS(inters_PSP,'results/PSP_intersection.rds')

### Mouse

##  FAD5X

FAD5X <- read.csv('archives/FAD5X_sig.csv',header = T, stringsAsFactors = F,
         row.names = 1 )

FAD5X$Month <- as.character(FAD5X$Month)

FAD5X <- FAD5X %>% mutate(DEG = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))

FAD_4CCX <- FAD5X[FAD5X$Month == '4' & FAD5X$Area == 'Cortex',]
FAD_12CCX <- FAD5X[FAD5X$Month == '12' & FAD5X$Area == 'Cortex',]
FAD_18CCX <- FAD5X[FAD5X$Month == '18' & FAD5X$Area == 'Cortex',]

FAD_4HIP <- FAD5X[FAD5X$Month == '4' & FAD5X$Area == 'Hippocampus',]
FAD_12HIP <- FAD5X[FAD5X$Month == '12' & FAD5X$Area == 'Hippocampus',]
FAD_18HIP <- FAD5X[FAD5X$Month == '18' & FAD5X$Area == 'Hippocampus',]

FADCCX_list <- list(Month_4 = unique(FAD_4CCX$Gene_Symbol),
                  Month_12 = unique(FAD_12CCX$Gene_Symbol),
                  Month_18 = unique(FAD_18CCX$Gene_Symbol))


FADHIP_list <- list(Month_4 = unique(FAD_4HIP$Gene_Symbol),
                    Month_12 = unique(FAD_12HIP$Gene_Symbol),
                    Month_18 = unique(FAD_18HIP$Gene_Symbol))

intercep_FAD_CCX <- venn::venn(FADCCX_list,
        zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)
inters_FADCXX <-attr(intercep_FAD_CCX,"intersections")

inters_FADCXX <- plyr::ldply(inters_FADCXX, rbind)
inters_FADCXX <- t(inters_FADCXX)
inters_FADCXX <- data.frame(inters_FADCXX,stringsAsFactors = F)
colnames(inters_FADCXX) <- inters_FADCXX[1,]
inters_FADCXX <- inters_FADCXX[-1,]
inters_FADCXX[is.na(inters_FADCXX)] <-""

saveRDS(inters_FADCXX,'results/FAD5_CCX_intersection.rds')

intercep_FAD_HIP <- venn::venn(FADHIP_list,
zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)
inters_FADHIP <-attr(intercep_FAD_HIP,"intersections")

inters_FADHIP <- plyr::ldply(inters_FADHIP, rbind)
inters_FADHIP <- t(inters_FADHIP)
inters_FADHIP <- data.frame(inters_FADHIP,stringsAsFactors = F)
colnames(inters_FADHIP) <- inters_FADHIP[1,]
inters_FADHIP <- inters_FADHIP[-1,]
inters_FADHIP[is.na(inters_FADHIP)] <-""

saveRDS(inters_FADHIP,'results/FAD5_HIP_intersection.rds')

## Tau

Tau <- read.csv('archives/Tau_sig.csv',header = T, stringsAsFactors = F,
        row.names = 1)

Tau_4 <-Tau[Tau$Month == '4',]
Tau_17 <- Tau[Tau$Month == '17',]

Tau_list <- list(Month_4 = unique(Tau_4$Gene_Name),
                 Month_17=unique(Tau_17$Gene_Name))


intercept_Tau <- venn::venn(Tau_list,
zcolor = 'style',cexil = 1.5,cexsn = 1,ilabels = F,borders = F,size = 5)

inters_Tau <-attr(intercept_Tau,"intersections")

inters_Tau <- plyr::ldply(inters_Tau, rbind)
inters_Tau <- t(inters_Tau)
inters_Tau <- data.frame(inters_Tau,stringsAsFactors = F)
colnames(inters_Tau) <- inters_Tau[1,]
inters_Tau <- inters_Tau[-1,]
inters_Tau[is.na(inters_Tau)] <-""

saveRDS(inters_Tau,'results/Tau_intersection.rds')
