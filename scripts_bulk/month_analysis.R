#### Month and Region analysis
library('dplyr')
library('stringr')
library('ggplot2')
library('cowplot')

## Human

ALZ <- read.csv('archives/ALZ_all.csv',header = T, 
                 stringsAsFactors = F,
                 row.names = 1)

names(ALZ)[names(ALZ) == 'AgeAtDeath'] <- 'Age_Group'

PSP <- read.csv('archives/PSP_all.csv',header = T, 
                 stringsAsFactors = F,
                 row.names = 1)

names(PSP)[names(PSP) == 'AgeAtDeath'] <- 'Age_Group'

### GGplot2 theme

theme_classic2 <- function(base_size = 15, base_family = ""){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border     = element_blank(),
      panel.background = element_rect(color = 'black',size = 0.5),
      axis.line        = element_line(colour = "black"),
      panel.grid.minor.y = element_line(size = 0.5,color='black',linetype = 'solid'),
      panel.grid.minor = element_line(size = 0.5,colour = 'black',linetype = 'solid'),
      panel.grid.minor.x = element_line(size=0.5,colour = 'black',linetype = 'solid'),
      strip.background = element_blank(),
      legend.key       = element_blank()
      
    )
}


j = ggplot2::ggplot(ALZ, ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

g = j + ggrepel::geom_text_repel(data=ALZ,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(legend.position = 'none',
    axis.text.x = element_text(face = "bold", color = "black",size = 15),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
    
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-8,8))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](padj))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
g

h = ggplot2::ggplot(PSP, ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

i = h + ggrepel::geom_text_repel(data=PSP,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 15),
    axis.text.x = element_text(face = "bold", color = "black",size = 15),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
    
  )   +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-8,8))+
  labs(x= expression(log[2](FC)),y= '') +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
i


png(file= 'results/figures/month_human.png',height = 1500,width = 3500,res = 200)

plot_grid(g,i,labels = c('ALZ','PSP'),label_size = 10)

dev.off()

### Mouse

FAD5X <- read.csv('archives/FAD5X_all.csv',header = T, 
         stringsAsFactors = F,
         row.names = 1)
Tau35 <- read.csv('archives/Tau_all.csv',header = T,
         stringsAsFactors = F,
         row.names = 1)


Tau35$Month <- factor(Tau35$Month,levels = c('4','17'))

Tau35 <- Tau35 %>% mutate(UP_DOWN = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))
Tau35 <- Tau35[Tau35$DEG != 'Non_Significant',]


FAD5X$Month <- factor(FAD5X$Month, levels = c('4','12','18'))
FAD5X_CCX <- FAD5X[FAD5X$Area == 'Cortex',]

FAD5X_CCX <- FAD5X_CCX %>% mutate(UP_DOWN = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))
FAD5X_CCX <- FAD5X_CCX[FAD5X_CCX$DEG != 'Non_Significant',]

FAD5X_HIP <- FAD5X[FAD5X$Area == 'Hippocampus',]

FAD5X_HIP <- FAD5X_HIP %>% mutate(UP_DOWN = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))
FAD5X_HIP <- FAD5X_HIP[FAD5X_HIP$DEG != 'Non_Significant',]

j = ggplot2::ggplot(FAD5X_CCX, ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

g = j + ggrepel::geom_text_repel(data=FAD5X_CCX,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(legend.position = 'none',
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Month,ncol = 3) + ylim(0,20)+ xlim(c(-8,8))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](padj))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
g

h = ggplot2::ggplot(FAD5X_HIP, ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

i = h + ggrepel::geom_text_repel(data=FAD5X_HIP,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 15),
    axis.text.x = element_text(face = "bold", color = "black",size = 15),
    axis.title = element_text(face='bold',color='black',size=15),
    axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
    
  )   +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Month,ncol = 3) + ylim(0,20)+ xlim(c(-8,8))+
  labs(x= expression(log[2](FC)),y= '') +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
i

png(file= 'results/figures/Month_FAD5X.png',height = 1500,width = 3500,res = 200)

plot_grid(g,i,labels = c('Cortex','Hippocampus'),label_size = 10)

dev.off()


j = ggplot2::ggplot(Tau35, ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

g = j + ggrepel::geom_text_repel(data=,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Month,ncol = 2) + ylim(0,20)+ xlim(c(-8,8))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](padj))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
g

png(file= 'results/figures/Month_TAU.png',height = 1500,width = 3500,res = 200)

plot_grid(g,labels = 'Hippocampus',label_size = 10)

dev.off()
