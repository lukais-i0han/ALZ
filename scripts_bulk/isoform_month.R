### Region analysis
library('stringr')
library('dplyr')
library('ggplot2')
library('cowplot')

### Human

mayo <- readRDS('results/isoform_results/mayo_ISW.rds')

mayo <- mayo %>% mutate(Group = case_when(
  condition_1 == 'AD' & condition_2 == 'Control' ~ 'ALZ',
  TRUE ~ 'PSP'
))

mayo_ALZ <- mayo[mayo$Group == 'ALZ',]
mayo_PSP <- mayo[mayo$Group == 'PSP',]

theme_classic2 <- function(base_size = 15, base_family = ""){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border     = element_blank(),
      panel.background = element_rect(color = 'black',size = 0.5),
      axis.line        = element_line(colour = "black"),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor = element_line(size = 0.5,colour = 'black',linetype = 'solid'),
      panel.grid.minor.x = element_blank(),
      strip.background = element_blank(),
      legend.key       = element_blank()
      
    )
}


j = ggplot2::ggplot(mayo_ALZ, ggplot2::aes(dIF,-log10(iso_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

g = j + ggrepel::geom_text_repel(data=mayo_ALZ,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-0.05,0.05),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(legend.position = 'none',
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF(ALZ/Control)),y= expression(-log[10](FDR))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
g

h = ggplot2::ggplot(mayo_PSP, ggplot2::aes(dIF,-log10(iso_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

i = h + ggrepel::geom_text_repel(data=mayo_PSP,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-0.05,0.05),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Age_Group,ncol = 2) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF(PSP/Control)),y= '') +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
i

plot_grid(g,i,labels = c('ALZ','PSP'),label_size = 10)

png(file= 'results/figures/Month_human_isoform.png',height = 1500,width = 3500,res = 200)

plot_grid(g,i,labels = c('ALZ','PSP'),label_size = 10)

dev.off()

### Mouse - FAD5X

FAD5X <- readRDS('results/isoform_results/FAD5X_ISW.rds')
FAD5X$Age_Group <- factor(FAD5X$Age_Group,levels= c('4','12','18'))

FAD5X_CCX <- FAD5X[FAD5X$Region == 'Cortex',]
FAD5X_HIP <- FAD5X[FAD5X$Region == 'Hippocampus',]


j = ggplot2::ggplot(FAD5X_CCX, ggplot2::aes(dIF,-log10(iso_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

g = j + ggrepel::geom_text_repel(data=FAD5X_CCX,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-0.05,0.05),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(legend.position = 'none',
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF(FAD5X/Control)),y= expression(-log[10](FDR))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
g

h = ggplot2::ggplot(FAD5X_HIP, ggplot2::aes(dIF,-log10(iso_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

i = h + ggrepel::geom_text_repel(data=FAD5X_HIP,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-0.05,0.05),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 15),
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Age_Group,ncol = 3) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF(PSP/Control)),y= '') +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
i

plot_grid(g,i,labels = c('Cortex','Hippocampus'),label_size = 10)


png(file= 'results/figures/Month_mouse_isoform.png',height = 1500,width = 3500,res = 200)

plot_grid(g,i,labels = c('Cortex','Hippocampus'),label_size = 10)

dev.off()


### Mouse - Tau

Tau <- readRDS('results/isoform_results/Tau_ISW.rds')
Tau$Age_Group <- factor(Tau$Age_Group,levels = c('4','17'))


j = ggplot2::ggplot(Tau, ggplot2::aes(dIF,-log10(iso_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

g = j + ggrepel::geom_text_repel(data=Tau,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("")  + theme_classic2()+
  geom_vline(xintercept =c(-0.05,0.05),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.01),linetype='dashed',color='red')+ 
  theme(legend.position = 'right',
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~Age_Group,ncol = 2) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF(TauD35/Control)),y= expression(-log[10](FDR))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
g

png(file= 'results/figures/Month_Tau_isoform.png',height = 1500,width = 3500,res = 200)

g

dev.off()
