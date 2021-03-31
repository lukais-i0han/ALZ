#### Plot of genes 
library('ggplot2')
library('stringr')
library('ggrepel')
library('cowplot')

### Import of table

FAD5X <- read.csv('archives/FAD5X_all.csv',row.names = 1,stringsAsFactors = F,
            header = T)

theme_classic2 <- function(base_size = 15, base_family = ""){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border     = element_blank(),
      panel.background = element_rect(color = 'black',size = 0.5),
      axis.line        = element_line(colour = "black"),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.background = element_blank(),
      legend.key       = element_blank()
      
    )
}

FAD5X_CCX <- FAD5X[FAD5X$Area == 'Cortex',]
FAD5X_CCX <- FAD5X_CCX[FAD5X_CCX$Gene_Symbol %in% c('Stat3','Tcf4','Ptk2b'),]

FAD5X_HIP <- FAD5X[FAD5X$Area == 'Hippocampus',]
FAD5X_HIP <- FAD5X_HIP[FAD5X_HIP$Gene_Symbol %in% c('Stat3','Tcf4','Ptk2b'),]


j = ggplot2::ggplot(FAD5X_CCX, 
  ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

g = j +geom_text_repel(data=FAD5X_CCX,aes(label=Gene_Symbol),nudge_x = 0.3,nudge_y = 0.3)+
  ggtitle("")  + theme_minimal()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='solid')+
  geom_hline(yintercept = -log10(0.01),linetype='solid',color='red')+ 
  theme(legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 15),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right') +ylim(0,30)+ xlim(c(-1.5,1.5))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](padj))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
g<-g+facet_wrap(~Month,ncol = 3)

h = ggplot2::ggplot(FAD5X_HIP, 
                    ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

i = h +geom_text_repel(data=FAD5X_HIP, aes(label=Gene_Symbol),
                  nudge_x = 0.3,nudge_y = 0.3)+
  ggtitle("")  + theme_minimal()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='solid')+
  geom_hline(yintercept = -log10(0.01),linetype='solid',color='red')+ 
  theme(legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 15),
        axis.text.x=element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right') +ylim(0,30)+ xlim(c(-1.5,1.5))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](padj))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
i<-i+facet_wrap(~Month,ncol = 3)

plot_grid(g,i,nrow = 2,labels = c('Cortex','Hippocampus'))

png(file= 'results/figures/FAD5X',height = 1500,width = 3500,res = 200)

plot_grid(g,i,nrow = 2,labels = c('Cortex','Hippocampus'))

dev.off()

#### Tau

Taud35 <- read.csv('archives/Tau_all.csv',header = T,
                   stringsAsFactors = F, row.names = 1)

Taud35 <- Taud35[Taud35$Gene_Name %in% c('Stat3','Tcf4','Ptk2b'),]

j = ggplot2::ggplot(Taud35, 
                    ggplot2::aes(log2FoldChange,-log10(padj))) +
  ggplot2::geom_point(ggplot2::aes(color = DEG),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('red','black')) 

l = j +geom_text_repel(data=Taud35, aes(label=Gene_Name),
                       nudge_x = 0.3,nudge_y = 0.3)+
  ggtitle("TauD35")  + theme_minimal()+
  geom_vline(xintercept =c(-log2(1.3),log2(1.3)),color = "red",linetype='solid')+
  geom_hline(yintercept = -log10(0.01),linetype='solid',color='red')+ 
  theme(legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 15),
        axis.text.x=element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right') +ylim(0,30)+ xlim(c(-1.5,1.5))+
  labs(x= expression(log[2](FC)),y= expression(-log[10](padj))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
l<-l+facet_wrap(~Month,ncol = 2)

png(file= 'results/figures/Taud35',height = 1500,width = 3500,res = 200)

l

dev.off()



