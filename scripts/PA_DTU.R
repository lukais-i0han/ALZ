### PA DTU

### DTU analysis
library('IsoformSwitchAnalyzeR')
library('dplyr')
library('stringr')
## Human


mayo_A_ISW <- readRDS('archives/ISW/MAYO/Age_A.fullAnalysis.rds')
mayo_B_ISW <- readRDS('archives/ISW/MAYO/Age_B.fullAnalysis.rds')
mayo_C_ISW <- readRDS('archives/ISW/MAYO/Age_C.fullAnalysis.rds')

## Human
mayo_A_features <- mayo_A_ISW$isoformFeatures

mayo_A_analysis <- mayo_A_ISW$isoformSwitchAnalysis
mayo_A_analysis <- mayo_A_analysis[order(mayo_A_analysis$iso_ref),]

mayo_A_features$gene_q_value <- mayo_A_analysis$gene_q_value
mayo_A_features$iso_q_value <- mayo_A_analysis$iso_q_value
mayo_A <- mayo_A_features[,c(3:9,16,22:23,27:28)]
mayo_A$'AgeAtDeath' <- 'A' 

mayo_B_features <- mayo_B_ISW$isoformFeatures

mayo_B_analysis <- mayo_B_ISW$isoformSwitchAnalysis
mayo_B_analysis <- mayo_B_analysis[order(mayo_B_analysis$iso_ref),]

mayo_B_features$gene_q_value <- mayo_B_analysis$gene_q_value
mayo_B_features$iso_q_value <- mayo_B_analysis$iso_q_value
mayo_B <- mayo_B_features[,c(3:9,16,22:23,27:28)]
mayo_B$'AgeAtDeath' <- 'B'

mayo_C_features <- mayo_C_ISW$isoformFeatures

mayo_C_analysis <- mayo_C_ISW$isoformSwitchAnalysis
mayo_C_analysis <- mayo_C_analysis[order(mayo_C_analysis$iso_ref),]

mayo_C_features$gene_q_value <- mayo_C_analysis$gene_q_value
mayo_C_features$iso_q_value <- mayo_C_analysis$iso_q_value
mayo_C <- mayo_C_features[,c(3:9,16,22:23,27:28)]
mayo_C$'AgeAtDeath' <- 'C'

mayo_ISW <- rbind(mayo_A,mayo_B,mayo_C)

mayo_ISW <- mayo_ISW[mayo_ISW$condition_1 == 'Control' &
                    mayo_ISW$condition_2 == 'Pathologic_Aging',]

mayo_ISW <- mayo_ISW  %>% mutate(DTU = case_when(
  abs(mayo_ISW$dIF) > 0.1 & mayo_ISW$iso_q_value < 0.05 ~ 'DTU',
  TRUE ~ 'Non-DTU'
))
mayo_ISW <- mayo_ISW %>% group_by(gene_name) %>% 
  filter(n()>1)

mayo_ISW_sig <- mayo_ISW[mayo_ISW$DTU == 'DTU',]
mayo_ISW_sig <- mayo_ISW_sig[order(mayo_ISW_sig$AgeAtDeath),]

write.csv(mayo_ISW_sig,'marcos_arquivos/DTU_PA.csv')

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


j = ggplot2::ggplot(mayo_ISW, ggplot2::aes(dIF,-log10(iso_q_value))) +
  ggplot2::geom_point(ggplot2::aes(color = DTU),size=2) +
  ggplot2::scale_shape_manual(values = c(20,20))+
  ggplot2::scale_color_manual(values = c('blue','black')) 

g = j + ggrepel::geom_text_repel(data=mayo_ISW,point.padding = 20)+
  ggplot2::aes(label="")+
  ggtitle("Pathologic Aging")  + theme_classic2()+
  geom_vline(xintercept =c(-0.1,0.1),color = "red",linetype='dashed')+
  geom_hline(yintercept = -log10(0.05),linetype='dashed',color='red')+ 
  theme(
        axis.text.x = element_text(face = "bold", color = "black",size = 15),
        axis.title = element_text(face='bold',color='black',size=15),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 15)
        
  )  +
  theme(legend.justification = 'top',
        legend.box.just = 'right')+
  facet_wrap(~AgeAtDeath,ncol = 3) + ylim(0,20)+ xlim(c(-0.35,0.35))+
  labs(x= expression(dIF),y= expression(-log[10](iso_q_value))) +
  theme(
    strip.text.x = element_text(
      size = 15, color = "black", face = "bold"
    ),
    strip.text.y = element_text(
      size = 15, color = "black", face = "bold"
    )
  )
g
