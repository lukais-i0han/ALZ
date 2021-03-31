#### Terms and DEG analysis

library('ggplot2')
library('dplyr')
library('stringr')
library('RedeR')

theme_classic2 <- function(base_size = 15, base_family = ""){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border     = element_blank(),
      panel.background = element_rect(color = 'black',size = 0.5),
      axis.line        = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.minor.x = element_line(element_blank()),
      strip.background = element_blank(),
      legend.key       = element_blank(),
      panel.grid.major.x = element_blank()
      
    )
}

### Intercept Human-Mouse

ALZ_PSP_DEG <- readRDS('trans_results/ALZ_PSP_DEG.rds')
ALZ_PSP_Terms <- readRDS('trans_results/ALZ_PSP_Terms.rds')
colnames(ALZ_PSP_Terms)[colnames(ALZ_PSP_Terms) == 'Group'] <- 'Region'
colnames(ALZ_PSP_Terms)[colnames(ALZ_PSP_Terms) == 'Age_Group'] <- 'Group'

CCX_HIP_Terms <- readRDS('trans_results/CCX_HIP_Terms.rds')
colnames(CCX_HIP_Terms)[colnames(CCX_HIP_Terms) == 'Month_Group'] <- 'Group'

MAYO_FAD_terms <- rbind(ALZ_PSP_Terms,CCX_HIP_Terms)

MAYO_FAD_terms_neural <- MAYO_FAD_terms[str_detect(MAYO_FAD_terms$Term,'neur|syna'),]

MAYO_FAD_terms_neural <- MAYO_FAD_terms_neural %>% group_by(Group) %>%
  slice_max(pAdj,n = 5)

MAYO_FAD_terms_neural <- MAYO_FAD_terms_neural[order(MAYO_FAD_terms_neural$Term,
                        decreasing = T),]

MAYO_FAD_terms_neural$Region <- factor(MAYO_FAD_terms_neural$Region,
                    levels = c('ALZ','PSP','Cortex','Hippocampus'))

MAYO_FAD_terms_neural$Group <- factor(MAYO_FAD_terms_neural$Group,
      levels = c('A_B_C','A_B','A','B','C','4_12_18','4','12','18'))

ggplot(MAYO_FAD_terms_neural,aes(x=Group, y=factor(Term, levels = unique(Term)), fill=-log10(as.numeric(pAdj)))) +
  geom_tile(size = 0.1, colour = "black") + scale_fill_gradient2(low="#ff8ca5", high="#f71d4d", na.value = "white") +# coord_equal() +
  labs(y=NULL, x=NULL, fill=expression(-Log[10](FDR))) + facet_grid(~Region, scales = "free", space = "free") +
  theme_classic2() + theme(
    strip.text.x = element_text(size = 10), plot.title = element_text(size = 20),
    strip.text.y = element_text(size = 10, color = "black", face = "bold"))

g <- ggplot(MAYO_FAD_terms_neural,aes(x=Group, y=factor(Term, levels = unique(Term)), fill=-log10(as.numeric(pAdj)))) +
  geom_tile(size = 0.1, colour = "black") + scale_fill_gradient2(low="#ff8ca5", high="#f71d4d", na.value = "white") +# coord_equal() +
  labs(y=NULL, x=NULL, fill=expression(-Log[10](FDR))) + facet_grid(~Region, scales = "free", space = "free") +
  theme_classic2() + theme(
    strip.text.x = element_text(size = 10), plot.title = element_text(size = 20),
    strip.text.y = element_text(size = 10, color = "black", face = "bold"))

png(file= 'trans_results/MAYO_5XFAD_neural.png',height = 1600,res = 300,
    width = 3500)

g

dev.off()

saveRDS(MAYO_FAD_terms_neural,'trans_results/MAYO_FAD_neural.rds')
### cluster visualization

