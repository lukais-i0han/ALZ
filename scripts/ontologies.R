### 
library('dplyr')
library('stringr')
library('RColorBrewer')
library('tidyr')

ALZ_DEG <- read.csv('marcos_arquivos/DEG_ALZ.csv',header = T, stringsAsFactors = F,
                    row.names = 1)

ALZ_DEG <- ALZ_DEG %>% mutate(UP_DOWN = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))

ALZ_DTU <- read.csv('marcos_arquivos/DTU_ALZ.csv',header = T, stringsAsFactors = F,
                    row.names = 1)

PA_DEG <- read.csv('marcos_arquivos/DEG_PA.csv',header = T, stringsAsFactors = F,
                   row.names = 1)
PA_DTU <- read.csv('marcos_arquivos/DTU_PA.csv',header = T, stringsAsFactors = F,
                   row.names = 1)

PSP_DEG <- read.csv('marcos_arquivos/DEG_PSP.csv',header = T, stringsAsFactors = F,
                    row.names = 1)

PSP_DEG <- PSP_DEG %>% mutate(UP_DOWN = case_when(
  log2FoldChange > log2(1.3) ~ 'Up',
  TRUE ~ 'Down'
))


PSP_DTU <- read.csv('marcos_arquivos/DTU_PSP.csv',header = T, stringsAsFactors = F,
                    row.names = 1)

#### Ontologies

library('gprofiler2')

## DEG
DEG <- gost(query = list(
 'ALZ' = ALZ_DEG$Gene_Symbol, 
 'PA' = PA_DEG$Gene_Symbol,
 'PSP' = PSP_DEG$Gene_Symbol),
organism = 'hsapiens',
significant = T,
multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
user_threshold = 0.01,
correction_method = 'fdr')
DEG_fprofiler <- DEG$result
DEG_fprofiler$'Query' <- 'DEG'


DTU  <- gost(query = list(
  'ALZ' = ALZ_DTU$gene_name, 
  'PA' = PA_DTU$gene_name,
  'PSP' = PSP_DTU$gene_name),
  organism = 'hsapiens',
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')
Dtu_fprofiler <- DTU$result
Dtu_fprofiler$'Query' <- 'DTU'

DGE_DTU <-gost(query = list(
  'ALZ' = c(ALZ_DEG$Gene_Symbol,ALZ_DTU$gene_name), 
  'PA' = c(PA_DEG$Gene_Symbol,PA_DTU$gene_name),
  'PSP' = c(PSP_DEG$Gene_Symbol,PSP_DTU$gene_name)),
  organism = 'hsapiens',
  significant = T,
  multi_query = F,sources = c('GO:BP','GO:MF','GO:CC'),
  user_threshold = 0.01,
  correction_method = 'fdr')
DGE_DTUfprofiler <- DGE_DTU$result
DGE_DTUfprofiler$'Query' <- 'DEG+DTU'


DEG_DTU_DEGDTU <- rbind(DEG_fprofiler,Dtu_fprofiler,DGE_DTUfprofiler)

df <- apply(DEG_DTU_DEGDTU,2,as.character)
write.csv(df,'marcos_arquivos/ontologias.csv')

ALZ <- DEG_DTU_DEGDTU[DEG_DTU_DEGDTU$query == 'ALZ',]
ALZ <- ALZ[str_detect(ALZ$term_name,'neur|syna'),]


PA <- DEG_DTU_DEGDTU[DEG_DTU_DEGDTU$query == 'PA',]
PA <- PA[str_detect(PA$term_name,'neur|syna'),]

PSP <- DEG_DTU_DEGDTU[DEG_DTU_DEGDTU$query == 'PSP',]
PSP <- PSP[str_detect(PSP$term_name,'neur|syna'),]

ALZ_PSP <- rbind(ALZ,PA,PSP)

ALZ_PSP_PA$Query <- factor(ALZ_PSP_PA$Query,levels =
            c('DEG','DEG+DTU','DTU'))

ALZ_PSP_PA <- ALZ_PSP_PA[ALZ_PSP_PA$source !='GO:MF',]

ALZ_PSP_PA_BP <- ALZ_PSP_PA[ALZ_PSP_PA$source == 'GO:BP',]
ALZ_PSP_PA_BP <- ALZ_PSP_PA_BP[order(ALZ_PSP_PA_BP$term_name,decreasing = T),]

ALZ_PSP_PA_CC <- ALZ_PSP_PA[ALZ_PSP_PA$source == 'GO:CC',]
ALZ_PSP_PA_CC <- ALZ_PSP_PA_CC[order(ALZ_PSP_PA_CC$term_name,decreasing = T),]


ggplot(ALZ_PSP_PA_BP,aes(x=query, y=factor(term_name, levels = unique(term_name)), fill=Query)) +
  geom_tile(size = 0.1, colour = "black") + 
  scale_fill_manual(values = c('yellow','lightblue','firebrick'),
                   na.value = 'white' ) +# coord_equal() +
  labs(y=NULL, x=NULL) + facet_grid(~source, scales = "free", space = "free") +
  theme_classic(base_size = 9) + theme(axis.text.x=element_text(angle = -45, hjust = 0),
  strip.text.x = element_text(size = 10), plot.title = element_text(size = 20),
  strip.text.y = element_text(size=6),
  axis.title.y = element_text(face='bold')) +
  theme(panel.background = element_rect(fill = 'white'))+
  theme(legend.position = 'right')



