#### Overlap between Human and Mouse
library('CEMiTool')
library('stringr')
library('ggplot2')
library('fgsea')

### DEGs

ALZ_DEG <- read.csv('archives/ALZ_sig.csv',header = T, stringsAsFactors = F,row.names = 1)
PSP_DEG <- read.csv('archives/PSP_sig.csv',header = T, stringsAsFactors = F, row.names = 1)
FAD5X_DEG <- read.csv('archives/FAD5X_sig.csv',header = T, stringsAsFactors = F, row.names = 1)
TAU_DEG <- read.csv('archives/Tau_sig.csv',header = T, stringsAsFactors = F, row.names = 1)

### DTUs

MAYO_DTU <- readRDS('results/isoform_results/mayo_ISW_sig.rds')
FAD5X_DTU <- readRDS('results/isoform_results/FAD5X_ISW_sig.rds')
TAU_DTU <- readRDS('results/isoform_results/Tau_ISW_sig.rds')

Human_genes <- c(ALZ_DEG$Gene_Symbol,PSP_DEG$Gene_Symbol,MAYO_DTU$gene_name)


Mouse_genes <- c(FAD5X_DEG$Gene_Symbol,TAU_DEG$Gene_Name,FAD5X_DTU$gene_name,
               TAU_DTU$gene_name)



### normalized matrix

vsd_mayo <- readRDS('CEMtool/vsd_files/vsd_mayo.rds')
vsd_FAD <- readRDS('CEMtool/vsd_files/vsd_FAD.rds')
vsd_TAU <- readRDS('CEMtool/vsd_files/vsd_TAU.rds')

#### metadado

mayo_metadado <- read.csv('CEMtool/metadado/mayo_metadado.csv',header = T,stringsAsFactors = F,
                        row.names = 1)

FAD_metadado <- read.csv('CEMtool/metadado/FAD_metadado.csv',header = T, 
                         stringsAsFactors = F,row.names = 1)

TAU_metadado <- read.csv('CEMtool/metadado/Tau_metadado.csv',header = T,
                         stringsAsFactors = F, row.names = 1)


### Human

vsd_mayo <- assay(vsd_mayo)
vsd_mayo <- as.data.frame(vsd_mayo)

mayo_meta <- data.frame(SampleName = mayo_metadado$ID,
                        Class = mayo_metadado$Group)

mayo_meta <- mayo_meta %>% mutate(Class = case_when(
  Class == 'AD' ~ 'ALZ',
  Class == 'Control' ~ 'CON',
  TRUE ~ 'PSP'
))

mayo_meta$Class <- factor(mayo_meta$Class,levels = c('ALZ','PSP','CON'))


### Human Pathway and network

human_GO<-read_gmt('CEMtool/gmt_files/Human_GO_AllPathways_no_GO_iea_February_05_2021_symbol.gmt')

int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)

mouse_GO <- read_gmt('CEMtool/gmt_files/Mouse_GO_AllPathways_no_GO_iea_February_05_2021_symbol.gmt')
dump_human <- str_detect(mouse_GO$term,'HUMAN')
mouse_GO <- mouse_GO[!dump_human,]

### Create cem archive


Human_cem <- cemitool(expr = vsd_mayo,annot = mayo_meta,
      interactions = int_df,filter_pval = 0.05,ora_pval = 0.01,plot_diagnostics = T)



### Module enrichment

Human_cem <- mod_gsea(Human_cem)
Human_cem <- plot_gsea(Human_cem)


### Expression patterns in modules

Human_cem <- plot_profile(Human_cem)


####  Over Representation Analysis
Human_cem <- mod_ora(Human_cem, human_GO)
Human_cem<- plot_ora(Human_cem)


### Interactions

int_df <- int_df[int_df$gene1symbol %in% Human_genes &
                 int_df$gene2symbol %in% Human_genes,]

interactions_data(Human_cem) <- int_df
Human_cem <- plot_interactions(Human_cem)

write_files(Human_cem, directory="CEMtool/human/Tables",force=T)

### Mouse FAD

vsd_FAD <- assay(vsd_FAD)
vsd_FAD <- as.data.frame(vsd_FAD)

FAD_meta <- data.frame(SampleName = FAD_metadado$SpecimenID,
                           Class = paste(FAD_metadado$Region,FAD_metadado$Group,sep='-'))


### Create cem archive

Mouse_cem <- cemitool(vsd_FAD,FAD_meta,force_beta = T,
            filter = T, filter_pval = 0.05,ora_pval = 0.01)

### Module enrichment

Mouse_cem <- mod_gsea(Mouse_cem)
Mouse_cem <- plot_gsea(Mouse_cem)


### Expression patterns in modules

Mouse_cem <- plot_profile(Mouse_cem)


####  Over Representation Analysis

Mouse_cem <- mod_ora(Mouse_cem,mouse_GO)
Mouse_cem <- plot_ora(Mouse_cem,"ora")


### Interactions

mouse_PPI <- readRDS('CEMtool/PPI/mouse_PPI.rds')

interactions_data(Mouse_cem) <- mouse_PPI
Mouse_cem <- plot_interactions(Mouse_cem)


write_files(Mouse_cem, directory="CEMtool/Tables",force=T)


### TAU
### Using the same PPI and mouse_GO 

vsd_TAU <- assay(vsd_TAU)
vsd_TAU <- as.data.frame(vsd_TAU)

TAU_meta <- data.frame(SampleName = TAU_metadado$id,
                       Class = TAU_metadado$mutation)


### Create cem archive

TAU_cem <- cemitool(vsd_TAU,TAU_meta,filter = T, filter_pval = 0.05,
                    ora_pval = 0.01)

### Module enrichment

TAU_cem <- mod_gsea(TAU_cem)
TAU_cem <- plot_gsea(TAU_cem)


### Expression patterns in modules

TAU_cem<- plot_profile(TAU_cem)


####  Over Representation Analysis

TAU_cem <- mod_ora(TAU_cem,mouse_GO)
TAU_cem<- plot_ora(TAU_cem,"ora")


### Interactions

mouse_PPI <- readRDS('CEMtool/PPI/mouse_PPI.rds')

interactions_data(TAU_cem) <- mouse_PPI
TAU_cem <- plot_interactions(TAU_cem)

write_files(TAU_cem, directory="CEMtool/Tables",force=T)


### Save all
saveRDS(Human_cem,'CEMtool/human_cem.rds')
saveRDS(Mouse_cem,'CEMtool/FAD_cem.rds')
saveRDS(TAU_cem,'CEMtool/Tau_cem.rds')

### Graphics 

Human_cem <- readRDS('CEMtool/results/human_cem.rds')
FAD_cem <- readRDS('CEMtool/results/FAD_cem.rds')
TAU_cem <- readRDS('CEMtool/results/Tau_cem.rds')

gsa_human <- show_plot(Human_cem,'gsea')
gsa_FAD <- show_plot(FAD_cem,'gsea')
gsa_TAU <- show_plot(TAU_cem,'gsea')

patternts_Human <- show_plot(Human_cem, "profile")
patternts_FAD <- show_plot(FAD_cem, "profile")
patternts_TAU <- show_plot(TAU_cem, "profile")

for(i in 1:length(patternts_TAU)){
  png(file=paste0('CEMtool//results/mouse_tau/module_',i,'.png'),height = 1500,width = 3500,res = 200)
  print(patternts_TAU[[i]])
  dev.off()
}
