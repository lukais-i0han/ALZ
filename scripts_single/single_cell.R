### raw data
library('stringr')
library('dplyr')
library('Seurat')
library('patchwork')
library('ggplot2')


# Turn off scientific notation
options(scipen=999)

## Cleaning of the data

table_1 <- read.csv('marie_data/raw-data-2019.csv',header = T, stringsAsFactors = F)
table_1 <- table_1[str_detect(table_1$geneID,'ENSMUSG'),]
table_1 <- table_1[!duplicated(table_1$geneID),]
rownames(table_1) <- table_1$geneID
table_1 <- table_1[,-1]
table_1 <- table_1[,!(str_detect(colnames(table_1),'bulk'))]
table_1$'X' <- rownames(table_1)


table_2 <- read.csv('marie_data/2020-08-14_rawdata.csv',header = T, stringsAsFactors = F)
table_2 <- table_2[str_detect(table_2$geneID,'ENSMUSG'),]
table_2 <- table_2[,1:179]
rownames(table_2) <- table_2$geneID

table_2_genes <- table_2[,c(1,178:179)]

table_2 <- table_2[,-c(1,178:179)]
table_2$'X' <- rownames(table_2)


dat <- full_join(table_1,table_2,by= "X")
rownames_dat <- dat$X
dat <- dat[,colnames(dat) != 'X']
rownames(dat) <- rownames_dat
dat[is.na(dat)] <-0

dat <- sapply(dat,as.numeric)
rownames(dat) <- rownames_dat
dat <- as.data.frame(dat)

table_2_genes <- table_2_genes[!duplicated(table_2_genes$symbol),]
table_2_genes <- table_2_genes[complete.cases(table_2_genes[,1:2]),]

dat <- dat[rownames(dat) %in% table_2_genes$geneID,]

mito_genes <- read.csv('marie_data/mito_genes.csv',header = T, stringsAsFactors = F)
mito_genes <- mito_genes[mito_genes$Tissues != '',]
mito_genes <- mito_genes[str_detect(mito_genes$Tissues,
                                    "all 14|cerebrum|cerebellum"),]


table_2_genes$'row_mito' <- table_2_genes$geneID %in% mito_genes$EnsemblGeneID
table_2_genes$'MT-genes' <-  ifelse(table_2_genes$row_mito == TRUE,
                                    paste('MT',table_2_genes$symbol,sep = '-'),
                                    table_2_genes$symbol)
table_2_genes <- table_2_genes[,c(1,2,3,5)]
table_2_genes <- table_2_genes[!duplicated(table_2_genes$`MT-genes`),]

write.csv(table_2_genes,'marie_data/results/gene_table.csv')

dat <- dat[order(rownames(dat)),]
table_2_genes <- table_2_genes[order(table_2_genes$geneID),]

dat <- dat[rownames(dat) %in% table_2_genes$geneID,]

rownames(dat) <- table_2_genes$`MT-genes`


saveRDS(dat,'marie_data/results/marie_data_clean.rds')


#### Single cell analysis

dat <- readRDS('marie_data/results/marie_data_clean.rds')

mouse1 <- rep('mouse1',181)
mouse2 <- rep('mouse2',176)


cell_ann <- data.frame(mouse=c(mouse1,mouse2),
                       cell_type = colnames(dat))

cell_ann <- cell_ann %>% mutate(cell_type = case_when(
  str_detect(cell_type,'Ascl1_Dlx2') ~ 'Ascl1_Dlx2',
  str_detect(cell_type,'Ascl1_S') ~ 'Ascl1',
  str_detect(cell_type,'Control') ~ 'Control',
  TRUE ~ 'Dlx2'
  
))
rownames(cell_ann) <- colnames(dat)

write.csv(cell_ann,'marie_data/results/metadata.csv')

### Construction of the Seurat Object
marie_hip <- CreateSeuratObject(counts = dat,min.cells = 3,min.features = 200,
              project = 'Hipp',meta.data = cell_ann,assay = 'RNA')

marie_hip_SC <- Seurat::as.SingleCellExperiment(marie_hip)

marie_hip[["percent.mt"]] <- PercentageFeatureSet(marie_hip, pattern = "^MT-")


plot1 <- FeatureScatter(marie_hip, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(marie_hip, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 +plot2

### subseting

marie_hip <- subset(marie_hip,nFeature_RNA > 300  & percent.mt < 15)

marie_hip <- NormalizeData(marie_hip, normalization.method = "LogNormalize", 
                           scale.factor = 10000)

marie_hip <- FindVariableFeatures(marie_hip, 
                                  selection.method = "vst", nfeatures = 2000)               

top10 <- head(VariableFeatures(marie_hip), 10)

plot1 <- VariableFeaturePlot(marie_hip)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(marie_hip)
marie_hip <- ScaleData(marie_hip, features = all.genes)

##PCA

marie_hip <- RunPCA(marie_hip, features = VariableFeatures(object = marie_hip))

print(marie_hip[["pca"]], dims = 1:5, nfeatures = 5)

DimHeatmap(marie_hip, dims = 1, cells = 10, balanced = TRUE)

ElbowPlot(marie_hip)

marie_hip <- FindNeighbors(marie_hip, dims = 1:15)

### CLuster the cells

marie_fc <- FindClusters(
  marie_hip,
  dims.use = 1:15,
  force.recalc = TRUE,
  print.output = TRUE,
  resolution = 0.6,
  save.SNN = TRUE)


marie_tsne <- RunTSNE(marie_fc,dims = 1:15)

DimPlot(marie_tsne,
        label.size = 5,
        reduction = 'tsne',label = T) + labs(title = 'Clustering of 317 cells')


DimPlot(marie_tsne,
        label.size = 3,split.by = 'cell_type',
        reduction = 'tsne',label = T) + labs(title = "FAC's markers")

DimPlot(marie_tsne,
        label.size = 5,split.by = 'mouse',
        reduction = 'tsne',label = T)+labs(title = "Clusters by mouse")

brain_markers <- c('Gfap','Olig1','C1qa','Cd68','Gad1','Pvalb',
                   'Sox1','Dcx','Map2')  

VlnPlot(tsne, features = brain_markers,group.by = 'cell_type',
        slot='counts')

FeaturePlot(marie_tsne,features = 
              c('Gfap','Olig1','C1qa','Gad1',
                'Dcx','Ccnb1'),label = T)          


new_cluster_ids <- c('Mix_Cells','Microglia','Induced_GABA',
                     'Astrocyte','Microglia','Astrocyte',
                     'Neuroblast')

names(new_cluster_ids) <- levels(marie_tsne)
marie_tsne <- RenameIdents(marie_tsne,new_cluster_ids)

DimPlot(marie_tsne, reduction = "tsne", 
        label = TRUE, pt.size = 1) 

saveRDS(marie_tsne,'marie_data/results/seurat_tsne.rds')

####

marie_tsne <- readRDS('marie_data/results/seurat_tsne.rds')

levels(marie_tsne) <- c('Astrocyte','Induced_GABA','Microglia','Mix_Cells',
                  'Neuroblast')

marie_hip_markers <- FindAllMarkers(marie_tsne, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)


marie_hip_10  <- marie_hip_markers  %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


DotPlot(marie_tsne, features = c('Gfap','Olig1','C1qa','Gad1',
                                 'Dcx','Ccnb1'), 
        cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()



png(file= 'marie_data/results/figures/dot_319.png',height = 1500,width = 3500,res = 200)

DoHeatmap(marie_tsne, features = marie_hip_10$gene,
          angle = 45,size = 3.5) + NoLegend()

dev.off()

DimPlot(marie_tsne,
        label.size = 4,
        reduction = 'tsne',label = T,sizes.highlight = T) + labs(title = 'Clustering of 317 cells')

png(file= 'marie_data/results/figures/FACs_317.png',
    height = 1500,width = 3500,res = 200)

DimPlot(marie_tsne,
        label.size = 3,split.by = 'cell_type',
        reduction = 'tsne',label = T) + labs(title = "FAC's markers")

dev.off()

DimPlot(marie_tsne,
        label.size = 5,split.by = 'mouse',
        reduction = 'tsne',label = T)+labs(title = "Clusters by mouse")
