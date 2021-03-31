#### Integration pipeline

marie_tsne <- readRDS('marie_data/results/seurat_tsne.rds')

dat <- readRDS('marie_data/results/marie_data_clean.rds')

metadata <- read.csv('marie_data/results/metadata.csv',header = T,
                     stringsAsFactors = F, row.names = 1)

row_marie <- data.frame(cell_type=marie_tsne@active.ident,stringsAsFactors = F)
row_marie$'Cells' <- rownames(row_marie)
row_marie <- row_marie[row_marie$cell_type != 'Microglia',]
idx <- str_detect(row_marie$Cells,'Ascl1_SC')
row_marie <- row_marie[!idx,]

dat <- dat[,colnames(dat) %in% row_marie$Cells,]
metadata <- metadata[rownames(metadata) %in% colnames(dat),]


marie_sub <- CreateSeuratObject(counts = dat,min.cells = 3,min.features = 200,
                                project = 'Hipp',meta.data = metadata,assay = 'RNA')

marie_sub[["percent.mt"]] <- PercentageFeatureSet(marie_sub, pattern = "^MT-")

marie_sub <- subset(marie_sub,nFeature_RNA > 300  & percent.mt < 15)

marie_list <- SplitObject(marie_sub, split.by = "mouse")

marie_list <- lapply(X = marie_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = marie_list)
marie_list <- lapply(X = marie_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE,approx=FALSE)
})

marie_anchors <- FindIntegrationAnchors(object.list = marie_list, 
      anchor.features = features, reduction = "rpca",k.anchor = 20)

marie_combined <- IntegrateData(anchorset = marie_anchors,k.weight = 20)

DefaultAssay(marie_combined) <- "integrated"

marie_combined <- ScaleData(marie_combined, verbose = FALSE)
marie_combined <- RunPCA(marie_combined, npcs = 30, verbose = FALSE)
marie_combined <- RunUMAP(marie_combined, reduction = "pca", dims = 1:30)
marie_combined <- FindNeighbors(marie_combined, reduction = "pca", dims = 1:30)
marie_combined <- FindClusters(marie_combined, resolution = 0.7)

p1 <- DimPlot(marie_combined, reduction = "umap", group.by = "mouse")
p2 <- DimPlot(marie_combined, reduction = "umap", label = TRUE, 
              repel = TRUE)
p1+p2
p2

markers.to.plot <- c('C1qa','Ccnb1','Dcx','Gad1','Gfap','Olig1')


FeaturePlot(marie_combined,features = markers.to.plot,label = T)   
DimPlot(marie_combined,label = T,split.by  = "cell_type")


saveRDS(marie_combined,'marie_data/results/seurat-combined.rds')




