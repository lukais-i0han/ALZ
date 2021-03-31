#### Slingshot with subset of Seurat Data

library('Seurat')
library('stringr')
marie_tsne <- readRDS('marie_data/results/seurat_tsne.rds')



### Import of  data to subset
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


### Construct of Seurat object
marie_sub <- CreateSeuratObject(counts = dat,min.cells = 3,min.features = 200,
            project = 'Hipp',meta.data = metadata,assay = 'RNA')

saveRDS(marie_sub,'marie_data/results/marie_no_ASCL1.rds')

