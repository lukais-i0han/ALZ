library('stageR')
library('DRIMSeq')

### ALZ

ALZ_A <- readRDS('stage_archives/refs/Human//Age_A.fullAnalysis.rds')
ALZ_A <- ALZ_A$isoformFeatures
ALZ_A <- ALZ_A[ALZ_A$condition_1 == 'AD' & ALZ_A$condition_2 == 'Control',]

ALZ_B <- readRDS('stage_archives/refs/Human/Age_B.fullAnalysis.rds')
ALZ_B <- ALZ_B$isoformFeatures
ALZ_B <- ALZ_B[ALZ_B$condition_1 == 'AD' & ALZ_B$condition_2 == 'Control',]

ALZ_C <- readRDS('stage_archives/refs/Human/Age_C.fullAnalysis.rds')
ALZ_C <- ALZ_C$isoformFeatures
ALZ_C <- ALZ_C[ALZ_C$condition_1 == 'AD' & ALZ_C$condition_2 == 'Control',]

#####
tx2gene <- read.table('stage_archives/refs/Human//tg2.txt',header = T)

counts_ALZ <- read.csv('stage_archives/refs/Human/counts_ALZ.csv',header = T,
                       stringsAsFactors = F,row.names = 1)

metadado_ALZ <- read.csv('stage_archives/refs/Human/metadado_ALZ.csv',header = T, stringsAsFactors = F,
                         row.names = 1)

colnames(counts_ALZ) <- metadado_ALZ$ID

#### ALZ_A

metadado_ALZ_A <- metadado_ALZ[metadado_ALZ$AgeAtDeath == 'A',]

counts_ALZ_A <- counts_ALZ[,metadado_ALZ_A$ID]
counts_ALZ_A <- counts_ALZ_A[rownames(counts_ALZ_A) %in% ALZ_A$isoform_id,]

tx2gene_A <- tx2gene[tx2gene$transcript %in% rownames(counts_ALZ_A),]

counts_ALZ_A$'feature_id' <- rownames(counts_ALZ_A)
counts_ALZ_A$'gene_id' <- tx2gene_A$ENSG

ALZ_A_samples <- data.frame(sample_id = metadado_ALZ_A$ID,
                  group = metadado_ALZ_A$condition)

d <- dmDSdata(counts = counts_ALZ_A,samples = ALZ_A_samples)

d <- dmFilter(d, min_samps_gene_expr = 34, min_samps_feature_expr = 3,
              min_gene_expr = 10, min_feature_expr = 10)
design_full <- model.matrix(~ condition,data=metadado_ALZ_A)

d <- dmPrecision(d,design = design_full)

d <- dmFit(d, design = design_full, verbose = 1)

d <- dmTest(d, coef = 'conditionControl',verbose = 1)


pScreen <- d@results_gene$pvalue
names(pScreen) <- d@results_gene$gene_id

pConfirmation <-matrix(d@results_feature$pvalue)
rownames(pConfirmation) <- d@results_feature$feature_id

tx2gene_A <- tx2gene_A[tx2gene_A$transcript %in% row_names(pConfirmation),]

colnames(tx2gene_A) <- c('feature_id','gene_id','gene_symbol')

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                      pScreenAdjusted = FALSE, tx2gene = tx2gene_A)

stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                 alpha = 0.01)
getSignificantGenes(stageRObj)
getSignificantTx(stageRObj)

padj_A <- getAdjustedPValues(stageRObj, order = TRUE,
                           onlySignificantGenes = T)
padj_A <- padj

#### ALZ_B

metadado_ALZ_B <- metadado_ALZ[metadado_ALZ$AgeAtDeath == 'B',]

counts_ALZ_B <- counts_ALZ[,metadado_ALZ_B$ID]
counts_ALZ_B <- counts_ALZ_B[rownames(counts_ALZ_B) %in% ALZ_B$isoform_id,]

tx2gene_B <- tx2gene[tx2gene$transcript %in% rownames(counts_ALZ_B),]

counts_ALZ_B$'feature_id' <- rownames(counts_ALZ_B)
counts_ALZ_B$'gene_id' <- tx2gene_B$ENSG

ALZ_B_samples <- data.frame(sample_id = metadado_ALZ_B$ID,
                            group = metadado_ALZ_B$condition)

d <- dmDSdata(counts = counts_ALZ_B,samples = ALZ_B_samples)

d <- dmFilter(d, min_samps_gene_expr = 71, min_samps_feature_expr = 3,
              min_gene_expr = 10, min_feature_expr = 10)
design_full <- model.matrix(~condition,data=metadado_ALZ_B)

d <- dmPrecision(d,design = design_full)

d <- dmFit(d, design = design_full, verbose = 1)

d <- dmTest(d, coef = 'conditionControl',verbose = 1)


pScreen <- d@results_gene$pvalue
names(pScreen) <- d@results_gene$gene_id

pConfirmation <-matrix(d@results_feature$pvalue)
rownames(pConfirmation) <- d@results_feature$feature_id

tx2gene_B <- tx2gene_B[tx2gene_B$transcript %in% row_names(pConfirmation),]

colnames(tx2gene_B) <- c('feature_id','gene_id','gene_symbol')

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                      pScreenAdjusted = FALSE, tx2gene = tx2gene_B)

stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                 alpha = 0.01)
getSignificantGenes(stageRObj)
getSignificantTx(stageRObj)

padj_B <- getAdjustedPValues(stageRObj, order = TRUE,
                           onlySignificantGenes = T)

### ALZ-C

metadado_ALZ_C <- metadado_ALZ[metadado_ALZ$AgeAtDeath == 'C',]

counts_ALZ_C <- counts_ALZ[,metadado_ALZ_C$ID]
counts_ALZ_C <- counts_ALZ_C[rownames(counts_ALZ_C) %in% ALZ_C$isoform_id,]

tx2gene_C <- tx2gene[tx2gene$transcript %in% rownames(counts_ALZ_C),]

counts_ALZ_C$'feature_id' <- rownames(counts_ALZ_C)
counts_ALZ_C$'gene_id' <- tx2gene_C$ENSG

ALZ_C_samples <- data.frame(sample_id = metadado_ALZ_C$ID,
                            group = metadado_ALZ_C$condition)

d <- dmDSdata(counts = counts_ALZ_C,samples = ALZ_C_samples)

d <- dmFilter(d, min_samps_gene_expr = 40, min_samps_feature_expr = 3,
              min_gene_expr = 10, min_feature_expr = 10)
design_full <- model.matrix(~condition,data=metadado_ALZ_C)

d <- dmPrecision(d,design = design_full)

d <- dmFit(d, design = design_full, verbose = 1)

d <- dmTest(d, coef = 'conditionControl',verbose = 1)


pScreen <- d@results_gene$pvalue
names(pScreen) <- d@results_gene$gene_id

pConfirmation <-matrix(d@results_feature$pvalue)
rownames(pConfirmation) <- d@results_feature$feature_id

tx2gene_C <- tx2gene_C[tx2gene_C$transcript %in% rownames(pConfirmation),]

colnames(tx2gene_C) <- c('feature_id','gene_id','gene_symbol')

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                      pScreenAdjusted = FALSE, tx2gene = tx2gene_C)

stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                 alpha = 0.01,allowNA=T)
getSignificantGenes(stageRObj)
getSignificantTx(stageRObj)

padj_C <- getAdjustedPValues(stageRObj, order = TRUE,
                             onlySignificantGenes = T)


write.csv(padj_A,'stage_archives/results/Human/stage_ALZ_A.csv')
write.csv(padj_B,'stage_archives/results/Human/stage_ALZ_B.csv')
write.csv(padj_C,'stage_archives/results/Human/stage_ALZ_C.csv')

#####


tx2gene <- read.table('stage_archives/refs/Human/tg2.txt',header = T)

counts_PSP <- read.csv('stage_archives/refs/Human/counts_PSP.csv',header = T,
                       stringsAsFactors = F,row.names = 1)

metadado_PSP <- read.csv('stage_archives/refs/Human/metadado_PSP.csv',header = T, stringsAsFactors = F,
                         row.names = 1)

colnames(counts_PSP) <- metadado_PSP$ID

### PSP_A

PSP_A <- readRDS('stage_archives/refs/Human/Age_A.fullAnalysis.rds')
PSP_A <- PSP_A$isoformFeatures
PSP_A <- PSP_A[PSP_A$condition_1 == 'Control' & PSP_A$condition_2 == 'PSP',]

metadado_PSP_A <- metadado_PSP[metadado_PSP$AgeAtDeath == 'A',]

counts_PSP_A <- counts_PSP[,metadado_PSP_A$ID]
counts_PSP_A <- counts_PSP_A[rownames(counts_PSP_A) %in% PSP_A$isoform_id,]

tx2gene_A <- tx2gene[tx2gene$transcript %in% rownames(counts_PSP_A),]

counts_PSP_A$'feature_id' <- rownames(counts_PSP_A)
counts_PSP_A$'gene_id' <- tx2gene_A$ENSG



PSP_A_samples <- data.frame(sample_id = metadado_PSP_A$ID,
                            group = factor(metadado_PSP_A$condition,levels = c('PSP','Control'))
                            )

d <- dmDSdata(counts = counts_PSP_A,samples = PSP_A_samples)

d <- dmFilter(d, min_samps_gene_expr = 66, min_samps_feature_expr = 3,
              min_gene_expr = 10, min_feature_expr = 10)


metadado_PSP_A$condition <- factor(metadado_PSP_A$condition,levels = c('PSP','Control'))

design_full <- model.matrix(~ condition,data=metadado_PSP_A)

d <- dmPrecision(d,design = design_full)

d <- dmFit(d, design = design_full, verbose = 1)

d <- dmTest(d, coef = 'conditionControl',verbose = 1)


pScreen <- d@results_gene$pvalue
names(pScreen) <- d@results_gene$gene_id

pConfirmation <-matrix(d@results_feature$pvalue)
rownames(pConfirmation) <- d@results_feature$feature_id

tx2gene_A <- tx2gene_A[tx2gene_A$transcript %in% rownames(pConfirmation),]

colnames(tx2gene_A) <- c('feature_id','gene_id','gene_symbol')

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                      pScreenAdjusted = FALSE, tx2gene = tx2gene_A)

stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                 alpha = 0.01)
getSignificantGenes(stageRObj)
getSignificantTx(stageRObj)

padj_A <- getAdjustedPValues(stageRObj, order = TRUE,
                             onlySignificantGenes = T)

### PSP_B

PSP_B <- readRDS('stage_archives/refs/Human/Age_B.fullAnalysis.rds')
PSP_B <- PSP_B$isoformFeatures
PSP_B <- PSP_B[PSP_B$condition_1 == 'Control' & PSP_B$condition_2 == 'PSP',]

metadado_PSP_B <- metadado_PSP[metadado_PSP$AgeAtDeath == 'B',]

counts_PSP_B <- counts_PSP[,metadado_PSP_B$ID]
counts_PSP_B <- counts_PSP_B[rownames(counts_PSP_B) %in% PSP_B$isoform_id,]

tx2gene_B <- tx2gene[tx2gene$transcript %in% rownames(counts_PSP_B),]

counts_PSP_B$'feature_id' <- rownames(counts_PSP_B)
counts_PSP_B$'gene_id' <- tx2gene_B$ENSG



PSP_B_samples <- data.frame(sample_id = metadado_PSP_B$ID,
                            group = factor(metadado_PSP_B$condition,levels = c('PSP','Control'))
)

d <- dmDSdata(counts = counts_PSP_B,samples = PSP_B_samples)

d <- dmFilter(d, min_samps_gene_expr = 46, min_samps_feature_expr = 3,
              min_gene_expr = 10, min_feature_expr = 10)


metadado_PSP_B$condition <- factor(metadado_PSP_B$condition,levels = c('PSP','Control'))

design_full <- model.matrix(~ condition,data=metadado_PSP_B)

d <- dmPrecision(d,design = design_full)

d <- dmFit(d, design = design_full, verbose = 1)

d <- dmTest(d, coef = 'conditionControl',verbose = 1)


pScreen <- d@results_gene$pvalue
names(pScreen) <- d@results_gene$gene_id

pConfirmation <-matrix(d@results_feature$pvalue)
rownames(pConfirmation) <- d@results_feature$feature_id

tx2gene_B <- tx2gene_B[tx2gene_B$transcript %in% rownames(pConfirmation),]

colnames(tx2gene_B) <- c('feature_id','gene_id','gene_symbol')

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                      pScreenAdjusted = FALSE, tx2gene = tx2gene_B)

stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                 alpha = 0.01)

padj_B <- getAdjustedPValues(stageRObj, order = TRUE,
                             onlySignificantGenes = T)

write.csv(padj_A,'stage_archives/results/Human/stage_PSP_A.csv')
write.csv(padj_B,'stage_archives/results/Human/stage_PSP_B.csv')
