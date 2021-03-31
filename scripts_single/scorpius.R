### Trajectories with Seurat and SCORPIUS

library(SCORPIUS)
library(Seurat)

marie_sub <- readRDS('marie_data/results/seurat-combined.rds')

counts <- marie_sub@assays[["RNA"]]@counts


sample_info <- data.frame(clusters=marie_sub@active.ident)


srt <- CreateSeuratObject(counts = counts, 
                meta.data = sample_info)

expression <- t(as.matrix(srt@assays$RNA@data))
group_name <- srt@meta.data$clusters

space <- reduce_dimensionality(expression, 
        dist = 'spearman',ndim = 3)
draw_trajectory_plot(space, progression_group = group_name, contour = TRUE)
traj <- infer_trajectory(space)

draw_trajectory_plot(
  space, 
  progression_group = group_name,
  path = traj$path,
  contour = T,path_size = 1
)

gimp <- gene_importances(
  expression,
  num_permutations = 10,  
  traj$time) 

gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
gene_sel <- gimp$gene[gimp$qvalue < 0.05]
expr_sel <- scale_quantile(expression[,gene_sel])

draw_trajectory_heatmap(expr_sel, traj$time, group_name,
                        show_labels_row = F)

modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
draw_trajectory_heatmap(expr_sel, traj$time,show_labels_row = F 
                        ,group_name, modules)
