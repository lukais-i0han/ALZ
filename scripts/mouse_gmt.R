### Script to get gmt file to Mouse and network of mouse

library('msigdb')
library('dplyr')
library('stringr')
library('RedeR')
library('igraph')

mouse_GO <- msigdb.genesets(sets=c('C2.CP'),
                            type='symbols',species = 'mouse')
write.gmt(mouse_GO,'CEMtool/gmt_files/mouse.gmt')

human_GO <- msigdb.genesets(sets=c('C5.BP','C5.CC', 'C5.MF'),
                            type='symbols',species = 'human')

write.gmt(human_GO,'CEMtool/gmt_files/human.gmt')

#### Mouse network

association_mouse <- read.table('results/tables/association_mouse.txt',header = F,
                                stringsAsFactors = F)
dictionary_FAD <- read.csv('results/tables/dictionary_FAD.csv',
                           header = T, stringsAsFactors = F, row.names = 1)
dictionary_FAD$ensembl_peptide_id <- gsub('^10090.','',dictionary_FAD$ensembl_peptide_id)

association_mouse$V1 <- gsub('^10090.','',association_mouse$V1)
association_mouse$V2 <- gsub('^10090.','',association_mouse$V2)

association_mouse <- association_mouse[
  association_mouse$V1 %in% dictionary_FAD$ensembl_peptide_id &
  association_mouse$V2 %in% dictionary_FAD$ensembl_peptide_id,]

mouse_network <- igraph::graph.data.frame(d=association_mouse,
                                          directed = F)
g_mouse <- RedeR::subg(g=mouse_network,dat=dictionary_FAD,refcol = 1)

g_mouse <- att.setv(g=g_mouse,from = 'external_gene_name','nodeAlias')

mouse_ev <- igraph::get.data.frame(g_mouse,what = 'both')

mouse_edges <- mouse_ev$edges
mouse_vertice <- mouse_ev$vertices

mouse_from <- data.frame(name = mouse_edges$from,stringsAsFactors = F)
mouse_from <- merge(mouse_from,mouse_vertice,by='name')
mouse_from <- mouse_from[match(mouse_edges$from,
                               mouse_from$name),]

mouse_to <- data.frame(name = mouse_edges$to,stringsAsFactors = F)
mouse_to <- merge(mouse_to,mouse_vertice,by='name')
mouse_to <- mouse_to[match(mouse_edges$to,
                               mouse_to$name),]
mouse_edges <- data.frame(gene1symbol = mouse_from$external_gene_name,
                          gene2symbol = mouse_to$external_gene_name,stringsAsFactors = F)

saveRDS(mouse_edges,'CEMtool/PPI/mouse_PPI.rds')
