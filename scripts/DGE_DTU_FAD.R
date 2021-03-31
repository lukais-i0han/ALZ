#### DGE - DTU 
library('dplyr')
library('stringr')
library('RColorBrewer')
library('ggplot2')

FAD5X_all <- read.csv('archives/FAD5X_all.csv',header = T,stringsAsFactors = F,
             row.names = 1)

FAD_DTU <- readRDS('results/isoform_results/FAD5X_ISW_sig.rds')
write.csv(FAD_DTU,'archives/FAD_DTU.csv')

FAD5X_all <- FAD5X_all[FAD5X_all$DEG != 'DEG',]


FAD5X_all <- FAD5X_all[FAD5X_all$Gene_Symbol %in% FAD_DTU$gene_name,]
FAD5X_all <- FAD5X_all[FAD5X_all$sig != 'FDR<0.01',]

FAD_NDEG_CCX <- FAD5X_all[FAD5X_all$Area == 'Cortex',]
FAD_NDEG_HIP <- FAD5X_all[FAD5X_all$Area == 'Hippocampus',]

### Hippocampus

FAD_DTU_HIP <- FAD_DTU[FAD_DTU$Region == 'Hippocampus',]

FAD_NDEG_HIP$'DTU_4' <- FAD_NDEG_HIP$Gene_Symbol %in% 
  FAD_DTU_HIP[FAD_DTU_HIP$Group =='4',5]

FAD_NDEG_HIP$'DTU_12' <- FAD_NDEG_HIP$Gene_Symbol %in% 
  FAD_DTU_HIP[FAD_DTU_HIP$Group =='12',5]

FAD_NDEG_HIP$'DTU_18' <- FAD_NDEG_HIP$Gene_Symbol %in% 
  FAD_DTU_HIP[FAD_DTU_HIP$Group =='18',5]

FAD_NDEG_HIP <- FAD_NDEG_HIP %>% mutate(DTU_column = case_when(
  DTU_4 == 'TRUE' & DTU_12 == 'FALSE' & DTU_18 == 'FALSE' ~ '4',
  DTU_4 == 'TRUE' & DTU_12== 'TRUE' & DTU_18 == 'FALSE' ~ '4_12',
  DTU_4 == 'TRUE' & DTU_12 == 'TRUE' & DTU_18 == 'TRUE' ~ '4_12_18',
  DTU_4 == 'TRUE' & DTU_12 == 'FALSE' & DTU_18 == 'TRUE'~ '4_18',
  DTU_4 == 'FALSE' & DTU_12 == 'TRUE' & DTU_18 == 'FALSE' ~ '12',
  DTU_4 == 'FALSE' & DTU_12 == 'TRUE' & DTU_18 == 'TRUE' ~ '12_18',
  DTU_4 == 'FALSE' & DTU_12 == 'FALSE' & DTU_18 =='TRUE' ~ '18',
  DTU_4 == 'FALSE' & DTU_12 == 'FALSE' & DTU_18 =='FALSE' ~ 'Z'
))
 
FAD_NDEG_HIP <- FAD_NDEG_HIP[FAD_NDEG_HIP$DTU_column != 'Z',]

### Cortex

FAD_DTU_CCX <- FAD_DTU[FAD_DTU$Region == 'Cortex',]

FAD_NDEG_CCX$'DTU_4' <- FAD_NDEG_CCX$Gene_Symbol %in% 
  FAD_DTU_CCX[FAD_DTU_CCX$Group =='4',5]

FAD_NDEG_CCX$'DTU_12' <- FAD_NDEG_CCX$Gene_Symbol %in% 
  FAD_DTU_CCX[FAD_DTU_CCX$Group =='12',5]

FAD_NDEG_CCX$'DTU_18' <- FAD_NDEG_CCX$Gene_Symbol %in% 
  FAD_DTU_CCX[FAD_DTU_CCX$Group =='18',5]

FAD_NDEG_CCX <- FAD_NDEG_CCX%>% mutate(DTU_column = case_when(
  DTU_4 == 'TRUE' & DTU_12 == 'FALSE' & DTU_18 == 'FALSE' ~ '4',
  DTU_4 == 'TRUE' & DTU_12== 'TRUE' & DTU_18 == 'FALSE' ~ '4_12',
  DTU_4 == 'TRUE' & DTU_12 == 'TRUE' & DTU_18 == 'TRUE' ~ '4_12_18',
  DTU_4 == 'TRUE' & DTU_12 == 'FALSE' & DTU_18 == 'TRUE'~ '4_18',
  DTU_4 == 'FALSE' & DTU_12 == 'TRUE' & DTU_18 == 'FALSE' ~ '12',
  DTU_4 == 'FALSE' & DTU_12 == 'TRUE' & DTU_18 == 'TRUE' ~ '12_18',
  DTU_4 == 'FALSE' & DTU_12 == 'FALSE' & DTU_18 =='TRUE' ~ '18',
  DTU_4 == 'FALSE' & DTU_12 == 'FALSE' & DTU_18 =='FALSE' ~ 'Z'
))

FAD_NDEG_CCX <- FAD_NDEG_CCX[FAD_NDEG_CCX$DTU_column != 'Z',]

FAD_NDEG_CCX_4 <- FAD_NDEG_CCX[FAD_NDEG_CCX$Month == '4',]
FAD_NDEG_CCX_12 <- FAD_NDEG_CCX[FAD_NDEG_CCX$Month == '12',]
FAD_NDEG_CCX_18  <- FAD_NDEG_CCX[FAD_NDEG_CCX$Month == '18',]

NDEG_DTU_CCX <- data.frame(source=c(rep('[4]',115),rep('[12]',102),rep('[18]',108)),
                          destination = c(rep('4',76),
                                    rep('4:12',1),
                                    rep('4:12:18',1),
                                    rep('12',27),
                                    rep('12:18',6),
                                    rep('18',4),
                                    rep('4',73),
                                    rep('4:12',1),
                                    rep('4:12:18',1),
                                    rep('12',24),
                                    rep('12:18',1),
                                    rep('18',2),
                                    rep('4',78),
                                    rep('4:12',1),
                                    rep('4:12:18',1),
                                    rep('12',26),
                                    rep('12:18',1),
                                    rep('18',1)),stringsAsFactors = F)

NDEG_DTU_CCX$source <- factor(NDEG_DTU_CCX$source,levels=
                    c('[4]','[12]','[18]'))

sources <- NDEG_DTU_CCX  %>%
  distinct(source) %>% rename(label = source)

destinations <-NDEG_DTU_CCX  %>%
  distinct(destination) %>% rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <- NDEG_DTU_CCX  %>%
  group_by(source,destination) %>%
  summarise(weight = n()) %>%
  ungroup()


edges <- per_route %>%
  left_join(nodes,by = c("source" = "label")) %>%
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("destination" = "label")) %>% 
  rename(to = id)


edges <- select(edges, from, to, weight)

routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes,directed =T )

routes_igraph_tidy <- as_tbl_graph(routes_igraph)
routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)

arc_color <- rep(c('4-conection','12-conection','18-conection'),c(6,6,6))
arc_color <- factor(arc_color,levels=
                      c('4-conection','12-conection','18-conection'))

routes_tidy <- routes_tidy %>%
  activate(edges) %>% mutate(arc_color = arc_color)

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour= arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.5, 1.5)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "No-DGE - DTU") +
  theme_graph() + ggtitle('NDEG - CCX')



#####


FAD_NDEG_HIP_4 <- FAD_NDEG_HIP[FAD_NDEG_HIP$Month == '4',]

FAD_NDEG_HIP_12 <- FAD_NDEG_HIP[FAD_NDEG_HIP$Month == '12',]

FAD_NDEG_HIP_18  <- FAD_NDEG_HIP[FAD_NDEG_HIP$Month == '18',]

NDGE_DTU_HIP <- data.frame(source=c(rep('[4]',65),rep('[12]',52),rep('[18]',43)),
                          destination = c(rep('4',10),
                                    rep('4:12',1),
                                    rep('4:12:18',1),
                                    rep('12',16),
                                    rep('12:18',10),
                                    rep('18',27),
                                    rep('4',9),
                                    rep('4:12',1),
                                    rep('12',14),
                                    rep('12:18',4),
                                    rep('18',24),
                                    rep('4',9),
                                    rep('12',13),
                                    rep('12:18',1),
                                    rep('18',20)),stringsAsFactors = F)

NDGE_DTU_HIP$source <- factor(NDGE_DTU_HIP$source,levels=
                  c('[4]','[12]','[18]'))

sources <- NDGE_DTU_HIP  %>%
  distinct(source) %>% 
  rename(label = source)

destinations <-NDGE_DTU_HIP  %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <- NDGE_DTU_HIP  %>%
  group_by(source,destination) %>%
  summarise(weight = n()) %>%
  ungroup()


edges <- per_route %>%
  left_join(nodes,by = c("source" = "label")) %>%
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("destination" = "label")) %>% 
  rename(to = id)


edges <- select(edges, from, to, weight)

routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes,directed =T )

routes_igraph_tidy <- as_tbl_graph(routes_igraph)
routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)


arc_color_HIP <- rep(c('4-conection','12-conection','18-conection'),c(6,5,4))
arc_color_HIP <- factor(arc_color_HIP,levels = 
                        c('4-conection','12-conection','18-conection'))

routes_tidy <- routes_tidy %>%
  activate(edges) %>% mutate(arc_color = arc_color_HIP)

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.1, 0.5)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "No-DGE - DTU") +
  theme_graph() + ggtitle('NDEG - HIP')


#### Sig Analysis

FAD5X_sig <- read.csv('archives/FAD5X_sig.csv',header = T, stringsAsFactors = F,
                      row.names = 1)

FAD_DEG_CCX <- FAD5X_sig[FAD5X_sig$Area == 'Cortex',]
FAD_DEG_HIP <- FAD5X_sig[FAD5X_sig$Area == 'Hippocampus',]

FAD_DTU <- readRDS('results/isoform_results/FAD5X_ISW_sig.rds')

FAD_DTU_CCX <- FAD_DTU[FAD_DTU$Region == 'Cortex',]
FAD_DTU_HIP <- FAD_DTU[FAD_DTU$Region == 'Hippocampus',]

### Hippocampus

FAD_DEG_HIP$'DTU_4' <- FAD_DEG_HIP$Gene_Symbol %in% 
  FAD_DTU_HIP[FAD_DTU_HIP$Group =='4',5]

FAD_DEG_HIP$'DTU_12' <- FAD_DEG_HIP$Gene_Symbol %in% 
  FAD_DTU_HIP[FAD_DTU_HIP$Group =='12',5]

FAD_DEG_HIP$'DTU_18' <- FAD_DEG_HIP$Gene_Symbol %in% 
  FAD_DTU_HIP[FAD_DTU_HIP$Group =='18',5]

FAD_DEG_HIP <- FAD_DEG_HIP %>% mutate(DTU_column = case_when(
  DTU_4 == 'TRUE' & DTU_12 == 'FALSE' & DTU_18 == 'FALSE' ~ '4',
  DTU_4 == 'TRUE' & DTU_12== 'TRUE' & DTU_18 == 'FALSE' ~ '4_12',
  DTU_4 == 'TRUE' & DTU_12 == 'TRUE' & DTU_18 == 'TRUE' ~ '4_12_18',
  DTU_4 == 'TRUE' & DTU_12 == 'FALSE' & DTU_18 == 'TRUE'~ '4_18',
  DTU_4 == 'FALSE' & DTU_12 == 'TRUE' & DTU_18 == 'FALSE' ~ '12',
  DTU_4 == 'FALSE' & DTU_12 == 'TRUE' & DTU_18 == 'TRUE' ~ '12_18',
  DTU_4 == 'FALSE' & DTU_12 == 'FALSE' & DTU_18 =='TRUE' ~ '18',
  DTU_4 == 'FALSE' & DTU_12 == 'FALSE' & DTU_18 =='FALSE' ~ 'Z'
))

FAD_DEG_HIP <- FAD_DEG_HIP[FAD_DEG_HIP$DTU_column != 'Z',]

### Cortex

FAD_DEG_CCX$'DTU_4' <- FAD_DEG_CCX$Gene_Symbol %in% 
  FAD_DTU_CCX[FAD_DTU_CCX$Group =='4',5]

FAD_DEG_CCX$'DTU_12' <- FAD_DEG_CCX$Gene_Symbol %in% 
  FAD_DTU_CCX[FAD_DTU_CCX$Group =='12',5]

FAD_DEG_CCX$'DTU_18' <- FAD_DEG_CCX$Gene_Symbol %in% 
  FAD_DTU_CCX[FAD_DTU_CCX$Group =='18',5]

FAD_DEG_CCX <- FAD_DEG_CCX %>% mutate(DTU_column = case_when(
  DTU_4 == 'TRUE' & DTU_12 == 'FALSE' & DTU_18 == 'FALSE' ~ '4',
  DTU_4 == 'TRUE' & DTU_12== 'TRUE' & DTU_18 == 'FALSE' ~ '4_12',
  DTU_4 == 'TRUE' & DTU_12 == 'TRUE' & DTU_18 == 'TRUE' ~ '4_12_18',
  DTU_4 == 'TRUE' & DTU_12 == 'FALSE' & DTU_18 == 'TRUE'~ '4_18',
  DTU_4 == 'FALSE' & DTU_12 == 'TRUE' & DTU_18 == 'FALSE' ~ '12',
  DTU_4 == 'FALSE' & DTU_12 == 'TRUE' & DTU_18 == 'TRUE' ~ '12_18',
  DTU_4 == 'FALSE' & DTU_12 == 'FALSE' & DTU_18 =='TRUE' ~ '18',
  DTU_4 == 'FALSE' & DTU_12 == 'FALSE' & DTU_18 =='FALSE' ~ 'Z'
))

FAD_DEG_CCX <- FAD_DEG_CCX[FAD_DEG_CCX$DTU_column != 'Z',]


### Graphics

###Hippocampus
FAD_DEG_HIP_4 <- FAD_DEG_HIP[FAD_DEG_HIP$Month == '4',]

FAD_DEG_HIP_12 <- FAD_DEG_HIP[FAD_DEG_HIP$Month == '12',]

FAD_DEG_HIP_18  <- FAD_DEG_HIP[FAD_DEG_HIP$Month == '18',]


DGE_DTU_HIP <- data.frame(source=c(rep('[4]',6),rep('[12]',19),rep('[18]',28)),
                          destination = c(rep('4',1),
                                    rep('4:12:18',2),
                                    rep('12',2),
                                    rep('18',1),
                                    rep('4',2),
                                    rep('4:12:18',3),
                                    rep('12',4),
                                    rep('12:18',6),
                                    rep('18',4),
                                    rep('4',2),
                                    rep('4:12',1),
                                    rep('4:12:18',3),
                                    rep('12',5),
                                    rep('12:18',9),
                                    rep('18',8)),stringsAsFactors = F)


DGE_DTU_HIP$source <- factor(DGE_DTU_HIP$source,levels=
                                c('[4]','[12]','[18]'))


sources <- DGE_DTU_HIP  %>%
  distinct(source) %>% 
  rename(label = source)

destinations <- DGE_DTU_HIP  %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <- DGE_DTU_HIP  %>%
  group_by(source,destination) %>%
  summarise(weight = n()) %>%
  ungroup()


edges <- per_route %>%
  left_join(nodes,by = c("source" = "label")) %>%
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("destination" = "label")) %>% 
  rename(to = id)


edges <- select(edges, from, to, weight)

routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes,directed =T )

routes_igraph_tidy <- as_tbl_graph(routes_igraph)
routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
routes_tidy <- routes_tidy %>%
  activate(edges) %>% mutate(arc_color = 
          rep(c('4-conection','12-conection','18-conection'),c(4,5,6)))

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "DGE - DTU") +
  theme_graph() + ggtitle('DEG-HIP')





### Cortex

FAD_DEG_CCX_4 <- FAD_DEG_CCX[FAD_DEG_CCX$Month == '4',]
FAD_DEG_CCX_12 <- FAD_DEG_CCX[FAD_DEG_CCX$Month == '12',]
FAD_DEG_CCX_18  <- FAD_DEG_CCX[FAD_DEG_CCX$Month == '18',]

DGE_DTU_CCX <- data.frame(source=c(rep('[4]',4),rep('[12]',12),rep('[18]',11)),
                          destination = c(rep('4',2),
                                    rep('4:12:18',1),
                                    rep('12',1),
                                    rep('4',1),
                                    rep('4:12:18',1),
                                    rep('12',3),
                                    rep('12:18',5),
                                    rep('18',2),
                                    rep('4:12:18',1),
                                    rep('12',2),
                                    rep('12:18',5),
                                    rep('18',3)),stringsAsFactors = F)


DGE_DTU_CCX$source <- factor(DGE_DTU_CCX$source,levels=
                               c('[4]','[12]','[18]'))


sources <- DGE_DTU_CCX %>%
  distinct(source) %>% 
  rename(label = source)

destinations <- DGE_DTU_CCX  %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <- DGE_DTU_CCX  %>%
  group_by(source,destination) %>%
  summarise(weight = n()) %>%
  ungroup()


edges <- per_route %>%
  left_join(nodes,by = c("source" = "label")) %>%
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("destination" = "label")) %>% 
  rename(to = id)


edges <- select(edges, from, to, weight)

routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes,directed =T )

routes_igraph_tidy <- as_tbl_graph(routes_igraph)
routes_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
routes_tidy <- routes_tidy %>%
  activate(edges) %>% mutate(arc_color = 
          rep(c('4-conection','12-conection','18-conection'),c(3,5,4)))

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "DGE - DTU") +
  theme_graph() +ggtitle('DEG-CCX')


### Tables DGE DTU

FAD_NDEG <- rbind(FAD_NDEG_CCX_4,FAD_NDEG_CCX_12,
                  FAD_NDEG_CCX_18,FAD_NDEG_HIP_4,
                  FAD_NDEG_HIP_12,FAD_NDEG_HIP_18)

FAD_DEG <- rbind(FAD_DEG_CCX_4,FAD_DEG_CCX_12,
                  FAD_DEG_CCX_18,FAD_DEG_HIP_4,
                  FAD_DEG_HIP_12,FAD_DEG_HIP_18)                 
FAD_DEG_DTU <- list(FAD_NDEG = FAD_NDEG,
                    FAD_DEG = FAD_DEG)

saveRDS(FAD_DEG_DTU,'results/FAD_DEG_DTU.rds')
