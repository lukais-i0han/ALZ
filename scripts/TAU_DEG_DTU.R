####

Tau_all <- read.csv('archives/Tau_all.csv',header = T,stringsAsFactors = F,
                    row.names = 1)

Tau_DTU <- readRDS('results/isoform_results/Tau_ISW_sig.rds')
write.csv(Tau_DTU,'archives/TAU_DTU.csv')

Tau_all <- Tau_all[Tau_all$DEG != 'DEG',]

Tau_all <- Tau_all[Tau_all$Gene_Name %in% Tau_DTU$gene_name,]

Tau_4 <- Tau_all[Tau_all$Month == '4',]
Tau_17 <- Tau_all[Tau_all$Month == '17',]


Tau_4$'DTU_4' <- Tau_4$Gene_Name %in%
  Tau_DTU[Tau_DTU$Group == '4',5]

Tau_4$'DTU_17' <- Tau_4$Gene_Name %in%
  Tau_DTU[Tau_DTU$Group == '17',5]

Tau_4$'DTU_column' <- ifelse(Tau_4$DTU_4 == 'FALSE' & Tau_4$DTU_17 == 'TRUE',
                      '17','4')


Tau_17$'DTU_4' <- Tau_17$Gene_Name %in%
  Tau_DTU[Tau_DTU$Group == '4',5]

Tau_17$'DTU_17' <- Tau_17$Gene_Name %in%
  Tau_DTU[Tau_DTU$Group == '17',5]

Tau_17$'DTU_column' <- ifelse(Tau_17$DTU_4 == 'FALSE' & Tau_17$DTU_17 == 'TRUE',
                        '17','4')

### Graphics

Tau_4_17_NDEG <- data.frame(source = c(rep('[4]',9),rep('[17]',80)),
                            destination = c(rep('17',9),rep('4',5),rep('17',75)),
              stringsAsFactors = F)

Tau_4_17_NDEG$source <- factor(Tau_4_17_NDEG$source,
                                  levels = c('[4]','[17]'))


sources <- Tau_4_17_NDEG  %>%
  distinct(source) %>% 
  rename(label = source)

destinations <-Tau_4_17_NDEG  %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <-Tau_4_17_NDEG  %>%
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
    rep(c('4-conection','17-conection'),c(1,2)))

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "No-DGE - DTU") +
  theme_graph() + ggtitle('NDEG - TauD35')



#### 

Tau_sig <- read.csv('archives/Tau_sig.csv',header = T,stringsAsFactors = F,
           row.names = 1) 

Tau_DTU <- readRDS('results/isoform_results/Tau_ISW_sig.rds')

Tau_sig_4 <- Tau_sig[Tau_sig$Month == '4',]
Tau_sig_17 <- Tau_sig[Tau_sig$Month == '17',]

Tau_sig_4$'DTU_4' <- Tau_sig_4$Gene_Name %in%
  Tau_DTU[Tau_DTU$Group == '4',5]

Tau_sig_4$'DTU_17' <- Tau_sig_4$Gene_Name %in%
  Tau_DTU[Tau_DTU$Group == '17',5]

Tau_sig_17$'DTU_4' <- Tau_sig_17$Gene_Name %in%
  Tau_DTU[Tau_DTU$Group == '4',5]

Tau_sig_17$'DTU_17' <- Tau_sig_17$Gene_Name %in%
  Tau_DTU[Tau_DTU$Group == '17',5]


Tau_4_17_DEG <- data.frame(source = c(rep("[4]",1),rep('[17]',10)),
                           destination = c('4:17',rep('17',10)),stringsAsFactors = F
                           )

Tau_4_17_DEG$source <- factor(Tau_4_17_DEG$source,
                                   levels = c("[4]",'[17]'))

sources <- Tau_4_17_DEG  %>%
  distinct(source) %>% 
  rename(label = source)

destinations <-Tau_4_17_DEG  %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <-Tau_4_17_DEG  %>%
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
          rep(c('4-conection','17-conection'),c(1,1)))

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "No-DGE - DTU") +
  theme_graph() + ggtitle('DEG - TauD35')
