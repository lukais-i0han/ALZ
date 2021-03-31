#### DGE - DTU 
library('dplyr')
library('stringr')
library('RColorBrewer')
library('ggplot2')
library('tidygraph')
library('ggraph')
library('igraph')
library('tibble')
library('tidyr')
library('networkD3')

ALZ_all <- read.csv('archives/ALZ_all.csv', header = T, stringsAsFactors = F,
                    row.names = 1)


PA_all <- read.csv('archives/PA_all.csv',header = T, stringsAsFactors = F,
                   row.names = 1)

PSP_all <- read.csv('archives/PSP_all.csv', header = T, stringsAsFactors = F,
                    row.names = 1)

MAYO_DTU <- readRDS('results/isoform_results/mayo_ISW_sig.rds')

ALZ_all <- ALZ_all %>% dplyr::select(ENSG,Gene_Symbol,DEG,Age_Group,Group) %>%
  group_by(Age_Group)

ALZ_all <- ALZ_all[ALZ_all$DEG != 'DEG',]

ALZ_DTU <- MAYO_DTU[MAYO_DTU$condition_1 == 'AD',] 

ALZ_all <- ALZ_all[ALZ_all$Gene_Symbol %in% ALZ_DTU$gene_name,]
ALZ_NDGE_A <- ALZ_all[ALZ_all$Age_Group == 'A',]
ALZ_NDGE_B <- ALZ_all[ALZ_all$Age_Group == 'B',]
ALZ_NDGE_C <- ALZ_all[ALZ_all$Age_Group == 'C',]

ALZ_DTU_A <- ALZ_DTU[ALZ_DTU$AgeAtDeath == 'A',]
ALZ_DTU_B <- ALZ_DTU[ALZ_DTU$AgeAtDeath == 'B',]
ALZ_DTU_C <- ALZ_DTU[ALZ_DTU$AgeAtDeath == 'C',]

ALZ_NDGE_A$'A_DTU' <- ALZ_NDGE_A$Gene_Symbol %in% ALZ_DTU_A$gene_name
ALZ_NDGE_A$'B_DTU' <- ALZ_NDGE_A$Gene_Symbol %in% ALZ_DTU_B$gene_name
ALZ_NDGE_A$'C_DTU' <- ALZ_NDGE_A$Gene_Symbol %in% ALZ_DTU_C$gene_name
ALZ_NDGE_A <- ALZ_NDGE_A %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'AB',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'ABC',
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'TRUE'~ 'AC',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'B',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'BC',
  A_DTU == 'FALSE' & B_DTU == 'FALSE' & C_DTU =='TRUE' ~ 'C'
))


ALZ_NDGE_B$'A_DTU' <- ALZ_NDGE_B$Gene_Symbol %in% ALZ_DTU_A$gene_name
ALZ_NDGE_B$'B_DTU' <- ALZ_NDGE_B$Gene_Symbol %in% ALZ_DTU_B$gene_name
ALZ_NDGE_B$'C_DTU' <- ALZ_NDGE_B$Gene_Symbol %in% ALZ_DTU_C$gene_name

ALZ_NDGE_B <- ALZ_NDGE_B %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'AB',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'ABC',
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'TRUE'~ 'AC',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'B',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'BC',
  A_DTU == 'FALSE' & B_DTU == 'FALSE' & C_DTU =='TRUE' ~ 'C'
))

ALZ_NDGE_C$'A_DTU' <- ALZ_NDGE_C$Gene_Symbol %in% ALZ_DTU_A$gene_name
ALZ_NDGE_C$'B_DTU' <- ALZ_NDGE_C$Gene_Symbol %in% ALZ_DTU_B$gene_name
ALZ_NDGE_C$'C_DTU' <- ALZ_NDGE_C$Gene_Symbol %in% ALZ_DTU_C$gene_name

ALZ_NDGE_C <- ALZ_NDGE_C %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'AB',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'ABC',
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'TRUE'~ 'AC',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'B',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'BC',
  A_DTU == 'FALSE' & B_DTU == 'FALSE' & C_DTU =='TRUE' ~ 'C'
))

DGE_DTU_ABC <- data.frame(source = c(rep('[A]',1437),rep('[B]',1395),rep('[C]',1335)),
                        destination = c(rep('A',45),
                                  rep('AB',16),
                                  rep('ABC',22),
                                  rep('AC',10),
                                  rep('B',148),
                                  rep('BC',272),
                                  rep('C',924),
                                  rep('A',43),
                                  rep('AB',13),
                                  rep('ABC',18),
                                  rep('AC',7),
                                  rep('B',142),
                                  rep('BC',262),
                                  rep('C',910),
                                  rep('A',45),
                                  rep('AB',15),
                                  rep('ABC',23),
                                  rep('AC',10),
                                  rep('B',148),
                                  rep('BC',269),
                                  rep('C',825)),stringsAsFactors = F)

sources <- DGE_DTU_ABC %>%
  distinct(source) %>% 
  rename(label = source)
    
destinations <- DGE_DTU_ABC %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <- DGE_DTU_ABC %>%
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
  rep(c('A-conection','B-conection','C-conection'),c(7,7,7)))


ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "No-DGE - DTU") +
  theme_graph()


#### PA

PA_all <- PA_all %>% dplyr::select(ENSG,Gene_Symbol,DEG,Age_Group,Group) %>%
  group_by(Age_Group)

PA_all <- PA_all[PA_all$DEG != 'DEG',]

PA_DTU <- MAYO_DTU[MAYO_DTU$condition_1 == 'PA',] 

PA_all <- PA_all[PA_all$Gene_Symbol %in% PA_DTU$gene_name,]
PA_NDGE_A <- PA_all[PA_all$Age_Group == 'A',]
PA_NDGE_B <- PA_all[PA_all$Age_Group == 'B',]
PA_NDGE_C <- PA_all[PA_all$Age_Group == 'C',]

PA_DTU_A <- PA_DTU[PA_DTU$AgeAtDeath == 'A',]
PA_DTU_B <- PA_DTU[PA_DTU$AgeAtDeath == 'B',]
PA_DTU_C <- PA_DTU[PA_DTU$AgeAtDeath == 'C',]

PA_NDGE_A$'A_DTU' <- PA_NDGE_A$Gene_Symbol %in% PA_DTU_A$gene_name
PA_NDGE_A$'B_DTU' <- PA_NDGE_A$Gene_Symbol %in% PA_DTU_B$gene_name
PA_NDGE_A$'C_DTU' <- PA_NDGE_A$Gene_Symbol %in% PA_DTU_C$gene_name
PA_NDGE_A <- PA_NDGE_A %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'AB',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'ABC',
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'TRUE'~ 'AC',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'B',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'BC',
  A_DTU == 'FALSE' & B_DTU == 'FALSE' & C_DTU =='TRUE' ~ 'C'
))


PA_NDGE_B$'A_DTU' <- PA_NDGE_B$Gene_Symbol %in% PA_DTU_A$gene_name
PA_NDGE_B$'B_DTU' <- PA_NDGE_B$Gene_Symbol %in% PA_DTU_B$gene_name
PA_NDGE_B$'C_DTU' <- PA_NDGE_B$Gene_Symbol %in% PA_DTU_C$gene_name

PA_NDGE_B <- PA_NDGE_B %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'AB',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'ABC',
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'TRUE'~ 'AC',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'B',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'BC',
  A_DTU == 'FALSE' & B_DTU == 'FALSE' & C_DTU =='TRUE' ~ 'C'
))

PA_NDGE_C$'A_DTU' <- PA_NDGE_C$Gene_Symbol %in% PA_DTU_A$gene_name
PA_NDGE_C$'B_DTU' <- PA_NDGE_C$Gene_Symbol %in% PA_DTU_B$gene_name
PA_NDGE_C$'C_DTU' <- PA_NDGE_C$Gene_Symbol %in% PA_DTU_C$gene_name

PA_NDGE_C <- PA_NDGE_C %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'AB',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'ABC',
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'TRUE'~ 'AC',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'B',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'BC',
  A_DTU == 'FALSE' & B_DTU == 'FALSE' & C_DTU =='TRUE' ~ 'C'
))

DGE_DTU_ABC <- data.frame(source = c(rep('[A]',1031),rep('[B]',1031),rep('[C]',555)),
                          destination = c(rep('A',27),
                                          rep('AC',10),
                                          rep('B',41),
                                          rep('BC',25),
                                          rep('C',928),
                                          rep('A',27),
                                          rep('AC',10),
                                          rep('B',41),
                                          rep('BC',25),
                                          rep('C',928),
                                          rep('A',20),
                                          rep('AC',9),
                                          rep('B',29),
                                          rep('BC',20),
                                          rep('C',477)),stringsAsFactors = F)

sources <- DGE_DTU_ABC %>%
  dplyr::distinct(source) %>% 
  rename(label = source)

destinations <- DGE_DTU_ABC %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <- DGE_DTU_ABC %>%
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
          rep(c('A-conection','B-conection','C-conection'),c(5,5,5)))


ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "Number of Transitions") +
  theme_graph()







### PSP

PSP_all <- PSP_all %>% select(ENSG,Gene_Symbol,DEG,Age_Group,Group) %>%
  group_by(Age_Group)

PSP_all <- PSP_all[PSP_all$DEG != 'DEG',]

PSP_DTU <- MAYO_DTU[MAYO_DTU$condition_1 == 'PSP',] 

PSP_all <- PSP_all[PSP_all$Gene_Symbol %in% PSP_DTU$gene_name,]
PSP_NDGE_A <- PSP_all[PSP_all$Age_Group == 'A',]
PSP_NDGE_B <- PSP_all[PSP_all$Age_Group == 'B',]

PSP_DTU_A <- PSP_DTU[PSP_DTU$AgeAtDeath == 'A',]
PSP_DTU_B <- PSP_DTU[PSP_DTU$AgeAtDeath == 'B',]

PSP_NDGE_A$'A_DTU' <- PSP_NDGE_A$Gene_Symbol %in% PSP_DTU_A$gene_name
PSP_NDGE_A$'B_DTU' <- PSP_NDGE_A$Gene_Symbol %in% PSP_DTU_B$gene_name
PSP_NDGE_A <- PSP_NDGE_A %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' ~ 'AB',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' ~ 'B'
))


PSP_NDGE_B$'A_DTU' <- PSP_NDGE_B$Gene_Symbol %in% PSP_DTU_A$gene_name
PSP_NDGE_B$'B_DTU' <- PSP_NDGE_B$Gene_Symbol %in% PSP_DTU_B$gene_name
PSP_NDGE_B <- PSP_NDGE_B %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' ~ 'AB',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' ~ 'B'
))



DGE_DTU_AB <- data.frame(source = c(rep('[A]',269),rep('[B]',311)),
                          destination = c(rep('A',173),
                                    rep('AB',13),
                                    rep('B',83),
                                    rep('A',204),
                                    rep('AB',15),
                                    rep('B',92)),stringsAsFactors = F)

sources <- DGE_DTU_AB %>%
  distinct(source) %>% 
  rename(label = source)

destinations <- DGE_DTU_AB %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <- DGE_DTU_AB %>%
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
      rep(c('A-conection','B-conection'),c(3,3)))

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "Number of Transitions") +
  theme_graph()


#### ALZ Sig

ALZ_DEG <- ALZ_DEG %>% select(ENSG,Gene_Symbol,DEG,Age_Group) %>%
  group_by(Age_Group)


ALZ_DTU <- MAYO_DTU[MAYO_DTU$condition_1 == 'AD',] 

ALZ_DEG <- ALZ_DEG[ALZ_DEG$Gene_Symbol %in% ALZ_DTU$gene_name,]
ALZ_DEG_A <- ALZ_DEG[ALZ_DEG$Age_Group == 'A',]
ALZ_DEG_B <- ALZ_DEG[ALZ_DEG$Age_Group == 'B',]
ALZ_DEG_C <- ALZ_DEG[ALZ_DEG$Age_Group == 'C',]

ALZ_DTU_A <- ALZ_DTU[ALZ_DTU$AgeAtDeath == 'A',]
ALZ_DTU_B <- ALZ_DTU[ALZ_DTU$AgeAtDeath == 'B',]
ALZ_DTU_C <- ALZ_DTU[ALZ_DTU$AgeAtDeath == 'C',]

ALZ_DEG_A$'A_DTU' <- ALZ_DEG_A$Gene_Symbol %in% ALZ_DTU_A$gene_name
ALZ_DEG_A$'B_DTU' <- ALZ_DEG_A$Gene_Symbol %in% ALZ_DTU_B$gene_name
ALZ_DEG_A$'C_DTU' <- ALZ_DEG_A$Gene_Symbol %in% ALZ_DTU_C$gene_name
ALZ_DEG_A <- ALZ_DEG_A %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'AB',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'ABC',
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'TRUE'~ 'AC',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'B',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'BC',
  A_DTU == 'FALSE' & B_DTU == 'FALSE' & C_DTU =='TRUE' ~ 'C'
))


ALZ_DEG_B$'A_DTU' <- ALZ_DEG_B$Gene_Symbol %in% ALZ_DTU_A$gene_name
ALZ_DEG_B$'B_DTU' <- ALZ_DEG_B$Gene_Symbol %in% ALZ_DTU_B$gene_name
ALZ_DEG_B$'C_DTU' <- ALZ_DEG_B$Gene_Symbol %in% ALZ_DTU_C$gene_name
ALZ_DEG_B <- ALZ_DEG_B %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'AB',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'ABC',
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'TRUE'~ 'AC',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'B',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'BC',
  A_DTU == 'FALSE' & B_DTU == 'FALSE' & C_DTU =='TRUE' ~ 'C'
))

ALZ_DEG_C$'A_DTU' <- ALZ_DEG_C$Gene_Symbol %in% ALZ_DTU_A$gene_name
ALZ_DEG_C$'B_DTU' <- ALZ_DEG_C$Gene_Symbol %in% ALZ_DTU_B$gene_name
ALZ_DEG_C$'C_DTU' <- ALZ_DEG_C$Gene_Symbol %in% ALZ_DTU_C$gene_name
ALZ_DEG_C <- ALZ_DEG_C %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'AB',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'ABC',
  A_DTU == 'TRUE' & B_DTU == 'FALSE' & C_DTU == 'TRUE'~ 'AC',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'FALSE' ~ 'B',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' & C_DTU == 'TRUE' ~ 'BC',
  A_DTU == 'FALSE' & B_DTU == 'FALSE' & C_DTU =='TRUE' ~ 'C'
))

DGE_DTU_ALZ <- data.frame(source = c(rep('[A]',47),rep('[B]',88),rep('[C]',149)),
                          destination = c(rep('A',6),
                                    rep('ABC',2),
                                    rep('B',3),
                                    rep('BC',8),
                                    rep('C',28),
                                    rep('A',8),
                                    rep('AB',3),
                                    rep('ABC',6),
                                    rep('AC',3),
                                    rep('B',9),
                                    rep('BC',18),
                                    rep('C',41),
                                    rep('A',6),
                                    rep('AB',1),
                                    rep('ABC',1),
                                    rep('B',3),
                                    rep('BC',11),
                                    rep('C',127)),stringsAsFactors = F)

sources <- DGE_DTU_ALZ %>%
  distinct(source) %>% 
  rename(label = source)

destinations <- DGE_DTU_ALZ %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <- DGE_DTU_ALZ %>%
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
              rep(c('A-conection','B-conection','C-conection'),c(5,7,6)))

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "DGE - DTU") +
  theme_graph()

### PSP sig

PSP_DEG <- PSP_DEG %>% select(ENSG,Gene_Symbol,DEG,Age_Group) %>%
  group_by(Age_Group)

PSP_DTU <- MAYO_DTU[MAYO_DTU$condition_1 == 'PSP',] 

PSP_DEG <- PSP_DEG[PSP_DEG$Gene_Symbol %in% PSP_DTU$gene_name,]
PSP_DEG_A <- PSP_DEG[PSP_DEG$Age_Group == 'A',]
PSP_DEG_B <- PSP_DEG[PSP_DEG$Age_Group == 'B',]

PSP_DTU_A <- PSP_DTU[PSP_DTU$AgeAtDeath == 'A',]
PSP_DTU_B <- PSP_DTU[PSP_DTU$AgeAtDeath == 'B',]

PSP_DEG_A$'A_DTU' <- PSP_DEG_A$Gene_Symbol %in% PSP_DTU_A$gene_name
PSP_DEG_A$'B_DTU' <- PSP_DEG_A$Gene_Symbol %in% PSP_DTU_B$gene_name
PSP_DEG_A <- PSP_DEG_A %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' ~ 'AB',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' ~ 'B'
))


PSP_DEG_B$'A_DTU' <- PSP_DEG_B$Gene_Symbol %in% PSP_DTU_A$gene_name
PSP_DEG_B$'B_DTU' <- PSP_DEG_B$Gene_Symbol %in% PSP_DTU_B$gene_name
PSP_DEG_B <- PSP_DEG_B %>% mutate(DTU_column = case_when(
  A_DTU == 'TRUE' & B_DTU == 'FALSE' ~ 'A',
  A_DTU == 'TRUE' & B_DTU == 'TRUE' ~ 'AB',
  A_DTU == 'FALSE' & B_DTU == 'TRUE' ~ 'B'
))



DGE_DTU_PSP <- data.frame(source = c(rep('[A]',48),rep('[B]',6)),
                         destination = c(rep('A',34),
                                   rep('AB',2),
                                   rep('B',12),
                                   rep('A',3),
                                   rep('B',3)),stringsAsFactors = F)



sources <- DGE_DTU_PSP %>%
  distinct(source) %>% 
  rename(label = source)

destinations <- DGE_DTU_PSP %>%
  distinct(destination) %>%
  rename(label = destination)

nodes <- full_join(sources,destinations,by='label')
nodes <- nodes %>% rowid_to_column("id")

per_route <- DGE_DTU_PSP %>%
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
                               rep(c('A-conection','B-conection'),c(3,2)))

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = weight,colour=arc_color), alpha = 0.8,fold = T) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label)) +
  labs(edge_width = "DGE - DTU") +
  theme_graph()


##### tables of DEG-DTU

ALZ_NDEG <- rbind(ALZ_NDGE_A,ALZ_NDGE_B,ALZ_NDGE_C)
PSP_NDEG <- rbind(PSP_NDGE_A,PSP_NDGE_B)

ALZ_DEG <- rbind(ALZ_DEG_A,ALZ_DEG_B,ALZ_DEG_C)
PSP_DEG <- rbind(PSP_DEG_A,PSP_DEG_B)

mayo_DEG_DTU <- list(ALZ_NDEG = ALZ_NDEG,
                     PSP_NDEG = PSP_NDEG,
                     ALZ_DEG = ALZ_DEG,
                     PSP_DEG = PSP_DEG)
saveRDS(mayo_DEG_DTU,'results/mayo_DEG_DTU.rds')

PA_NDEG <- rbind(PA_NDGE_A,PA_NDGE_B,PA_NDGE_C)
write.csv(PA_NDEG,'results/PA_NDEG.csv')
