
# library
library(igraph)

### Mayo

links <- data.frame(
  source=c('A','A','A','B','B','B','C','C','C'),
  target=c('A:B','A:C','A:B:C','A:B','B:C','A:B:C','A:C','B:C','A:B:C'),
  importance=(c(1,1,2,1,1,2,1,1,2))
)


network <- graph_from_data_frame(d=links, directed=T) 

coul  <-c('#E41A1C','#E41A1C','#E41A1C','#FFFF33','#FFFF33','#FF7F00',
          '#FFFF33','#FFFF33','#FFFF33')

plot(network, vertex.color=coul,layout=layout.fruchterman.reingold)

### MSBB
links_MSBB <- data.frame(
  source = c('A','B'),
  target = c('A:B','A:B')
)

network_MSBB <- graph_from_data_frame(d=links_MSBB,directed = T)

coul_MSBB <- c('#E41A1C','#FFFF33')

plot(network_MSBB, vertex.color=coul_MSBB,layout=layout.fruchterman.reingold)
