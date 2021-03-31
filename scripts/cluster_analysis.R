#### Anaçysi of clustering of cytoscape
library('stringr')
library('dplyr')
library('reshape2')
library('RedeR')
library('RColorBrewer')

ALZ_PSP_V <- read.csv('results/transcriptogramer/clustered_fuzzy/ALZ_PSP_V_fuzzy.csv') 

ALZ_PSP_E <- read.csv('results/transcriptogramer/clustered_fuzzy/ALZ_PSP_clustered_edge.csv',
                      header = T, stringsAsFactors = F)

ALZ_PSP_E <-as.data.frame(apply(ALZ_PSP_E,2,function(x)gsub('\\s+', '',x)))

ALZ_PSP_E <- colsplit(ALZ_PSP_E$name,pattern = '\\s*\\([^\\)]+\\)',
              names = c('from','x','to'))
ALZ_PSP_E <- ALZ_PSP_E[,-3]
colnames(ALZ_PSP_E) <- c('protein1','protein2')

ALZ_PSP_V <- ALZ_PSP_V %>% mutate(Age_Group = case_when(
  Age_Group == 'ALZ-A_B_C' ~ 'ALZ[A-B-C]',
  Age_Group == 'ALZ-C' ~ 'ALZ[C]',
  Age_Group == 'PSP-A' ~ 'PSP[A]',
  TRUE ~ 'PSP[A-B]'
))

ALZ_PSP_3_7_6 <- ALZ_PSP_V[ALZ_PSP_V$Fuzzy_Cluster == '3' |
                          ALZ_PSP_V$Fuzzy_Cluster == '6'|
                          ALZ_PSP_V$Fuzzy_Cluster == '7',]

g <- igraph::graph.data.frame(d=ALZ_PSP_E,
                              directed = F)

g_human<- RedeR::subg(g=g, dat = ALZ_PSP_3_7_6,
                   refcol = 9,
        maincomp = F,connected = F)
g_human<- att.setv(g = g_human, from='name', to="nodeAlias")
g_human<- att.setv(g = g_human,from="Age_Group", to="nodeColor",
                  cols = c('#7FC97F','#FEDF8F','#F0027F','#6990AA'))
g_human <- att.setv(g_human,from='degree.layout',to='nodeSize')

saveRDS(g_human,'results/transcriptogramer/ALZ_PSP_igraph.rds')

### RedeR

rdp <- RedPort() 
calld(rdp)


addGraph(rdp,g_human,gcoord=c(10,25),isNest=TRUE, theme='tm1', zoom=30)


relax(rdp,p2=400,p5=30,ps=T)
scl <- g_human$legNodeColor$scale 
leg <- g_human$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")

scl <- g_human$legNodeSize$scale
leg <- g_human$legNodeSize$legend
addLegend.size(rdp, sizevec=scl, labvec=leg, title="")

resetd(rdp)


### Mouse
CCX_HIP_glastly <- read.csv('results/transcriptogramer/cytoscape/CCX_HIP_V_glusty.csv',
              header = T, stringsAsFactors = F)
CCX_HIP_E <- read.csv('results/transcriptogramer/igraph/CCX_HIP_E.csv',header = T,
                      stringsAsFactors = F)

CCX_HIP_glastly <- CCX_HIP_glastly %>% mutate(Month_Group = case_when(
  Month_Group == 'Cortex-12' ~ 'Cortex[12]',
  Month_Group == 'Cortex-18' ~ 'Cortex[18]',
  Month_Group == 'Cortex-4_12_18' ~ 'Cortex[4-12-18]',
  Month_Group == 'Hippocampus-4' ~ 'Hippocampus[4]',
  Month_Group == 'Hippocampus-12' ~ 'Hippocampus[12]',
  Month_Group == 'Hippocampus-18' ~ 'Hippocampus[18]',
  TRUE ~ 'Hippocampus[4-12-18]'
))


CCX_HIP_36 <- CCX_HIP_glastly[CCX_HIP_glastly$glayCluster == '3' |
                              CCX_HIP_glastly$glayCluster == '6',]

g <- igraph::graph.data.frame(d=CCX_HIP_E,
                              directed = F)

g_mouse<- RedeR::subg(g=g, dat = CCX_HIP_36,
                      refcol = 13,
                      maincomp = F,connected = F)
g_mouse<- att.setv(g = g_mouse, from='name', to="nodeAlias")
g_mouse<- att.setv(g = g_mouse,from="Month_Group", to="nodeColor",
                   cols = brewer.pal(n = 7, name = "Accent"))
g_mouse <- att.setv(g_mouse,from='Degree',to='nodeSize')

rdp <- RedPort() 
calld(rdp)


addGraph(rdp,g_mouse,gcoord=c(10,25),isNest=TRUE, theme='tm1', zoom=30)


relax(rdp,p2=400,p5=30,ps=T)
scl <- g_mouse$legNodeColor$scale 
leg <- g_mouse$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="")

scl <- g_human$legNodeSize$scale
leg <- g_human$legNodeSize$legend
addLegend.size(rdp, sizevec=scl, labvec=leg, title="")

resetd(rdp)
