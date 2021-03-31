library('RColorBrewer')
library('stringr')
library('ggplot2')
library('readr')
library('dplyr')
library('igraph')
library('RedeR')
### Analysis of DEG-DTU and CEMitool results (Ora plots)

H_DEG_DTU <- readRDS('CEMtool/results/human/H_DEG_DTU.rds')

F_DEG_DTU <- readRDS('CEMtool/results/mouse_fad/F_DEG_DTU.rds')

T_DEG_DTU <- readRDS('CEMtool/results/mouse_tau/T_DEG_DTU.rds')


# Human


human_cem <- list(modules= read_tsv('CEMtool/results/Tables/Human/module.tsv'),
                  ora = read_tsv('CEMtool/results/Tables/Human/ora.tsv'),
                  ppi = read_tsv('CEMtool/results/Tables/Human/interactions.tsv'))

human_modules <- human_cem$modules

H_DEG_DTU <- H_DEG_DTU[H_DEG_DTU$genes %in% human_modules$genes,]
human_modules <- merge(human_modules,H_DEG_DTU,by='genes')

human_modules <- human_modules[human_modules$modules != 'Not.Correlated',]

ALZ_DEG <- human_modules[human_modules$Group == 'ALZ' & 
                         human_modules$Class == 'DEG',]

ALZ_DTU <- human_modules[human_modules$Group == 'ALZ' & 
                           human_modules$Class == 'DTU',]

human_ppi <- human_cem$ppi
human_ppi <- human_ppi[,2:3]

g <- igraph::graph.data.frame(d=human_ppi,
                              directed = F)

g_ALZ_DEG <- RedeR::subg(g=g, 
  dat = ALZ_DEG,refcol = 1,
                       maincomp = F,connected = T)

g_ALZ_DEG <-att.setv(g=g_ALZ_DEG,from='Age_Group', to='nodeColor',
                   cols = brewer.pal(n = 7, name = "Accent"))

g_ALZ_DEG <-att.setv(g=g_ALZ_DEG,from='modules',to = 'nodeShape')


g_ALZ_DTU <- RedeR::subg(g=g, 
                         dat = ALZ_DTU,refcol = 1,
                         maincomp = F,connected = T)

g_ALZ_DTU <-att.setv(g=g_ALZ_DTU,from='Age_Group', to='nodeColor',
                     cols = brewer.pal(n = 7, name = "Accent"))

g_ALZ_DTU <-att.setv(g=g_ALZ_DTU,from='modules',to = 'nodeShape')

### RedeR

rdp <- RedPort() 
calld(rdp)


addGraph(rdp, g_ALZ_DTU,gscale=30,zoom=30)

relax(rdp,50,400)
scl <- g_ALZ_DEG$legNodeColor$scale 
leg <- g_ALZ_DEG$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, 
                title="Age_Groups")


resetd(rdp)
