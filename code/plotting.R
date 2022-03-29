library(igraph)
library(networkD3)

g <- readRDS('agraph.rds')
m5 <- readRDS('linegraph_5mp.rds')
m6 <- readRDS('linegraph_6mp.rds')

keep = which(V(m6)$mp)
g_sub <- subgraph.edges(g,keep)

library(tidyverse)
library(ggnetwork)


V(g_sub)$year <- str_extract(get.vertex.attribute(g_sub,'name'),'[0-9]{4}')
V(g_sub)$year[is.na(V(g_sub)$year)] <- 2021
V(g_sub)$order <- rank(V(g_sub)$year,ties.method = 'first')

get.vertex.attribute(g_sub,'name')
net <- ggnetwork(g_sub)
net$name[grepl('ANGST',net$name)] <- "ANGST M 2021 PUBLIC ADMIN REV"
net$name2 <- str_extract(net$name,'.*[0-9]{4}')
library(ggrepel)

ggplot(data = net,
  aes(x = x, y = y,  xend = xend, yend = yend))+
  geom_nodes()+
  geom_edges(lwd = 1.5) +
  geom_text_repel(data = net[!duplicated(net$name),],col = 'grey50',nudge_x = .1,
                  aes(label = name2,x = x,y = y),size = 2,min.segment.length = .25) +
  theme_blank()


?geom_edges
fortify.igraph(g_sub)
ggplot()



V(g_sub)$year

networkD3::
plot(g_sub)
subgraph(g,)
table(V(m5)$mp)
table(V(m6)$mp)
