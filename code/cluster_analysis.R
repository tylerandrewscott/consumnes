set.seed(24)

packs <- c('tidyverse','devtools','data.table','pbapply','stringr','igraph','bibliometrix','Matrix','Matrix.utils','parallel','criticalpath','slam','scales')
need <- packs[!packs %in% installed.packages()[,'Package']]
lapply(need,install.packages)
lapply(packs,require,character.only = T)

fls = list.files('input/wos/feb15_query2/',full.names = T)
M_list <- pblapply(fls,convert2df)

M_dt <- rbindlist(M_list,fill = T,use.names = T)
M_dt <- M_dt[M_dt$LA=='ENGLISH',]
M_dt <- M_dt[grepl('ARTICLE',M_dt$DT),]
M_dt$PY[is.na(M_dt$PY)] <- str_extract(M_dt[is.na(PY),]$EA,'[0-9]{4}')
summary(M_dt$PY)

#summary(M_dt$PY)
results <- biblioAnalysis(M_dt, sep = ";")
options(width=100)
S <- summary(object = results, k = 30, pause = FALSE)
library(htmlTable)
htmlTable(S$MostRelSources)


library(tidyr)
CR <- citations(M_dt, field = "article", sep = ";")
ref_list <- str_split(M_dt$CR,"\\s{0,};\\s{0,}")
A_citation <- cocMatrix(M_dt, Field = "CR", sep = "\\s{0,};\\s{0,}",type = 'sparse')

colnames(A_citation) <- str_remove_all(colnames(A_citation),'\\.|\\,|\\s$|^\\*')
colnames(A_citation) <- str_remove(colnames(A_citation),'\\sP[0-9]{1,}$')
rownames(A_citation) <- M_dt$SR
rownames(A_citation) <- str_remove_all(rownames(A_citation),'\\.|\\,|\\s$')
rowyears <- str_extract(rownames(A_citation),'[0-9]{4}')
colyears <- str_extract(colnames(A_citation),'[0-9]{4}')
#A_citation <- A_citation[,is.na(colyears) | colyears>=min(rowyears.na.rm = T)]
A_citation <- A_citation[,!grepl('^ORGANIZATION FOR ECONOMIC CO(-|)OPERATION|^OECD',colnames(A_citation))]

labels<-rownames(A_citation)
b<-aggregate.Matrix(A_citation,labels)
tb<-t(b)
labels2<-rownames(tb)
b2<-aggregate.Matrix(tb,labels2)
A_citation_2 <- t(b2)
A_citation_2 <- A_citation_2[,colnames(A_citation_2)!="NO TITLE CAPTURED"]
A_citation_2 <- A_citation_2[,grepl('[A-Z]',colnames(A_citation_2))]


A_tripmat_full <- as.simple_triplet_matrix(A_citation_2)
A_edgelist_full <- data.table(from = A_tripmat_full$i,to = A_tripmat_full$j)
A_edgelist_full$to <- colnames(A_citation_2)[A_edgelist_full$to]
A_edgelist_full$from <- rownames(A_citation_2)[A_edgelist_full$from]

A_tripmat <- as.simple_triplet_matrix(A_citation_2)
A_edgelist <- data.table(from = A_tripmat$i,to = A_tripmat$j)

A_edgelist$to <- colnames(A_citation_2)[A_edgelist$to]
A_edgelist$from <- rownames(A_citation_2)[A_edgelist$from]
A_edgelist <- A_edgelist[to!=from,]

while(any(paste(A_edgelist$from,A_edgelist$to) %in% paste(A_edgelist$to,A_edgelist$from))){
  for(i in which(paste(A_edgelist$from,A_edgelist$to) %in% paste(A_edgelist$to,A_edgelist$from))){
    A_edgelist <- A_edgelist[-i,]
  }
}

table(colSums(A_citation_2>0)>=5)
keep_tos <- A_edgelist[,.N,by=.(to)][N>=5,]
A_edgelist <- A_edgelist[to %in% keep_tos$to,]

A_graph <- graph_from_edgelist(as.matrix(A_edgelist),directed = T)
A_graph <- simplify(A_graph,remove.multiple = T,remove.loops = T)
V(A_graph)$Year <- str_extract(V(A_graph)$name,'[0-9]{4}')

# clu <- components(A_graph)
# connected_component <- names(clu$membership)[clu$membership==1]
# if(any(clu$membership>1)){
#   A_graph <- induced_subgraph(A_graph,v = which(V(A_graph)$name %in% connected_component))}

A_graph_ud <- as.undirected(A_graph,'collapse')

louv <- cluster_louvain(A_graph_ud)

# modularity measure
modularity(louv)
B = modularity_matrix(A_graph)

#summary(degree(A_graph,mode = 'out'))
#summary(degree(A_graph,mode = 'in'))
#summary(degree(A_graph))

V(A_graph)$louv_cluster <- louv$membership
prank <- page_rank(A_graph)
V(A_graph)$prank <- unlist(as.vector(prank$vector))

betw <- betweenness(A_graph)
V(A_graph)$betweenness <- as.vector(betw)
betw_group <- data.table(node = V(A_graph)$name,between = betw,group = V(A_graph)$louv_cluster)
betw_group[,grank:=rank(-between,ties.method = 'average'),by=.(group)]


V(A_graph)$high_between <- V(A_graph)$name %in% betw_group[grank<=4,]$node

#V(A_graph)[betw == max(betw)]
prank_max <- data.table(node = V(A_graph)$name,group = unlist(V(A_graph)$louv_cluster),pagerank = unlist(V(A_graph)$prank))
prank_max_by_group <- prank_max[,max(pagerank),by=.(group)]

max_by_group = left_join(prank_max,prank_max_by_group)[pagerank==V1,]

V(A_graph)$highest_in_group <- V(A_graph)$name %in% max_by_group$node
norder <- data.table(node = V(A_graph)$name,cluster = V(A_graph)$louv_cluster)[order(cluster),]$node

#bifab <- RBioFabric::bioFabric(inGraph = A_graph,dropNodeLabels = T)
library(ggnetwork)
library(ggthemes)

#table(V(A_graph)$louv_cluster)
#dropcluster <- which.min(table(V(A_graph)$louv_cluster))
#E(A_graph)$weight <- edge.weights(sping, A_graph)
#A_graph_2 <- induced_subgraph(A_graph,which(V(A_graph)$louv_cluster != dropcluster))

edge.weights <- function(community, network, weight.within = 100, weight.between = 1) {
  bridges <- crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights) 
}

M_dt$SR_clean <- str_remove_all(M_dt$SR_FULL,'\\,')


M_dt2 <- M_dt[M_dt$SR_clean %in% V(A_graph)$name,]
M_dt2$cluster <- V(A_graph)$louv_cluster[match(M_dt2$SR_clean,V(A_graph)$name)]

M_dt2$text <- tolower(M_dt2$TI)

#M_dt2$text <- str_replace_all(M_dt2$text,'; ',';')
#M_dt2$text <- str_replace_all(M_dt2$text,'\\s','_')
#M_dt2$text <- str_replace_all(M_dt2$text,';',' ')
library(udpipe)
library(countrycode)
cnames <- tolower(unique(countrycode::countryname_dict$country.name.en))
ud_model <- udpipe_download_model(language = "english")
ud_model <- udpipe_load_model(ud_model$file_model)
x <- udpipe_annotate(ud_model, x = tolower(M_dt2$AB))
x <- as.data.frame(x)

x$cluster <- M_dt2$cluster[as.numeric(str_extract(x$doc_id,'[0-9]{1,}'))]
#x <- x[x$upos !='PUNCT',]
#x <- x[!x$token %in% stopwords::stopwords(),]
cluster_keys <- lapply(sort(unique(M_dt2$cluster)),function(i) {
  print(i)
  nr <- 0;nmin <- 10
  while(nr == 0){
    rake <- keywords_rake(x = x[x$cluster==i,], 
                          term = "token", group = c('cluster'),
                          relevant = x$upos[x$cluster==i] %in% c("NOUN", "ADJ",'PROPN'),
                          ngram_max = 3,n_min = nmin)  %>%
      arrange(-rake) %>%  mutate(cluster = i) %>% head(.,3) 
    # filter(!keyword %in% c('network governance','governance networks','policy networks','policy network','social networks','network analysis',
    #                   'natural resource','environmental management','environmental policy','environmental decision')) %>%
    #  filter(!keyword %in% c('policy networks','policy network','social network analysis')) %>% 
    #  head(.,3) %>% 
    nr = nrow(rake)
    nmin = nmin -1
  }
  rake
})

cluster_rake_keys <- rbindlist(cluster_keys) %>% arrange(cluster) 



#lgl_lay<-igraph::layout.lgl(A_graph)
drl_layout <- layout.drl(A_graph,use.seed = T)
#opt_layout <- layout.graphopt(A_graph)
#fr_layout = layout_with_fr(A_graph)
#kk_layout = layout_with_kk(A_graph)
#mds_layout = layout_with_mds(A_graph)
#lgl_layout = layout_with_lgl(A_graph)
#eb_communities <- edge.betweenness.community(A_graph)


cores <- detectCores() * 0.8
g <- A_graph
make_new_search <- F
if(make_new_search){
  linegraph <- make_line_graph(g)
  sink_edges <- V(linegraph)[degree(linegraph, mode = "out") == 0]
  source_edges <- V(linegraph)[degree(linegraph, mode = "in") == 0]
  search_list <- pblapply(
    source_edges,
    #all paths between vertices in linegraph (i.e., edges in "real" graph)
    #... that go between two vertices without visiting any vertex more than once
    all_simple_paths,
    graph = linegraph,
    to = sink_edges,
    mode = "out",cutoff = 6,cl = 5)
  spc <- tabulate(unlist(search_list),nbins = vcount(linegraph))
  V(linegraph)$spc <- spc
  
  paths <- unlist(search_list,recursive = FALSE)
  linegraph_edgelist <- get.edgelist(linegraph)
  sch <- sch_new() %>%
    sch_add_activities(
      id = seq(vcount(linegraph)),
      name = seq(vcount(linegraph)),
      duration = V(linegraph)$spc) %>%
    sch_add_relations(
      to = as.integer(linegraph_edgelist[,2]),
      from   = as.integer(linegraph_edgelist[,1]))
  
  schplan <- sch %>% sch_plan()
  V(linegraph)$mp <- schplan$activities$critical
  #saveRDS(linegraph,'scratch/linegraph_saved.rds')
  #gr <- readRDS('graph_scratch.RDS')
  #lg <- readRDS('linegraph_scratch.rds')
  
  saveRDS(linegraph,'scratch/linegraph_saved.rds')
}else{linegraph <- readRDS('scratch/linegraph_saved.rds')}

mp <- which(V(linegraph)$mp)

main_path <- subgraph.edges(A_graph, mp, delete.vertices = TRUE)

V(A_graph)$mp <- V(A_graph)$name %in% V(main_path)$name
V(A_graph)$dead_end <- degree(A_graph,mode = 'out')==0
V(A_graph)$start_end <- degree(A_graph,mode = 'in')==0

gnet <- ggnetwork(x = A_graph,layout = drl_layout)
gnet$name2 <- str_extract(gnet$name,'^.+[0-9]{4}')

cluster_medoids = rbindlist(lapply(unique(gnet$louv_cluster),function(x) {
  temp = gnet[!duplicated(gnet$name)&gnet$louv_cluster==x,] %>% select(x,y)
  data.table(louv_cluster = x,cluster::pam(temp,k = 1)$medoids)
}))


cluster_rake_keys$keyword[cluster_rake_keys$keyword=='social network analysis']<-'SNA'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='advocacy coalition framework']<-'ACF'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='new public management']<-'NPM'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='policy network theory']<-'policy net. theory'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='network management strategies']<-'net. manage. strateg.'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='network administrative organization']<-'NAO'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='natural resource management']<-'nat. res. manage.'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='natural management strategies']<-'net. manage. strategies'
cluster_rake_keys$keyword <- str_replace(cluster_rake_keys$keyword,' -','-')
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='social-ecological systems']<-'SES'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='network governance']<-'network gov.'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='policy network']<-'policy networks'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='environmental management']<-'env. management'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='local governance']<-'local gov.'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='european union']<-'EU'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='local governments']<-'local govs'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='policy networks']<-'policy nets.'
cluster_rake_keys$keyword[cluster_rake_keys$keyword=='transnational policy networks']<-'transnat. policy nets.'

cluster_rake_keys <- cluster_rake_keys[!duplicated(paste(keyword,cluster))&!keyword%in%c('natural resource governance','governance','natural resources'),]
cluster_phrases <- cluster_rake_keys[,paste(keyword,collapse = '+'),by=.(cluster)]
setnames(cluster_phrases,'cluster','louv_cluster')
cluster_medoids <- left_join(cluster_medoids,cluster_phrases)
cluster_medoids$name2 <- cluster_medoids$V1
net_combo = full_join(gnet,cluster_medoids[order(louv_cluster),])
net_combo <- net_combo %>% rename(cluster = louv_cluster)
dropcluster = 0

label_key <- gnet[(gnet$high_between==T|gnet$mp==T)&!duplicated(gnet$name),] %>% select(name,Year,louv_cluster)
label_key$title <- M_dt$TI[match(label_key$name,str_replace_all(M_dt$SR,'\\,',''))]
names(label_key) <- c('Author/journal','Year','cluster','Title')
label_key$cluster <- cluster_phrases$V1[label_key$cluster]
label_key <- label_key[!is.na(label_key$Title),]

summary_cluster_tab <- (as.data.frame(table(V(A_graph)$louv_cluster)))
summary_cluster_tab$keywords <- cluster_medoids[order(louv_cluster),]$name2

V(A_graph)$total_citations <- M_dt2$TC[match(V(A_graph)$name,M_dt2$SR_clean)]

temp = data.table(name = V(A_graph)$name,cites = V(A_graph)$total_citations,cluster = V(A_graph)$louv_cluster)
temp[,citation_rank:=rank(-cites),by=.(cluster)]

summary_cluster_tab$top3_cited <- temp[citation_rank<=3,][order(cluster,citation_rank),][,paste(name,collapse = ','),by=.(cluster)]$V1

htmlTable(summary_cluster_tab)

htmlTable(label_key,rnames = F,file = 'output/paper_clustergraph_key.html')

A_graph_nobject <- intergraph::asNetwork(A_graph)
mmat <- network::mixingmatrix(A_graph_nobject,'louv_cluster')
mmat <- melt(mmat)
mmat <- mmat[mmat$value>0,]
mmat_edges <- as.matrix(mmat)

cluster_graph <- igraph::graph_from_edgelist(mmat_edges[,1:2],directed = T)
V(cluster_graph)$cluster <- V(cluster_graph)
E(cluster_graph)$weight <- mmat_edges[,3]
cluster_graph <- simplify(cluster_graph,remove.multiple = T,remove.loops = T)

V(cluster_graph)$paper_count <- tabulate(network::get.vertex.attribute(A_graph_nobject,'louv_cluster'))
V(cluster_graph)$keywords <- sapply(sort(unique(cluster_rake_keys$cluster)),function(x) paste(cluster_rake_keys$keyword[cluster_rake_keys$cluster==x],collapse = '+'))

cluster_net <- ggnetwork(cluster_graph,layout = as.matrix(cluster_medoids[order(louv_cluster),.(x,y)]),scale = F)

net_combo_dt <- data.table(net_combo)

net_combo_dt$name2 <- ifelse(net_combo_dt$mp & !net_combo_dt$dead_end,paste0('*',net_combo_dt$name2,'*'),net_combo_dt$name2)

net_combo_dt[((mp==T&!dead_end&!start_end)|high_between|(!is.na(V1)&!cluster%in%dropcluster))&!duplicated(name2),]

net_combo_dt$name2[!is.na(net_combo_dt$V1)] <- net_combo_dt$V1[!is.na(net_combo_dt$V1)]

tableau_cols <- tableau_color_pal(direction = 1)(max(louv$membership))
net_combo_dt$cluster <- as.factor(net_combo_dt$cluster)


high_between_entries = rbindlist(lapply(seq_along(cluster_keys),function(i){
  net_combo_dt[net_combo_dt$high_between==T&net_combo_dt$cluster==i&!duplicated(name2),]
}))
high_between_entries <- high_between_entries[,.(name2,cluster)]

M_dt$louv_cluster <- V(A_graph)$louv_cluster[match(M_dt$SR_clean,V(A_graph)$name)]
M_dt$cluster_phrase <- cluster_phrases$V1[match(M_dt$louv_cluster,cluster_phrases$louv_cluster)]

full_join(M_dt[,.N,by=.(PY,cluster_phrase)][order(-PY),][PY<2022,] %>% mutate(PY = as.numeric(PY)),
          data.table(PY = 1978:2021))

time_plot_df <- full_join(M_dt[,.N,by=.(PY,cluster_phrase,louv_cluster)][order(-PY),][PY<2022,] %>% mutate(PY = as.numeric(PY)) ,data.table(PY = 1978:2021)) %>%
  mutate(N = ifelse(!is.na(N),N,0))

time_plot_df <- time_plot_df[order(louv_cluster,PY),]

time_plot_df[,cumN:=cumsum(N),by=.(cluster_phrase)]

time_plot_df$cluster_phrase <- fct_inorder(time_plot_df$cluster_phrase)


gg_time <- ggplot(time_plot_df[!is.na(cluster_phrase),] %>% 
                    arrange(-PY),aes(x = PY,y = cumN,colour = cluster_phrase,group = cluster_phrase)) + 
  geom_point() + geom_path() +
  ggtitle('Cumulative count of articles by cluster in publication sample') + theme_bw() +
  scale_x_continuous(name = 'Publication year',breaks=seq(1978,2021,5)) +
  scale_y_continuous(name = '# articles published') +
  scale_color_manual(values = tableau_cols)+ 
  
  #  guides(colour= 'none') + 
  theme(legend.position = c(0.4,0.7),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA))
ggsave(plot = gg_time,filename = 'output/cumulative_articles_by_year.png',dpi = 500,units = 'in',width = 6, height = 4.5)

gg_time <- ggplot(time_plot_df[!is.na(cluster_phrase),] %>% arrange(-PY),aes(x = PY,y = N,colour = cluster_phrase,group = cluster_phrase)) + 
  geom_point() + geom_path() +
  ggtitle('Yearly count of articles by cluster in publication sample') + theme_bw() +
  scale_x_continuous(name = 'Publication year',breaks=seq(1978,2021,5)) +
  scale_y_continuous(name = '# articles published') +
  scale_color_manual(values = tableau_cols)+ 
  
  #  guides(colour= 'none') + 
  theme(legend.position = c(0.4,0.7),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA))

ggsave(plot = gg_time,filename = 'output/articles_by_year.png',dpi = 500,units = 'in',width = 6, height = 4.5)

(netplot <- ggplot() + 
    theme_blank() + 
    geom_edges(lwd = 0.02,alpha = 0.2,col = 'grey40',
               data = net_combo,aes(x = x,y = y,xend = xend,yend = yend))+
    scale_color_manual(values = tableau_cols)+
    # geom_nodes(data =cluster_net,aes(x = x,y =y),pch=5,size = 10)+
    geom_nodes(pch = 21,alpha = 0.3,data = net_combo_dt[!cluster%in%dropcluster & is.na(V1),],aes(col = cluster,x = x,y = y)) + 
    geom_nodelabel_repel(data = net_combo_dt[((mp==T&!dead_end&!start_end)|high_between|(!is.na(V1)&!cluster%in%dropcluster))&!duplicated(name2),],
                         aes(x = x,y = y,label=name2,size = is.na(V1),
                             col = cluster),max.overlaps = Inf,box.padding= 1.5,max.iter = 1e6,force = 1.1)+
    guides(colour = 'none',size = 'none') + 
    scale_size_manual(values = c(4,2.5)) + 
    #ggtitle('citation communities') +
    #scale_x_continuous(limits = c(.25,NA))+
    labs(caption = '*[]* = on critical path\nother labeled nodes = highest betweenness scores by cluster') +
    NULL)

ggsave(netplot,filename = 'output/paper_clustergraph.png',dpi = 500,units = 'in',height = 10,width = 10)

#net_combo_dt[!is.na(V1),][order(cluster),]
#net_cluster_combo_dt <- net_combo_dt[!is.na(V1)&!cluster%in%dropcluster,][order(-cluster),]
cluster_net <- data.table(cluster_net)[order(cluster),]

net_cluster_combo_dt <- rbind(net_combo_dt,cluster_net,use.names = T,fill = T)

M_dt2[,citation_rank:=rank(-TC,ties.method = 'first')]
M_dt2[,citation_rank_by_year:=rank(-TC,ties.method = 'first'),by=.(PY)]
M_dt2[,citation_rank_by_cluster:=rank(-TC,ties.method = 'first'),by=.(cluster)]

M_dt2$citation_per_year <- M_dt2$TC/(2022 - as.numeric(M_dt2$PY))
M_dt2$citation_per_year[M_dt2$citation_per_year==Inf] <- NA
M_dt2[,citation_per_year_rank:=rank(-citation_per_year,ties.method = 'average')]

(netplot2 <- ggplot() + 
    theme_blank() + 
    # geom_edges(lwd = 0.01,alpha = 0.2,col = 'grey40',data = net_combo,aes(x = x,y = y,xend = xend,yend = yend))+
    scale_color_manual(values = tableau_cols) + 
    geom_edges(data = net_cluster_combo_dt[!is.na(paper_count)&!is.na(weight)&weight>=mean(weight,na.rm = T),],
               aes(x = x,y = y,xend = xend,yend = yend,alpha = weight),
               arrow = arrow(length = unit(0.3, "lines"), type = "open"),
               curvature = .1,size = 1.25)+
    geom_nodes(pch = 21,alpha = 0.2,data = net_cluster_combo_dt[!cluster%in%dropcluster & is.na(V1),],aes(col = cluster,x = x,y = y)) + 
    
    geom_nodes(data = net_cluster_combo_dt[!is.na(paper_count),],
               aes(x = x,y = y,size = paper_count^2,
                   col = cluster),max.overlaps = Inf,pch = 19)+
    scale_size_binned(range = c(2,8)) +
    geom_nodelabel_repel(data = net_cluster_combo_dt[!is.na(paper_count)&!duplicated(keywords),],
                         aes(x = x,y = y,label=keywords,
                             col = cluster),max.overlaps = Inf,
                         box.padding= 1.6,max.iter = 2e6,force = 1.2)+
    guides(colour = 'none',size = 'none',alpha = 'none') + 
    ggtitle('Cluster-level citation network') + 
    labs(caption = '*node size = # papers in cluster\n**edge darkness = # citations between clusters\n***edges with # papers < mean edge value omitted')+
    NULL + theme(title = element_text(size = 14)))

ggsave(plot = netplot2,filename = 'output/cluster_graph.png',dpi = 500,units = 'in',height = 10,width = 10)


ed <- get.edgelist(main_path)
edg <- paste(ed[,1],ed[,2],sep = '|')
ag <- get.edgelist(A_graph)
eids <- which(paste(ag[,1],ag[,2],sep = '|') %in% edg)
mp_graph <- subgraph.edges(A_graph,eids,delete.vertices = T)
end_caps <- V(mp_graph)[degree(mp_graph,mode = 'out')==0]
start_caps <- V(mp_graph)[degree(mp_graph,mode = 'in')==0]

cap_ranks <- rank(-colSums(A_citation_2[,which(colnames(A_citation_2) %in% names(end_caps))]))<=3

keep_caps <- names(which(degree(mp_graph,mode = 'out')>0|V(mp_graph)$name %in% names(cap_ranks[cap_ranks==T])))

mp_graph_caps <- induced.subgraph(mp_graph,vids = which(V(mp_graph)$name %in% keep_caps))

#V(A_graph)$top20btw <- rank(-V(A_graph)$betweenness,ties.method = 'average') <= 20
names_on_path <- V(mp_graph)$name[!(V(mp_graph)$name %in% names(c(start_caps,end_caps)))]
A_edge <- igraph::get.edgelist(A_graph)

paper_on_path <- V(mp_graph)$name
paper_on_path <- paper_on_path[!(paper_on_path %in% names(end_caps) | paper_on_path %in% names(start_caps))]
path_papers <- data.table(paper = paper_on_path,
                          year = str_extract(paper_on_path,'[0-9]{4}'))
path_papers <- path_papers[order(year),]
empty_dt <- data.table(from = NULL,to = NULL)

top = seq(2)
for(i in 1:nrow(path_papers)){
  print(i)
  p <- path_papers$paper[i]
  cited_by <- as.vector(which(A_citation_2[,p]>0))  
  cited_names <- rownames(A_citation_2)[cited_by]
   top_2 <- M_dt[SR_clean %in% cited_names & !SR_clean %in% path_papers$paper & !SR_clean %in% empty_dt$from,][order(-TC),][,.(SR_clean,TC)][top,]
  print(p)
  print(top_2)
  empty_dt <- rbind(empty_dt,data.table(from = top_2$SR_clean,to = p))
}

extra_edges_dt <- empty_dt

extra_edges <- which(paste(A_edgelist$from,A_edgelist$to) %in% 
                       paste(extra_edges_dt$from,extra_edges_dt$to))

sort_set <- A_edgelist[extra_edges,] %>% 
  mutate(enames = paste(from,to,sep ='|'))
sort_set <- sort_set[!(sort_set$from %in% V(mp_graph_caps)$name & sort_set$to %in% V(mp_graph_caps)$name),]
extra_graph <-  subgraph.edges(A_graph,E(A_graph)[sort_set$enames])
extra_graph <- delete_vertex_attr(extra_graph,'mp')


#combo_graph <- union(mp_graph_caps,extra_graph)
combo_graph <- union(mp_graph_caps,extra_graph)
backbone_net <- ggnetwork(combo_graph) 
backbone_net$name_short <- str_extract(backbone_net$name,'^.+[0-9]{4}')
backbone_net$year <- str_extract(backbone_net$name_short,'[0-9]{4}')
backbone_net$first_name <- str_extract(backbone_net$name_short,'^[A-Z]+')
backbone_net$mp_1[is.na(backbone_net$mp_1)]<-0
backbone_net$mp_2[is.na(backbone_net$mp_2)]<-0
backbone_net <- data.table(backbone_net)

backbone_net$louv_cluster <- ifelse(!is.na(backbone_net$louv_cluster_1),backbone_net$louv_cluster_1,backbone_net$louv_cluster_2)
backbone_net$cluster_keywords <- cluster_phrases$V1[match(backbone_net$louv_cluster,cluster_phrases$louv_cluster)]


backbone_net$cluster_keywords <- gsub('policy networks|policy network','policy nets.',backbone_net$cluster_keywords)
backbone_net$mp[is.na(backbone_net$mp)]<-0
backbone_net$louv_cluster <- as.factor(backbone_net$louv_cluster)
backbone_net$cluster_keywords <- as.factor(backbone_net$cluster_keywords)
backbone_net$cluster_keywords <- fct_reorder(backbone_net$cluster_keywords,as.numeric(backbone_net$louv_cluster))
color_sub <- unique(as.numeric(as.character(backbone_net$louv_cluster)))
color_sub  <- sort(color_sub )
backbone_net$on_path <- backbone_net$name %in% path_papers$paper
backbone_net$good_label <- 
  paste0(ifelse(backbone_net$name %in% backbone_net$name[backbone_net$mp==1] & backbone_net$dead_end_1!=T & backbone_net$start_end_1!=T,'*',''),backbone_net$first_name,' ',backbone_net$year,ifelse(backbone_net$name %in% backbone_net$name[backbone_net$mp==1] & backbone_net$dead_end_1!=T & backbone_net$start_end_1!=T,'*',''))


#backbone_net$cluster_color <- tableau_cols[backbone_net$louv_cluster]
library(ggnewscale)
backbone_net2 <- backbone_net
backbone_net2$x <-backbone_net2$x*-1
backbone_net2$y <-backbone_net2$y*-1
backbone_net2$xend <-backbone_net2$xend *-1
backbone_net2$yend <-backbone_net2$yend *-1
gpath <- ggplot(data = backbone_net2,
                 aes(x = x,y = y, yend = yend,xend = xend,colour = cluster_keywords,
                     label = good_label)) + 
    geom_nodes(aes(size = as.factor(mp))) + 
    geom_edges(aes(size = as.factor(mp),linetype = as.factor(mp)),
               arrow = arrow(length = unit(0.3, "lines"), type = "closed")) + 
    scale_color_manual(name = 'community keywords',
                       values = tableau_cols[color_sub])+
    scale_linetype_manual(values = c(2,1),name = 'edge type',labels=c('cites path','critical path')) + 
    scale_size_manual(values = c(0.75,2.5),name = 'edge type',labels=c('cites path','critical path'))+
  theme_blank() + 
    theme(title = element_text(size = 20),text = element_text(family = 'Times'),
          legend.background = element_rect(fill = alpha('white',0)),
          legend.spacing = unit(.2,'cm'),legend.position = c(0.75,0.25)) + 
    ggtitle('Critical path + most cited papers citing the path') +
  guides(linetype = 'none',
         color = guide_legend(override.aes = list(size = 5))) +  
  new_scale("size")+  
  geom_nodelabel_repel(aes(size = as.factor(on_path+0)),
                       min.segment.length = 0.2,
                       max.overlaps = Inf,force = 12) +  
  scale_size_manual(values = c(3,4))+
  guides(size = 'none')+
  labs(caption = '*paper* = on critical path')
ggsave(filename = 'output/critical_path.png',plot = gpath,width = 10,height = 10,units = 'in',dpi = 500)

grep("WILLIAMSON",colnames(A_citation_2),value = T)
M_dt2[,.(TC,SR_FULL,TI)][order(-TC),][1:10]
label_key <- backbone_net
label_key <- label_key[,.(name_short,name,cluster_keywords,mp)]
label_key$title <- M_dt$TI[match(label_key$name,str_replace_all(M_dt$SR,'\\,',''))]
htmlTable(label_key,rnames = F,file = 'output/critical_path_key.html')



gpath
plot(subgraph.edges(gr,eids = which(V(lg)$mp)))


E(gr)[which(V(lg)$mp)]

library(tidytext)
library(tm)
title_words <- M_dt2 %>%
  unnest_tokens(word, text) %>%
  count(cluster, word, sort = TRUE) %>% 
  filter(!word %in% tidytext::stop_words$word)

total_words <- title_words %>% 
  group_by(cluster) %>% 
  summarize(total = sum(n))

title_words <- left_join(title_words, total_words)

freq_by_rank <- title_words %>% 
  group_by(cluster) %>% 
  mutate(rank = row_number(), 
         `term frequency` = n/total) %>%
  ungroup()
freq_by_rank

book_tf_idf <- title_words %>%
  bind_tf_idf(word, cluster, n)

book_tf_idf %>%
  select(-total) %>%
  arrange(desc(tf_idf))
library(udpipe)




install.packages('udpipe')
tdf <- data.table(book_tf_idf)
cluster_top <- tdf[,cluster_rank:=rank(-tf_idf,ties.method = 'first'),by=.(cluster)]
cluster_top <- cluster_top[cluster_rank<=2,]


cluster_top[cluster%in%c(8)&cluster_rank<=4,]$word
lubs <- grep('LUBELL',V(A_graph)$name)
V(A_graph)$louv_cluster[lubs]
V(cluster_graph)$keywords[c(1,8)]


















g <- A_graph
linegraph <- make_line_graph(g)
sink_edges <- V(linegraph)[degree(linegraph, mode = "out") == 0]
source_edges <- V(linegraph)[degree(linegraph, mode = "in") == 0]

search_list2 <- pblapply(
  source_edges,
  #all paths between vertices in linegraph (i.e., edges in "real" graph)
  #... that go between two vertices without visiting any vertex more than once
  all_simple_paths,
  graph = linegraph,
  to = sink_edges,
  mode = "out",cutoff = 5,cl = 7)


spc <- tabulate(unlist(search_list2),nbins = vcount(linegraph))
V(linegraph)$spc <- spc


paths <- unlist(search_list2,recursive = FALSE)
linegraph_edgelist <- get.edgelist(linegraph)
sch <- sch_new() %>%
  sch_add_activities(
    id = seq(vcount(linegraph)),
    name = seq(vcount(linegraph)),
    duration = V(linegraph)$spc) %>%
  sch_add_relations(
    to = as.integer(linegraph_edgelist[,2]),
    from   = as.integer(linegraph_edgelist[,1]))

schplan <- sch %>% sch_plan()
V(linegraph)$mp <- schplan$activities$critical
saveRDS(A_graph,'graph_scratch.RDS')
saveRDS(linegraph,'linegraph_scratch.rds')


#mp <- which(schplan$activities$critical)
mp_graph <- subgraph.edges(A_graph, eids = which(V(linegraph)$mp==T), delete.vertices = TRUE)
mp_net <- ggnetwork(mp_graph)
mp_net$name_short <- str_extract(mp_net$name,'^.+[0-9]{4}')

mp_net$year <- str_extract(mp_net$name_short,'[0-9]{4}')
mp_net$first_name <- str_extract(mp_net$name_short,'^[A-Z]+')

gpath <- ggplot(data = mp_net,
                aes(x = x,y = y, yend = yend,xend = xend,
                    label = paste(first_name,year))) + 
  geom_nodes() + geom_edges(arrow = arrow(length = unit(0.3, "lines"), type = "closed")) + 
  geom_nodelabel_repel(size = 3) + 
  theme_blank() +
  ggtitle('Critical path')

ggsave(plot = gpath,filename = 'output/critical_path.png',units = 'in',width = 10,height = 5)

E(A_graph)$mp <- (V(linegraph)$mp==T)+0
gnet <- ggnetwork(x = A_graph,layout = nice_layout)
gnet$name2 <- str_extract(gnet$name,'^.+[0-9]{4}')

netplot <- ggplot() + 
  theme_blank() + 
  geom_edges(data = gnet[gnet$mp==F,],aes(x = x,y = y,xend = xend,yend = yend),alpha=.2,lwd = 0.02)+
  geom_edges(data = gnet[gnet$mp==T,],aes(x = x,y = y,xend = xend,yend = yend),alpha=.8,lwd = .8,col = 'orange',arrow = arrow(length = unit(0.3, "lines"), type = "closed"))+
  geom_nodes(data = gnet[gnet$mp==T,],aes(x = x,y = y),col = 'orange')


V(A_graph)$name[1:10]


yrs <- c(1990,2000,2010,Inf)

#gdt <- data.table(gnet)
grobs <- lapply(yrs,function(year){
  A_graph_nobject <- intergraph::asNetwork(A_graph)
  node_year <- network::get.vertex.attribute(A_graph_nobject,'Year')
  node_year[is.na(node_year)]<-2021
  
  subnet <- network::get.inducedSubgraph(A_graph_nobject,v = which(node_year<year))
  mmat <- network::mixingmatrix(subnet,'louv_cluster')
  mmat_edges <- data.table::melt(mmat)
  mmat_edges <- mmat_edges[mmat_edges$value>0,]
  
  cluster_graph <- igraph::graph_from_edgelist(as.matrix(mmat_edges[,1:2]),directed = T)
  E(cluster_graph)$weight <- mmat_edges$value
  cluster_graph <- simplify(cluster_graph,remove.multiple = T,remove.loops = T)
  V(cluster_graph)$paper_count <- tabulate(network::get.vertex.attribute(A_graph_nobject,'louv_cluster'))
  V(cluster_graph)$cluster <- 1:length(V(cluster_graph))
  V(cluster_graph)$keywords <- sapply(sort(unique(cluster_top$cluster)),function(x) paste(cluster_top$word[cluster_top$cluster==x],collapse = '+'))
  cluster_net <- ggnetwork(cluster_graph,layout = circle)
  ggplot(cluster_net) + 
    geom_edges(aes(x = x,y = y,xend = xend,yend = yend,
                   size = sqrt(weight),alpha = weight),curvature = .1,
               arrow = arrow(length = unit(0.3, "lines"), type = "open")) + 
    geom_nodes(aes(x = x,y = y,size = paper_count,color = as.factor(cluster))) + 
    geom_nodelabel_repel(aes(label= keywords,x = x,y = y),size = 1.5)+
    theme_blank() + scale_color_tableau() +
    ggtitle(paste0('Citations between clusters through ',ifelse(year==Inf,'present',year)))+
    guides(color = 'none',size= 'none',alpha= 'none')
})

library(gridExtra)
cluster_by_time <- grid.arrange(grobs[[1]],grobs[[2]],grobs[[3]],grobs[[4]],ncol = 2)
ggsave(plot = cluster_by_time,filename = 'output/cluster_by_time.png',dpi = 500,units = 'in',height = 10,width = 10)




cluster_year <- data.table(V(A_graph)$Year,V(A_graph)$louv_cluster)
cluster_year$V1[is.na(cluster_year$V1)]<-2021
cluster_year <- cluster_year[order(V2,V1),]
cluster_year$n <- 1

cluster_year[,cumcount:=cumsum(n),by=.(V2)]
cluster_year <- left_join(cluster_year,data.table(V2 = 1:10,keywords = sapply(sort(unique(cluster_top$cluster)),function(x) paste(cluster_top$word[cluster_top$cluster==x],collapse = '+'))))
cluster_year$V1 <- as.numeric(cluster_year$V1)
ggplot(cluster_year[V1>=1980,]) + 
  scale_y_continuous(name = '# papers')+
  scale_x_continuous(name = 'year')+
  geom_path(aes(x = V1,y = cumcount,group = V2,col = as.factor(keywords))) + 
  scale_color_tableau(name = 'key terms') + theme_bw()

labs(caption = '*node size proportional to # papers\n**edge size and shading proportional to edge count\n***top 2 terms by tf-idf by cluster')


data.table(V(cluster_graph)$paper_count,V(cluster_graph)$keywords)
S$MostRelSources
V(cluster_graph)$paper_count
plot(cluster_graph)
as.matrix(mmat[-11,])

class(as.matrix(mmat[-11,]))

tdf[cluster_rank<=2,]
book_tf_idf 
geom_edges(lwd = 0.05,col = 'grey40',aes(col = as.factor(mp),alpha = mp))+
  scale_color_tableau(palette = 'Tableau 10')+
  
  
  geom_nodes(pch = 21,alpha = 0.3,
             aes(col = as.factor(louv_cluster))) + 
  #geom_nodelabel_repel(data = gnet[gnet$high_between==T&!duplicated(gnet$name),],size = 1.5,aes(x = x,y = y,label=name2,col = as.factor(louv_cluster)),max.overlaps = 100)+
  guides(colour = 'none') + ggtitle('Louvain clusters') +
  labs(caption = 'labeled nodes = top 3 betweenness scores by cluster')
netplot


str(paths)
path_lengths <- unlist(lapply(paths, function (x) sum(x$spc)))
vertex_attr(linegraph, "main_path") <- 0
vertex_attr(
  linegraph,
  "main_path",
  paths[[which(path_lengths == max(path_lengths))[[1]]]]) <- 1
V(linegraph)$main_path



imc<-cluster_infomap(A_graph)

table(as.vector(membership(imc)))


table(V(A_graph)$louv_cluster)




cmod <- fastgreedy.community(A_graph)

V(A_graph)



sort(unique(test$membership))

cluster_lov

pranks <- page_rank(A_graph)

summary(pranks$vector)


V(A_graph)

sort(unique(test$membership))

prank <- page_rank(cite_graph)

V(cite_graph)[prank$vector==max(prank$vector)]

ergm::summary_formula(net~edges)







cl <- cluster_edge_betweenness(cite_graph,weights = NULL,
                               bridges = F,merges = F,
                               edge.betweenness = F,directed = T)




##CURENT PARTITION
# run louvain with edge weights 
louvain_partition <- igraph::cluster_louvain(cite_graph, weights = E(cite_graph)$weight) 
# assign communities to graph 
cite_graph$community <- louvain_partition$membership 
# see how many communities there are 
unique(cite_graph$community) 


V(cite_graph)
communities <- data.frame() 

A_tripmat_p5 <- as.simple_triplet_matrix(A_citation_3)
A_edgelist_p5 <- data.table(from = A_tripmat_p5$i,to = A_tripmat_p5$j)





for (i in unique(cite_graph$community)) { 
  i = 10
  # create subgraphs for each community 
  subgraph <- induced_subgraph(cite_graph, v = which(cite_graph$community == i)) 
  # get size of each subgraph 
  size <- igraph::gorder(subgraph) 
  # get betweenness centrality 
  btwn <- igraph::betweenness(subgraph) 
  
  V(subgraph)[11]
  communities <- communities %>% 
    dplyr::bind_rows(data.frame(
      community = i, 
      n_characters = size, 
      most_important = names(which(btwn == max(btwn))) 
    ) 
    ) 
} 
knitr::kable(
  communities %>% 
    dplyr::select(community, n_characters, most_important)
)

cores <- detectCores() * 0.8
g <- A_graph

linegraph <- make_line_graph(g)
sink_edges <- V(linegraph)[degree(linegraph, mode = "out") == 0]
source_edges <- V(linegraph)[degree(linegraph, mode = "in") == 0]



search_list2 <- pblapply(
  source_edges,
  #all paths between vertices in linegraph (i.e., edges in "real" graph)
  #... that go between two vertices without visiting any vertex more than once
  all_simple_paths,
  graph = linegraph,
  to = sink_edges,
  mode = "out",cutoff = 6,cl = 30)

bads <- which(sapply(search_list2,class)!='list')
bads
search_list2[bads] <- pblapply(
  source_edges[bads],
  #all paths between vertices in linegraph (i.e., edges in "real" graph)
  #... that go between two vertices without visiting any vertex more than once
  all_simple_paths,
  graph = linegraph,
  to = sink_edges,
  mode = "out",cutoff = 6,cl = 1)

t2 <- pblapply(
  source_edges[bads,
               #all paths between vertices in linegraph (i.e., edges in "real" graph)
               #... that go between two vertices without visiting any vertex more than once
               all_simple_paths,
               graph = linegraph,
               to = sink_edges,
               mode = "out",cutoff = 6,cl = 30)

t2
all_simple_paths(from = source_edges[1203],to = sink_edges,graph = linegraph,mode = 'out',cutoff = -1)

table(sapply(sapply(search_list2,class),is.null))



str(search_list)
spc <- tabulate(unlist(search_list),nbins = vcount(linegraph))
V(linegraph)$spc <- spc
#Select edges constituting the main path 
#(e.g. starting from source nodes iteratively select 
# the edge with the highest weight 
# or select the path with the overall highest weight)  

paths <- unlist(search_list,recursive = FALSE)

linegraph_edgelist <- get.edgelist(linegraph)

sch <- sch_new() %>%
  sch_add_activities(
    id = seq(vcount(linegraph)),
    name = seq(vcount(linegraph)),
    duration = V(linegraph)$spc) %>%
  sch_add_relations(
    to = as.integer(linegraph_edgelist[,2]),
    from   = as.integer(linegraph_edgelist[,1]))

schplan <- sch %>% sch_plan()
V(linegraph)$mp <- schplan$activities$critical
saveRDS(linegraph,'linegraph_5mp.rds')
saveRDS(A_graph,'agraph.rds')


mp <- which(schplan$activities$critical)


mp_graph <- subgraph.edges(A_graph, mp, delete.vertices = TRUE)

mp_graph


plot(mp_graph)
grep('LUBELL',get.vertex.attribute(mp_graph,'name'),value = T)

plot(mp_graph)
table(str_extract(get.vertex.attribute(mp_graph,'name'),
                  '[0-9]{4}'))


str(V(mp_graph))
as.character(V(mp_graph))
grep('BIX',as.character(V(mp_graph)))



str_split(E(A_graph)[mp][1],'->')
str_split(E(A_graph)[mp]
          
          grep("ERNST",colnames(A_citation_2),value = T)
          plot(schplan)
          sch_validate(sch)
          criticalpath::sch_critical_activities(schplan)
          linegraph_edgelist %>% head()
          
          sch_validate(sch)
          path_lengths <- unlist(pblapply(paths, function (x) 
            sum(V(linegraph)$spc[x]),cl = 4))
          
          vertex_attr(linegraph, "main_path") <- 0
          vertex_attr(
            linegraph,
            "main_path",
            paths[[which(path_lengths == max(path_lengths))[[1]]]]) <- 1
          V(linegraph)$main_path
          
          
          which(V(linegraph)$main_path==1]

igraph::plot.igraph()

}


ms <- main_search(A_graph,cores = cores)

10375 * 12 

120,900

29800/100000
29800/ (120900 * 1.03)


124527/12
library(igraph)
#https://stackoverflow.com/questions/67792685/main-path-analysis-in-citation-network-using-igraph-in-r
cit_g <- graph_from_incidence_matrix(A_citation_2,directed = T)

length(igraph::vertices(cit_g)[[1]])

length(unique(c(rownames(A_citation_2),colnames(A_citation_2))))
is.bipartite(cit_g)
vcount(cit_g)

dim(A_citation_2)


x <- matrix(c(1, 0, 0, 2), nrow = 2)
s <- as.simple_triplet_matrix(x)
s



rowSums()


sink_edges
###
#The method proceeds in two steps:

#Compute an edge weight (e.g. Search Path Count [SPC] 
#counts the number of paths from source nodes to sink nodes that go 
#through an edge)



Select edges constituting the main path (e.g. starting from source nodes iteratively select the edge with the highest weight or select the path with the overall highest weight)



test <- main_search(cit_g)

igraph::betweenness(cit_g)

plot(cit_g,
     layout=layout_with_sugiyama(cit_g,
                                 layers = V(cit_g)$name)$layout)






test <- Matrix.utils::aggregate.Matrix(A_citation,groups = colnames(A_citation),fun = 'sum')
class(A_citation)
Matrix.utils::aggregate.Matrix()


y <- colnames(A_citation)
ymat <- model.matrix(~ 0 + y)
dim(ymat)

colnames(ymat) <- 1:10
head(ymat)


x1 <- m %*% ymat
x2 <- d %*% ymat
all.equal(as.matrix(x1), x2)



#https://stackoverflow.com/questions/36778166/how-to-combine-columns-with-identical-names-in-a-large-sparse-matrix
OP2 <- function(x) {
  nms <- colnames(x)
  uniquenms <- unique(nms)
  # build the sparseMatrix again: x's with same index values are automatically
  # added together, keeping in mind that indexes stored from 0 but built from 1
  sparseMatrix(i = x@i + 1, 
               j = match(nms, uniquenms)[x@j + 1],
               x = x@x,
               dimnames = list(rownames(x), uniquenms),
               giveCsparse = FALSE)
}

A_citation_2 <- OP2(A_citation)




grep('DOI',colnames(A_citation),value = T)


str(slam::as.simple_triplet_matrix(A_citation))

rownames(A_citation)[duplicated(rownames(A_citation))]


nms <- colnames(A_citation)
tf_mat_bind <- sapply(unique(nms), function(i)rowSums(A_citation[, nms==i, drop=FALSE]))




grep('BIXLER',rownames(A_citation),value = T)



dim(A_tripmatrix)



rownames(A_citation) %in% colnames(A_citation)

rownames(A_citation) 

grep("POLICY STUD J",colnames(A_citation), value = T)

class(A_citation)

grep("WASSERMAN",colnames(A_citation),value = T)
table(duplicated(colnames(A_citation)))
table()

dim(A_citation)

colnames(A_citation) <- str_remove(colnames(A_citation),'\\s{1,}$')
which(duplicated(colnames(A_citation)))


colnames(A_citation)[70]

dimnames(A_citation)[[2]][1:10]

table(A_citation)

grep('HILEMAN.*LUBELL',M_dt$AU,value = F)


?cocMatrix

M_dt[510,]$CR

?A_citation
CR$Source[1]
CR$Cited[1]
head(CR)
grep("LUBELL",dimnames(A_citation)[[2]])
M_dt[1,]$TI
M_dt[1,]$CR

str(CR)
CR[[1]][2]


dim(A_citation)
A_author <- cocMatrix(M_dt, Field = "AU", sep = ";")

A_citation <- cocMatrix(M_dt, Field = "CR", sep = ".  ")

dimnames(A_citation)[[1]][1:10]


library(stringr)
coln <- str_remove(dimnames(A_citation)[[2]],'\\s{1,}$')
M_dt$ROWNAME<-str_remove_all(M_dt$SR_FULL,'\\,')
coln[coln %in% M_dt$ROWNAME]

M_dt[ROWNAME=="LUBELL M 2020 SOC NATUR RESOUR",]

grep("LUBELL",coln,value = T)
grep("LUBELL",M_dt$ROWNAME,value = T)
head(coln)
table(coln %in% str_remove_all(M_dt$SR_FULL,'\\,'))

coln[1:10]
sort(M_dt$SR_FULL)[1:10]
M_dt$SR[1]
M_dt$SR_FULL[1]

table(M_dt$SR_FULL %in% dimnames(A_citation)[[2]])

any(dimnames(A_citation)[[1]] %in% dimnames(A_citation)[[2]])
dimnames(A_citation)[[2]][1:4]

dna_g <- graph_from_data_frame(dna_edges, directed=T)
plot(dna_g,
     layout=layout_with_sugiyama(dna_g,
                                 layers = V(dna_g)$name)$layout)


dna_g <- graph_from_data_frame(dna_edges, directed=T)
plot(dna_g,
     layout=layout_with_sugiyama(dna_g,
                                 layers = V(dna_g)$name)$layout)


M_dt[SO=='IEEE ACCESS']
cocite_network <- biblioNetwork(M_dt, analysis = "co-citation", network = "references", sep = ";")

cocite_network

options(width=130)
histResults <- histNetwork(M_dt, min.citations = 4, sep = ";")


# Plot a historical co-citation network
net <- histPlot(histResults, n=300, size = 2, labelsize=3)

CS <- conceptualStructure(M_dt,field="ID", method="CA", minDegree=4, clust=5, stemming=FALSE, labelsize=10, documents=10)

plot(CS)
str(cocite_network)


# Plot the network
net=networkPlot(NetMatrix, n = dim(NetMatrix)[1], Title = "Country Collaboration", type = "circle", size=TRUE, remove.multiple=FALSE,labelsize=0.7,cluster="none")



plot(A_citation)

M_dt$CR
head(M_dt)
M_dt$ID
M_dt$CR
dim(M_dt)
sum(sapply(M_list,nrow))

library(criticalpath)

install.packages('criticalpath')

