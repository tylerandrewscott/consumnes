

packs <- c('devtools','data.table','pbapply','stringr','igraph','bibliometrix','Matrix','Matrix.utils','parallel','criticalpath')
need <- packs[!packs %in% installed.packages()[,'Package']]
lapply(need,install.packages)
lapply(packs,require,character.only = T)

fls <- list.files('input/wos/feb14_query_2/',full.names = T)
M_list <- pblapply(fls,convert2df)

M_dt <- rbindlist(M_list,fill = T,use.names = T)

results <- biblioAnalysis(M_dt, sep = ";")
options(width=100)
S <- summary(object = results, k = 30, pause = FALSE)

CR <- citations(M_dt, field = "article", sep = ";")
library(slam)
ref_list <- str_split(M_dt$CR,"\\s{0,};\\s{0,}")
A_citation <- cocMatrix(M_dt, Field = "CR", sep = "\\s{0,};\\s{0,}",type = 'sparse')

colnames(A_citation) <- str_remove_all(colnames(A_citation),'\\.|\\,|\\s$|^\\*')
colnames(A_citation) <- str_remove(colnames(A_citation),'\\sP[0-9]{1,}$')
rownames(A_citation) <- M_dt$SR
rownames(A_citation) <- str_remove_all(rownames(A_citation),'\\.|\\,|\\s$')

labels<-rownames(A_citation)
b<-aggregate.Matrix(A_citation,labels)
tb<-t(b)
labels2<-rownames(tb)
b2<-aggregate.Matrix(tb,labels2)
A_citation_2 <- t(b2)
A_citation_2 <- A_citation_2[,colSums(A_citation_2)>=5&colnames(A_citation_2)!="NO TITLE CAPTURED",]


A_tripmat <- as.simple_triplet_matrix(A_citation_2)


A_edgelist <- data.table(from = A_tripmat$i,to = A_tripmat$j)

A_edgelist$to <- colnames(A_citation_2)[A_edgelist$to]
A_edgelist$from <- rownames(A_citation_2)[A_edgelist$from]

A_graph <- graph_from_edgelist(as.matrix(A_edgelist),directed = T)

A_graph <- simplify(A_graph,remove.multiple = T,remove.loops = T)

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

A_graph <- readRDS('graph_scratch.RDS')
linegraph <- readRDS('linegraph_scratch.rds')

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

?main_search
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

set.seed(1)
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

