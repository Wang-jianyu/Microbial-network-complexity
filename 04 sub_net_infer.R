# There is an outlier that we use the mean of the remaining 
# 19 sample-level network attributes to supplement.

require(igraph)
require(tidyverse)
rm(list=ls())
source('function.R')
graph.name <- list.files("trim_graph//", pattern = '.graphml')
otutb.name <- list.files("trim_otu//", pattern = '.csv')

net_topo <- function(net){ 
  require(igraph)
  topo <- c(
    Node.number = length(V(net)),
    
    Edge.number = length(E(net)),
    
    Positive.edge = sum(E(net)$dir>0), 
    
    Negative.edge = sum(E(net)$dir<0),
    
    Average.degree = mean(igraph::degree(net)),
    
    Diameter = diameter(net, directed = FALSE, unconnected = TRUE),
    
    Clustering.coefficient = transitivity(net),
    
    Average.path.length = mean_distance(net),
    
    Modularity = modularity(net,membership(cluster_fast_greedy(net))),
    
    Centralization.closeness = centralization.closeness(net)$centralization,
    
    Centralization.betweenness = centralization.betweenness(net)$centralization,
    
    Edge.connectivity = edge_connectivity(net),
    
    Average.neighborhood = mean(ego_size(net)),
    
    Connectance = edge_density(net,loops=FALSE),
    
    Nestedness = c(unodf(as.matrix(get.adjacency(net,attr = 'weight')))[1],use.names = F)
    
  )
  return(topo)
}




  ot <- paste("trim_otu/",otutb.name,sep = "") %>% 
    read.csv(header = TRUE,row.names = 1)
  net <- read.graph(paste("trim_graph/",graph.name,
                          sep = ""),format = "graphml")
  topo_trim <- NULL
  for (k in 1:20) {
    
    sub_net <- induced_subgraph(net, which(V(net)$name %in% 
                                             as.character(rownames(ot)[which(ot[,k]!=0)])),
                                impl = 'create_from_scratch')
    if(length(E(sub_net)) == 0) next
    
    topo<-net_topo(sub_net)
    
    # write.graph(sub_net,paste0("trim_graph/sub_graph/",
    #                            gsub(".csv",paste0('_sample_',k,".graphml"),otutb.name)),
    #             format = "graphml")
    topo<-c(gsub(".csv",paste0('_sample_',k),otutb.name),topo)
    topo_trim <- rbind(topo_trim,topo)
    
  }

write.csv(topo_trim,'trim_graph/sub_graph/topo_trim_sub.csv')