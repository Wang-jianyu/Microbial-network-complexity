require(igraph)
require(tidyverse)
rm(list=ls())
if(grepl('Windows',sessionInfo()$running)) 
  source('unodf.R',encoding = 'utf-8') else require(UNODF)


graph.name <- list.files("network_graph//", pattern = '.graphml')
otutb.name <- list.files("network_otu//", pattern = '.csv')

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




  ot <- paste("network_otu/",otutb.name,sep = "") %>% 
    read.csv(header = TRUE,row.names = 1)
  net <- read.graph(paste("network_graph/",graph.name,
                          sep = ""),format = "graphml")
  topo_network <- NULL
  for (k in 1:20) {
    
    sub_net <- induced_subgraph(net, which(V(net)$name %in% 
                                             as.character(rownames(ot)[which(ot[,k]!=0)])),
                                impl = 'create_from_scratch')
    if(length(E(sub_net)) == 0) next
    
    topo<-net_topo(sub_net)
    
    # write.graph(sub_net,paste0("network_graph/sub_graph/",
    #    gsub(".csv",paste0('_sample_',k,".graphml"),otutb.name)),format = "graphml")
    topo<-c(gsub(".csv",paste0('_sample_',k),otutb.name),topo)
    topo_network <- rbind(topo_network,topo)
    
  }

write.csv(topo_network,'network_graph/sub_graph/topo_sub_network_sub.csv')