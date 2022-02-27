#### network property ####
rm(list=ls())
require(igraph)
if(grepl('Windows',sessionInfo()$running)) 
  source('unodf.R',encoding = 'utf-8') else require(UNODF)

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

vn_property <- function(net,otutable){
  ot <- read.csv(otutable,header = TRUE,row.names = 1)
  
  require(igraph)
  
  V(net)$node.degree = degree(net)
  
  V(net)$node.norm.degree = degree(net, normalized=TRUE) 
  
  V(net)$neighborhood = ego_size(net)
  
  V(net)$node.clustering.coefficient = transitivity(net,type="local")
  
  V(net)$node.betweenness = betweenness(net)
  
  V(net)$node.centralization.betweenness = betweenness(net,normalized = TRUE)
  
  V(net)$node.closeness = closeness(net) 
  
  V(net)$node.centralization.closeness = closeness(net,normalized = TRUE)
  
  V(net)$abundance = rowSums(ot)[which(names(rowSums(ot)) %in% V(net)$name)]
  
  fc<-cluster_fast_greedy(net)  
  
  modularity = modularity(net,membership(fc))
  
  comps = membership(fc)
  
  modulenumber <- c(1:max(comps))
  
  V(net)$module = modulenumber[comps]
  
  E(net)$pn<-sign(E(net)$dir)
  
  V(net)$bf<-substr(V(net)$name,1,1)
  
  E(net)$edge_type <- apply(strtrim(get.edgelist(net),1),1,function(x) paste0(x[1],x[2]))
  
  return(net)
  
}

##### Calculate network property #######


  net <- read.graph(paste("network_graph/",
                          gsub("csv","graphml",otutb.name),
                          sep = ""),format = "graphml") %>% 
        vn_property(.,otutable=paste("network_otu/",otutb.name,sep = ""))


  
  write.graph(net,file = paste("network_graph/",
                               gsub("csv","graphml",otutb.name),
                               sep = ""),format = "graphml")

  topo <- net_topo(net)
  write.csv(topo,'network_graph/topo.csv')
  