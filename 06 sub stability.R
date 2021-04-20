rm(list=ls())

require(igraph)
graph.name <- list.files("network_graph/sub_graph", pattern = '.graphml')
otutb.name <- list.files("network_otu//", pattern = '.csv')
dir.create('sub_stability/')
source('function.R')
area1 <- area2 <- list()


for (percent in c(0.25,0.5,0.75)) {
  


  
  area1 <-lapply(1:20,function(k) {read.graph(paste0("network_graph/sub_graph/",
                             gsub(".csv",paste0('_sample_',k,".graphml"),otutb.name)),
                      format = "graphml")}) %>% 
      lapply(.,function(net) stability(net=net,percent = percent,decrease=T)) %>% 
      lapply(.,function(x) sum((x$extant)*(x$order[2]))) %>% do.call(cbind,.)
    
  area2 <-lapply(1:20,function(k) {read.graph(paste0("network_graph/sub_graph/",
                                                     gsub(".csv",paste0('_sample_',k,".graphml"),otutb.name)),
                                              format = "graphml")}) %>% 
      lapply(.,function(net) stability(net=net,percent = percent,decrease=F)) %>% 
      lapply(.,function(x) sum(x$extant)*(x$order[2])) %>% do.call(cbind,.)
    
  print(paste('percent = ',percent,date()))
       

  area3<-lapply(area1, function(x)  data.frame('Simulation scenarios'='Most-to-least abundant taxa removal','Stability'=x))
  area4<-lapply(area2, function(x)  data.frame('Simulation scenarios'='Least-to-most abundant taxa removal','Stability'=x))
  cbind (data.frame(do.call(rbind,area3)) ,
         data.frame(do.call(rbind,area4)) ) %>% write.csv(paste0('sub_stability/sub_',percent,'_area.csv'))



}

###############################################
rm(list=ls())

otutb.name <- list.files("network_otu/", pattern = '.csv')
require(igraph)
source('function.R')
for (percent in c(0.25,0.5,0.75)) {
  result1 <- NULL
  sub_area <- list()

  for (k in 1:20) {
    
    sub_net <- read.graph(paste0("network_graph/sub_graph/",
                                 gsub(".csv",paste0('_sample_',k,".graphml"),otutb.name)),
                          format = "graphml")
    #For simplicity, only random 10 times
    sub_area[[k]] <- lapply(1:10,function(times)rand_stability(times=times,net=sub_net,percent = percent)) %>%
      lapply(.,function(x) data.frame('Times' = x[1,'Times'],'Area' =sum(x$extant)*(x$order[2]))) %>% 
      lapply(.,function(x) data.frame('Name' = otutb.name,'Sample' = k,'Percent' = percent,x)) %>% 
      do.call(rbind,.)
    print(paste('Percent =',percent,"k =",k,date()))
    
  }
  result <- do.call(rbind,sub_area)
  result1 <- rbind(result1 , result)

write.csv(result1,paste0('sub_stability/rand_',percent,'_stability_area.csv'))

}
