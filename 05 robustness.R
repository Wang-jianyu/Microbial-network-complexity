rm(list=ls())
require(igraph)
require(tidyverse)
dir.create('robustness/')
source('function.R')
otutb.name <- list.files("network_otu//", pattern = '.csv')



net_list <- read.graph(paste("network_graph/",
              gsub("csv","graphml",otutb.name), sep = ""),format = "graphml")

for (percent in c(0.25,0.5,0.75)) {
  
  result <- net_list %>% robustness(percent = percent,decrease=T)
  area <- result %>%  sum((result$extant)*(result$order[2])) 
  names(area) <- paste0('Most-to-least_stablity',percent)
  print(area)
  
  result1 <- net_list %>% robustness(percent = percent,decrease=F)
  area1 <- result1 %>%  sum((result1$extant)*(result1$order[2])) 
  names(area1) <- paste0('Least-to-most_stablity',percent)
  print(area1)

  write.csv(result ,paste0('robustness/',percent,'_Most-to-least_xy.csv'))
  write.csv(result1 ,paste0('robustness/',percent,'_Least-to-most_xy.csv'))
  
}
###################Random############################


otutb.name <- list.files("network_otu/", pattern = '.csv')
area1 <- result1 <- NULL


require(igraph)
for (percent in c(0.25,0.5,0.75)){

    net <- read.graph(paste("network_graph/",
                            gsub("csv","graphml",otutb.name),
                            sep = ""),format = "graphml")
    
    result<-lapply(1:100,function(times)rand_robustness(times=times,net=net,percent = percent)) %>% 
      lapply(.,function(.) data.frame('Name' = otutb.name,.))
    
    area <- lapply(result,function(a) sum((a$extant))*(a$order[2])) %>% do.call(cbind,.)
    
    result <- do.call(rbind,result)
    area1 <- rbind(area1,area)
    print(paste('percent = ',percent,date()))
    result1<- rbind(result1,result)
  


 
}
write.csv(area1,paste0('robustness/rand_robustness_area.csv'),row.names = F)
write.csv(result1,paste0('robustness/rand_xy_robustness.csv'))
