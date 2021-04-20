removeSpeciesMultilayer <- function(network,species){
  network[species,] <- 0
  network[,species] <-0
  return(network)
}

nc <- function (matrix){
  #matrix: the adjacent matrix of association network
  eig <- eigen(matrix,only.values = TRUE)
  exp.eig <- exp(eig$values)
  natural_connectivity <- log(mean(exp.eig))
  return(natural_connectivity)
}

stability <- function(net,percent,decrease){

  dataMain <- as.matrix(get.adjacency(net, attr = 'weight'))
  dataMain[dataMain != 0] <- 1
  
  dataMainOrg <- dataMain
  
  aboundance <- c(V(net)$abundance)
  names(aboundance) <- V(net)$name
  
  aboundance <- aboundance   %>%  sort(decreasing = decrease)   %>%  names()
  
  network_connetivity<-list()
  extinctionOrder<-list()
  network_connetivity[[1]] <- nc(dataMain)
  x=0
  repeat{  
    x=x+1
    if(length(aboundance)==0)    break
    extinctionOrder[[x]]<-aboundance[1]
    dataMain<-removeSpeciesMultilayer(dataMain,aboundance[1])
    
    second<-colSums(dataMain) <= colSums(dataMainOrg)*(1-percent)
    
    dataMain<-which(second) %>% names() %>% removeSpeciesMultilayer(dataMain,.)
    
    aboundance<-aboundance   %in% colnames(dataMain)[-which(second)]   %>% which(.) %>%  aboundance[.]
    network_connetivity[[x+1]]<-nc(dataMain)
  }
  rm(x)
  network_connetivity <- do.call(rbind,network_connetivity)
  extinctionOrder <- do.call(c,extinctionOrder)
  
  rownames(network_connetivity) <- c('Initial',extinctionOrder)
  
  d <- data.frame(order = rep(0:length(extinctionOrder),1),
                  extant = c(network_connetivity),
                  Percentage=percent)
  d$extant <- as.numeric(d$extant)
  d$order <- as.numeric(d$order)
  d$order <- d$order/ncol(dataMainOrg)
  return(d)
  
}

rand_stability <- function(times,net,percent){

  dataMain <- as.matrix(get.adjacency(net, attr = 'weight'))
  dataMain[dataMain != 0] <- 1

  dataMainOrg <- dataMain
  
  aboundance <- V(net)$name %>% sample()
  
  network_connetivity<-list()
  extinctionOrder<-list()
  network_connetivity[[1]] <- nc(dataMain)
  x=0
  repeat{  
    x=x+1
    if(length(aboundance)==0)    break
    extinctionOrder[[x]]<-aboundance[1]
    dataMain<-removeSpeciesMultilayer(dataMain,aboundance[1])
    second<-colSums(dataMain) <=colSums(dataMainOrg)*(1-percent)
    dataMain<-which(second) %>% names() %>% removeSpeciesMultilayer(dataMain,.)
    
    aboundance<-aboundance   %in% colnames(dataMain)[-which(second)]   %>% which(.) %>%  aboundance[.] %>% sample()
    network_connetivity[[x+1]] <- nc(dataMain)
  }
  rm(x)
  network_connetivity <- do.call(rbind,network_connetivity)
  extinctionOrder <- do.call(c,extinctionOrder)
  
  rownames(network_connetivity) <- c('Initial',extinctionOrder)
  
  d <- data.frame(order = rep(0:length(extinctionOrder),1),
                  extant = c(network_connetivity),
                  Percentage=percent,Times = times)
  d$extant <- as.numeric(d$extant)
  d$order <- as.numeric(d$order)
  d$order <- d$order/ncol(dataMainOrg)
  return(d)
  
}