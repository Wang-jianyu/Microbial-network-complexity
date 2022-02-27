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