rm(list = ls())
require(vegan)
require(tidyverse)
require(igraph)
require(reshape2)

#Create new folders
dir.create('trim_otu/');dir.create('trim_graph/');
dir.create('trim_otu/perm');dir.create('trim_otu/pmat');    
dir.create('trim_otu/boot_mean');dir.create('trim_otu/boot_sd');
dir.create('trim_graph/sub_graph/')

#Elevation only for 3755m
otutb.name <- list.files("trim_otu//", pattern = '.csv')

####Jaccard dissimilarity####
###permutation procedure
jac_p_matrix <- function(otutable=paste("trim_otu/",otutb.name,sep = "")){
  perm <- function(x){return(sample(x))}
  ot <- read.csv(otutable,header = TRUE,row.names = 1)
  jac.mean.perm <- matrix(0, nrow=nrow(ot), ncol=nrow(ot))
  for(l in 1:1000){
    ot.perm <- apply(ot, MARGIN = 1, FUN = perm) 
    jac.perm <- vegdist(t(ot.perm),method = "jaccard")
    jac.mean.perm  <- jac.mean.perm + (1-as.matrix(jac.perm)) 
    # print(x = paste("l =",l))
  }
  jac.mean.perm <- round(jac.mean.perm/1000,4)
  return(jac.mean.perm)
}
###bootstrap procedure
jac_boot_matrix <- function(otutable=paste("trim_otu/",otutb.name,sep = "")){
  ot <- read.csv(otutable,header = TRUE,row.names = 1)
  jac.sd.boot <- matrix(0, nrow=nrow(ot), ncol=nrow(ot))
  jac.mean.boot <- matrix(0, nrow=nrow(ot), ncol=nrow(ot))

  for(j in 1:1000){
    jac.boot <- vegdist(ot[,sample(ncol(ot),replace = TRUE)],method = "jaccard")
    ###Prevent NA value
    if(length(which(is.na(as.matrix(jac.boot))))>0){
      repeat{ 
        jac.boot <- vegdist(ot[,sample(ncol(ot),replace = TRUE)],method = "jaccard")
        if(length(which(is.na(as.matrix(jac.boot))))==0) break }
    }
    jac.mean.boot  <- jac.mean.boot + (1-as.matrix(jac.boot))
    # print(x = paste("j =",j))
  }
  for(k in 1:1000){
    jac.boot <- vegdist(ot[,sample(ncol(ot),replace = TRUE)],method = "jaccard")
    ###Prevent NA value
    if(length(which(is.na(as.matrix(jac.boot))))>0){
      repeat{ 
        jac.boot <- vegdist(ot[,sample(ncol(ot),replace = TRUE)],method = "jaccard")
        if(length(which(is.na(as.matrix(jac.boot))))==0) break }
    }
    
    jac.sd.boot  <- jac.sd.boot+((1-as.matrix(jac.boot))- jac.mean.boot/1000)^2
    # print(x = paste("k =",k))
  }

  jac.mean.boot <- round(jac.mean.boot/1000,4)
  jac.sd.boot <- round(sqrt(jac.sd.boot/1000),4)

  return(list(mean = jac.mean.boot,
              sd = jac.sd.boot))
}



  ot.trim.jac.boot <- jac_boot_matrix(otutable=paste("trim_otu/",otutb.name,sep = ""))
  ot.trim.jac.perm <- jac_p_matrix(paste("trim_otu/",otutb.name,sep = ""))
  write.csv(ot.trim.jac.perm, paste("trim_otu/perm/jac_",otutb.name,sep = ""))
  write.csv(ot.trim.jac.boot[[1]],paste("trim_otu/boot_mean/jac_",otutb.name,sep = ""))
  write.csv(ot.trim.jac.boot[[2]],paste("trim_otu/boot_sd/jac_",otutb.name,sep = ""))




####Spearman correlation####
###permutation procedure
sp_p_matrix <- function(otutable=paste("trim_otu/",otutb.name,sep = "")){
  perm <- function(x){return(sample(x))}
  ot <- read.csv(otutable,header = TRUE, row.names = 1)
  sp.mean.perm <- matrix(0, nrow=nrow(ot), ncol=nrow(ot))
  for(l in 1:1000){
    ot.perm <- apply(ot, MARGIN = 1, FUN = perm)
    sp.perm <- cor(ot.perm,method = "spearman")
    sp.mean.perm  <- sp.mean.perm + as.matrix(sp.perm)  
    # print(x = paste("l =",l))
  }
  sp.mean.perm <- round(sp.mean.perm/1000,4)
  return(sp.mean.perm)
}

###bootstrap procedure
sp_boot_matrix <- function(otutable=paste("trim_otu/",otutb.name,sep = "")){
  ot <- t(read.csv(otutable,header = TRUE,row.names = 1))
  sp.sd.boot <- sp.mean.boot <- matrix(0, nrow=ncol(ot), ncol=ncol(ot))
  for(j in 1:1000){
    sp.boot <- cor(ot[sample(nrow(ot),replace = TRUE),],
                     method = "spearman")
    sp.mean.boot  <- sp.mean.boot + as.matrix(sp.boot)  
    # print(x = paste("j =",j))
  }
  for(k in 1:1000){
    sp.boot <- cor(ot[sample(nrow(ot),replace = TRUE),],
                     method = "spearman")
    sp.sd.boot  <- sp.sd.boot+(as.matrix(sp.boot) - sp.mean.boot/1000)^2
    # print(x = paste("k =",k))
  }

  sp.mean.boot <- round(sp.mean.boot/1000,4)
  sp.sd.boot <- round(sqrt(sp.sd.boot/1000),4)
  return(list(mean = sp.mean.boot,
              sd = sp.sd.boot))
}


  ot.trim.sp.boot <- sp_boot_matrix(paste("trim_otu/",otutb.name,sep = ""))
  ot.trim.sp.perm <- sp_p_matrix(paste("trim_otu/",otutb.name,sep = ""))
  ot.trim.sp.perm[is.na(ot.trim.sp.perm)] <- 0
  ot.trim.sp.boot[[1]][is.na(ot.trim.sp.boot[[1]])] <- 0
  ot.trim.sp.boot[[2]][is.na(ot.trim.sp.boot[[2]])] <- 0
  write.csv(ot.trim.sp.perm, paste("trim_otu/perm/sp_",otutb.name,sep = ""))
  write.csv(ot.trim.sp.boot[[1]],paste("trim_otu/boot_mean/sp_",otutb.name,sep = ""))
  write.csv(ot.trim.sp.boot[[2]],paste("trim_otu/boot_sd/sp_",otutb.name,sep = ""))

### estimate p values####
# The P value was then obtained as the probability of 
# the null value under a Gaussian curve fitted to the 
# mean and standard deviation of the bootstrap distribution.

  sp.mean.perm <- read.csv(paste("trim_otu/perm/sp_",otutb.name,sep=""),row.names = 1)
  sp.mean.boot <- read.csv(paste("trim_otu/boot_mean/sp_",otutb.name,sep=""),row.names = 1)
  sp.sd.boot <- read.csv(paste("trim_otu/boot_sd/sp_",otutb.name,sep=""),row.names = 1)

  cor.p <- matrix(0, nrow=nrow(sp.mean.boot), ncol=ncol(sp.mean.boot))
  for(i in 1:nrow(sp.mean.boot)) {
    for (j in 1:nrow(sp.mean.boot)){        
      p <- ks.test(x = abs(sp.mean.perm[i,j]),
                   "pnorm",
                   mean = abs(sp.mean.boot[i,j]),
                   sd = sp.sd.boot[i,j])
      cor.p[i,j] <- p$p.value
    }
  }
  write.csv(cor.p,paste("trim_otu/pmat/sp_",otutb.name,sep=""))



  jac.mean.perm <- read.csv(paste("trim_otu/perm/jac_",otutb.name,sep=""),row.names = 1)
  jac.mean.boot <- read.csv(paste("trim_otu/boot_mean/jac_",otutb.name,sep=""),row.names = 1)
  jac.sd.boot <- read.csv(paste("trim_otu/boot_sd/jac_",otutb.name,sep=""),row.names = 1)
  
  cor.p <- matrix(0, nrow=nrow(jac.mean.perm), ncol=ncol(jac.mean.perm))
  for(i in 1:nrow(jac.mean.perm)) {
    for (j in 1:nrow(jac.mean.perm)){        
      p <- ks.test(x = abs(jac.mean.perm[i,j]),
                   "pnorm",
                   mean = abs(jac.mean.boot[i,j]),
                   sd = jac.sd.boot[i,j])
      cor.p[i,j] <- p$p.value
    }
  }
  write.csv(cor.p,paste("trim_otu/pmat/jac_",otutb.name,sep=""))