unodf=function(A=matrix, selfloop=FALSE){
  
  if(nrow(A)!=ncol(A)){ stop("The matrix is non-square (m-by-n for which m-n). Please input a square matrix.")}
  
  if(selfloop==FALSE & sum(diag(A))>0){
    warning("Please note that the network contains self-loops and these were discarded. To consider them, please set 'selfloop=TRUE'.") 
    diag(A)=0
    }
  if(selfloop==TRUE & sum(diag(A))==0){
    warning("Please note that 'selfloop=TRUE' but this network does not contain self-loops. If self-loops are possible in your system, ignore this message and proceed; otherwise please set 'selfloop=FALSE' and run it again.")
    }

  m=dim(A)[1]
  n=dim(A)[2]
  
  cmarg=apply(A,2,sum) 
  nest.m.c=matrix(0,n,n)
  
  for (i in 1:n){
    
    for (j in 1:n){
      
      if (cmarg[i]<cmarg[j]){ 
        
        den=cmarg[i] 
        if(selfloop==FALSE){ if(A[j,i]==1){den=den-1} }  			
        if (den>0){
          sum.r=A[,i] + A[,j]
          nest=sum(sum.r==2)/den
          nest.m.c[i,j]=nest 
        }	
      }
    }
  }
  
  rmarg=apply(A,1,sum) 
  nest.m.r=matrix(0,m,m)
  
  for (i in 1:m){
    
    for (j in 1:m){
      
      if (rmarg[i]<rmarg[j]){ 
        
        den=rmarg[i]
        if(selfloop==FALSE){ if(A[i,j]==1){den=den-1} }
        
        if (den>0){
          sum.r=A[i,] + A[j,]
          nest=sum(sum.r==2)/den
          nest.m.r[i,j]=nest
        }	
      }
    }
  }
  
  aux1=(sum(cmarg==0)*(n-1))
  aux2=(sum(rmarg==0)*(m-1))
  
  network.nest.c=(2*sum(nest.m.c, na.rm=TRUE))/((n*(n-1)-aux1))
  network.nest.r=(2*sum(nest.m.r,na.rm=TRUE))/((m*(m-1)-aux2)) 
  
  result=data.frame(network.nest.c,network.nest.r)
  names(result)=c("UNODFc","UNODFr")
  return(result)
}

