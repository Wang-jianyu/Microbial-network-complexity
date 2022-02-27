#' Unipartite Nestedness metric
#' @description Unipartite Nestedness metric, \code{UNODF}, for unipartite graphs
#' @param mat Binary matrix depicting interactions between elements of a biological system
#' @param selfloop Boolean. \code{selfloop=TRUE} if there are self-loops in the network (i.e. a node can be connected to itself; there are non-zeroed element in the diagonal of the input matrix \code{mat}); \code{selfloop=FALSE} otherwise.
#' @usage unodf(x)
#' @return Returns a data frame with nestedness degree for the columns \code{(UNODFc)} and rows \code{(UNODFr)} of the matrix
#' @details The Unipartite Nestedness metric, UNODF, is the Nestedness metric based on Overlap and Decreasing Fill for two-mode networks (NODF, Almeida-Neto et al. 2008) generalized for one-mode networks (see also  the metric proposed by Lee et al. 2012; see Cantor et al. for details).
#' @author Mathias Pires
#' @examples # binary square matrix
#' mat = rbind(c(1,1,1,1,1),
#'             c(1,1,1,0,1),
#'             c(1,1,0,1,0),
#'             c(1,0,1,0,0),
#'             c(1,1,0,0,0))
#'             colnames(mat) = paste("element", seq(1:5), sep="")
#'             
#'# Nestedness in rows and columns
#'unodf(mat, selfloop=TRUE)
#' @references Almeida-Neto M, Guimaraes P, Guimaraes PR Jr, Loyola RD, Ulrich W. 2008 A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement. Oikos 117:1227–1239. doi:10.1111/j.2008.0030-1299.16644.x.
#' 
#' Lee DS, Maeng SE, Lee JW. 2012. Scaling of nestedness in complex networks. J Korean Phys Soc 60:648-656. doi: 10.3938/jkps.60.648
#' 
#' Cantor M, Pires MM, Marquitti FDM, Raimundo RLG, Sebastian-Gonzalez E, Coltri P, Perez I, Barneche D, Brandt DYC, Nunes K, Daura-Jorge FG, Floeter SR, Meyer D, Guimaraes PR Jr. 2017 Nestedness across biological scales. PLoS ONE 12(2): e0171691. doi:10.1371/journal.pone.0171691
#' @export
unodf=function(A=matrix, selfloop=FALSE){
  
  if(nrow(A)!=ncol(A)){ stop("The matrix is non-square (m-by-n for which m ≠ n). Please input a square matrix.")}
  
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








#' Null model 2: Randomize binary matrices based on rows and columns totals
#' @description This null model creates random binary matrices by resorting the 1's among the matrix cells according to marginal totals of rows and columns.
#' @param mat Binary matrix
#' @param iter Number of iterations
#' @param symmetric Boolean, \code{default=FALSE}. If \code{TRUE}, it assumes the input matrix \code{mat} is symmetrical and so the random ones will also be (with diagonal=0 and lower triangle = upper triangle)
#' @usage null2(mat, iter)
#' @return mat.t array of randomized matrices
#' @details The function creates random binary matrices by randomly resorting the 1's among the matrix cells according to marginal totals of rows and columns (Bascompte et al. 2003, Pires et al. 2011). See details in Cantor et al.
#' @references Bascompte J, Jordano P, Melian CJ and Olesen JM (2003) The nested assembly of plant-animal mutualistic networks. PNAS, 100(16):9383-9387
#' 
#' Cantor M, Pires MM, Marquitti FDM, Raimundo RLG, Sebastian-Gonzalez E, Coltri P, Perez I, Barneche D, Brandt DYC, Nunes K, Daura-Jorge FG, Floeter SR, Meyer D, Guimaraes PR Jr. Nestedness across biological scales.
#' @author Mathias Pires, Mauricio Cantor
#' @examples # Creating a binary matrix
#' mat =  rbind(c(1,1,0,0,0),
#'               c(1,0,1,0,1),
#'               c(0,1,0,1,1),
#'               c(0,0,0,1,1),
#'               c(0,1,1,0,0))
#' 
#' null2(mat, 100) 
#' @export
null2 <- function(mat, iter=1000, symmetric=FALSE){
  
  if(symmetric != isSymmetric(mat)){ stop("Please double-check if the matrix is symmetrical; either input a symmetrical matrix or change the 'symmetric' argument")}

  nR <- nrow(mat)   ; nC <- ncol(mat)
  mR <- rowSums(mat); mC <- colSums(mat)
  
  # if matrix is empty, skip the loop and return an empty array of equal size
  if(sum(mat)!=0){ 
    
    #generating a matrix with probabilities for each cell
    tR   <- rep(mR, each=nC)
    tC   <- rep(mC, nR)
    prob <- matrix(((tR/nC)+(tC/nR))/2, nR, nC, byrow=TRUE)
    
    #filling theoretical matrices
    mat.t=array(0,c(nR,nC,iter))
    
    for(s in 1:iter){
      rand=matrix(runif(nR*nC),nR,nC)
      aux=mat.t[,,s]
      #fill empty matrix whenever empirical probabilities are higher than random ones
      aux[which(rand<prob)]=1
      
      # keeping symmetry, if desired
      if(symmetric==TRUE){
        tri <- lower.tri(aux)
        aux[tri]<- t(aux)[tri]
        diag(aux)=0
      }
      
      #Making sure random and empirical matrices are of equal size
      #empirical and random number of interactions
      fill <- sum(mat) ; randfill <- sum(aux)
      
      #if random interactions < empirical, select random rows and cols and place remaining interactions
      if(randfill < fill){
        remaining = fill-randfill
        emptycells = which(aux==0, arr.ind = T)                        #empty cells in the rand matrix
        
        if(symmetric==TRUE){                                           # for symmetric matrices:
          emptycells = emptycells[-(which(emptycells[,1]==emptycells[,2])),]# do not consider the diagonal, which is zeroed
          remaining=remaining/2                                        # we need half of the empty matrices (lower tri)
        }
        
        pick = sample(nrow(emptycells), size=remaining, replace=FALSE) #pick empty row/col at random
        for(k in 1:length(pick)){
          aux[ emptycells[pick[k],][1], emptycells[pick[k],][2] ] = 1  #fill them in the rand matrix
          if(symmetric==TRUE){
            aux[ emptycells[pick[k],][2], emptycells[pick[k],][1] ] = 1#both triangles for symmetric matrix
          }
        }
      }
      
      #if random interactions > empirical, select random rows and cols and remove remaining interactions
      if(randfill > fill){
        excess = randfill-fill
        filledcells = which(aux==1, arr.ind = T)                        #filled cells in the rand matrix
        
        if(symmetric==TRUE){                                            #for symmetric matrices:
          excess=excess/2                                               #we need half of filled matrice (lower tri)
        }
        pick2 = sample(nrow(filledcells), size=excess, replace=FALSE)   #pick filled row/col at random
        for(k in 1:length(pick2)){
          aux[ filledcells[pick2[k],][1], filledcells[pick2[k],][2] ] = 0         #remove them in the rand matrix
          if(symmetric==TRUE){
            aux[ filledcells[pick2[k],][2], filledcells[pick2[k],][1] ] = 0       #both triangles for symmetric matrix
          }
        }
      }
      
      #storing matrices
      mat.t[,,s]=aux
    }  
    return(mat.t)
    
    # if input matrix was empty in the first place, it returns an empty array
  } else {
    mat.t=array(0,c(nR,nC,iter))
    mat.t
  }
}




#' @title Importing multiple files into a list
#' @description The function reads all files of a given extension in a given directory path into a list of objects. 
#' @param dir Full path to the target directory
#' @param pattern A string with characters that give the pattern of the file name (default: \code{*.txt})
#' @param length Number of characters for the name of each file in the list
#' @details This function will apply the function \code{read.table()} to read the files in the directory, usually text files of the \code{.txt} extension.
#' @return A list of objects; each object is a file.
#' @author Mauricio Cantor
#' @export
import <- function(dir="data", pattern="*.txt", length=15){
  
  filename <- list.files(dir, pattern, full.names=TRUE)
  auxlist  <- lapply(filename, read.table)
  #  names(auxlist) <- substr(filename, nchar(dir)+2 , length)
  aux <- gsub(paste(dir,'/', sep=""), "\\1 ", filename)
  aux <- gsub(".txt", "\\1 ", aux)
  names(auxlist) <- gsub("^\\s+|\\s+$", "", aux)
  
  auxlist
}


#' @title Binarizing quantitative matrices
#' @description Transforms quantitative network adjacency matrices into binary matrices using different cut-offs to define an interaction
#' @param data Quatitative matrix
#' @param cut vector giving the levels of cut-off. Default is \code{[0.1, 0.9]}
#' @return A list with binary matrices using different cut-offs to define a presence (i.e. 1 in the new matrix)
#' @author Mauricio Cantor
#' @export
cutoff <- function(data, cut = seq(0,0.9, length.out=10)){
  
  auxlist <- vector("list", length(cut))
  
  for(i in 1:length(cut)){
    auxmat <- data
    if(cut[i]==0){ auxmat[which(auxmat > cut[i])] = 1
    } else {    auxmat[which(auxmat >= cut[i])] = 1 }
    auxmat[which(auxmat <  cut[i])] = 0
    auxlist[[i]] <- auxmat
  }
  auxlist
}





#'@title From network edge list to adjacency matrix
#'@description This function transforms a network edge list into its square adjacency matrix couterpart
#'@param edgelist A 2-column matrix with node in, node out; or a a 3-column matrix or data frame with node in, node out and weight of interaction
#'@param type \code{binary} returns a binary matrix and \code{weighted} returns a weighted matrix 
#'@param symmetric Boolean. Is the link weight a symmetric measure (i.e. if A interacts with B, so B interacts with A with the same strength) AND does the edge list contains only one copy of the interactions (say A-B, but not B-A)? For the output matrix be also symmetric (squared and with lower triangle values same as upper triangle values), we need to double the edgelist. \code{symmetric=TRUE} does it.
#'@param numeric Boolean. Are the node labels numeric and sequential (e.g. 1,2,3...)? PLEASE NOTE: the function works better for sequential nodes (see examples) and you may have problems if the node labels are alphanumerical and non-sequential. In this case, please double check if the output makes sense.
#'@details An edge list is a two- or three-column matrix representing the nodes and the edges (links) of a network. The first two columns represent interacting nodes; the third column gives the strength of interaction for weighted networks and is absent in binary networks. Some functions require the full adjacency matrix as input, so this the \code{el_mat()} function transforms the edge list into a square matrix.
#'@return The adjacency matrix of a given network
#'@examples #Creating a binary network edge list
#'elb =  rbind(c(1,2),
#'             c(1,3),
#'             c(1,4),
#'             c(2,3),
#'             c(3,4))
#'
#'# Transforming it into a binary and symmetric adjacency matrix
#'el_mat(elb, type="binary", symmetric=TRUE, numeric=TRUE)
#'
#'# Creating a weighted edge list
#'elw =  rbind(c(1,2,0.2),
#'             c(1,3,0.3),
#'             c(1,4,0.4),
#'             c(2,3,0.5),
#'             c(3,4,0.6))
#'
#'# Transforming it into a weighted edge list
#'el_mat(elw, type="weighted", symmetric=TRUE, numeric=TRUE)
#'@author Mauricio Cantor
#'@export
el_mat <- function(edgelist, type='weighted', symmetric=FALSE, numeric=FALSE){
  
  if(symmetric==TRUE){
    aux <- edgelist
    aux[,1] <- edgelist[,2]
    aux[,2] <- edgelist[,1]
    edgelist <- rbind(edgelist, aux)
  }
  
  if(numeric==TRUE){
    edgelist=as.matrix(edgelist)
    mat=matrix(0, max(edgelist[,1]),max(edgelist[,2]))
    } else {
    mat = matrix(0, length(unique(edgelist[,1])), length(unique(edgelist[,2])))
  }
  mat.w=mat
  
  if(type=="binary"){ 
    for (k in 1:nrow(edgelist)) {
      aux1=edgelist[k,1]
      aux2=edgelist[k,2]
      mat[aux1,aux2]=1
    }
    return(mat)
  
  } else {
    
    for (k in 1:nrow(edgelist)) {
      aux1=edgelist[k,1]
      aux2=edgelist[k,2]
      mat.w[aux1,aux2]=edgelist[k,3]
    }
    return(mat.w)
  }
  
}






#'@title From network adjacency matrix to edge list
#'@description This function transforms a network square adjacency matrix into its edge list counterpart
#'@param mat matrix of interactions
#'@param type \code{binary} returns a binary edgelist and \code{weighted} returns a weighted edgelist
#'@details The  adjacecy matrix represents the nodes and interactions (links, edges) of a network. An edge list is a two- or three-column matrix representing the nodes and the edges (links) of a network. The first two columns represent interacting nodes; the third column gives the strength of interaction for weighted networks and is absent in binary networks. Some functions require the edge list as input, so this the \code{mat_el()} function transforms the adjacency matrix into the edge list.
#'@seealso \code{el_mat()}
#'@return A two-column matrix with the nodes of a network if \code{type="binary"}; a three-column matrix with the nodes of a network if \code{type="weighted"}.
#'@examples # Creating a binary network adjacency matrix
#'matb =  rbind(c(1,1,0,0,0,1,1),
#'              c(1,0,1,0,1,0,1),
#'              c(0,1,0,1,1,1,1),
#'              c(0,0,0,1,1,0,1),
#'              c(0,1,1,0,0,1,0))
#'
#'# Transforming it into a binary edge list
#'mat_el(matb, type="binary")
#'
#'# Creating a weighted matrix
#'matw =  rbind(c(0.1,0.9,0,0,0,0.3,1),
#'              c(0.6,0,0.7,0,0.5,0,0.3),
#'              c(0,0.2,0,0.4,0.4,0.4,0.2),
#'              c(0,0,0,0.6,0.7,0,0.1),
#'              c(0,0.2,0.6,0,0,0.8,0))
#'
#'# Transforming it into a weighted edge list
#'mat_el(matw, type="weighted")
#'@author Mauricio Cantor 
#'@export
mat_el <- function(mat, type='weighted'){
  
  el <- reshape2::melt(mat)
  el <- el[which(el[,3]!=0),]
  el <- as.matrix(el)
  colnames(el) = rownames(el) <- NULL
  
  if(type=="binary"){
    return(el[,1:2])
  } else {
    return(el)
  }
}




#' @title Removing zeroed columns and/or zeroed rows in matrices
#' @description The function takes a binary matrix and remove rows and/or columns which contains only zeroes
#' @param mat Binary matrix
#' @param what \code{row}: remove zeroed rows; \code{col}: remove zeroed columns; \code{both}: remove zeroed columns and rows 
#' @return A submatrix with no zeroed columns and/or rows
#' @examples
#' # Binary matrix with last row and last column zeroed
#'matb =  rbind(c(1,1,0,0,0,1,1,0),
#'              c(1,0,1,0,1,0,1,0),
#'              c(0,1,0,1,1,1,1,0),
#'              c(0,0,0,1,1,0,1,0),
#'              c(0,1,1,0,0,1,0,0),
#'              c(0,0,0,0,0,0,0,0))
#' # Removing zeroed row
#' rmzero(matb, what="row")
#' # Removing zeroed column
#' rmzero(matb, what="col")
#' # Removing zeroed row and column
#' rmzero(matb, what="both")
#' @author Mauricio Cantor
#' @export
rmzero <- function(mat, what="both"){
  if(what %in% c("both", "row", "col") == FALSE) stop("choose 'both', 'row' or 'col'")
  if(what=="both"){mat = mat[rowSums(abs(mat))>0 & rowSums(abs(mat))>0, colSums(abs(mat))>0 & colSums(abs(mat))>0]}
  if(what=="col"){mat = mat[, colSums(abs(mat))>0 & colSums(abs(mat))>0]}
  if(what=="row"){mat = mat[rowSums(abs(mat))>0 & rowSums(abs(mat))>0, ]}
  mat
}




#' @title Printing results summary
#' @description This is an internal function that organizes the output of the analysis in the script \code{analysis.R}
#' @param emp List with empirical \code{unodf()} results for each cut-off
#' @param rand List with random results for each cut-off
#' @param iter Number of iterations used to create random matrix through \code{null2()}
#' @param cuts Number of cut-offs used. Usually 10 [0 : 0.9]
#' @param save Boolean. Would you like to export a txt file with the results?
#' @param ... Save details" other arguments to be passed to \code{write.table()}
#' @author Mauricio Cantor
#' @export
printresult <- function(emp = uniind[[1]], rand = unirandind[[1]], iter=iter, cuts=10, save=FALSE, ...){
  
  result <- matrix(NA, cuts, 9)
  colnames(result) <- c("cutoff", "Nc", "P1", "2.5", "97.5", "Nr", "P2", "2.5", "97.5")
  
  result[,1] <- seq(0.0, 0.9, length.out=cuts)
  for(i in 1:cuts){
    result[i,2] = as.numeric(emp[i,1])
    result[i,6] = as.numeric(emp[i,2])
    result[i,3] = sum(rand[,i     ] > as.numeric(emp[i,1]))/iter
    result[i,7] = sum(rand[,i+cuts] > as.numeric(emp[i,2]))/iter
    result[i,4:5] = quantile(rand[,  i]   , probs=c(0.025, 0.975), type=2, na.rm=T)
    result[i,8:9] = quantile(rand[,i+cuts], probs=c(0.025, 0.975), type=2, na.rm=T)
  }
  
  if(save==TRUE){ write.table(result, ...)}
  result
}



#' @title Plotting binary and quantitative matrices
#' @description The function plots a given binary matrix with one colour for the matrix cells with presences and a another color for the cells with absences, or set a color code to plot quantitative matrices (adapted from a function of the package \code{phaget4})
#' @param x matrix
#' @param ... additional parameters to be passe to the \code{plot()} function
#' @param scale. Boolean, TRUE if you want to print a color code scale
#' @return a plot of a matrix with a colour code
#' @details original function can be retrieved at 'source("http://www.phaget4.org/R/myImagePlot.R")'
#' @author Chris Seidel
#' @references citation(package = "phaget4")
#' @export
plotmat <- function(x, scale=T, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  if(scale==TRUE){
    layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  } #else { 
    #layout(matrix(data=c(1), nrow=1, ncol=1), widths=c(4), heights=c(1))
  #}

  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  if(scale==TRUE){
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
  }
  
}


#' @title Plotting binary matrices
#' @description The function plots a given binary matrix, with one colour for the matrix cells with presences and a another color for the cells with absences
#' @param x matrix
#' @param colpres Color of cells with presence
#' @param colabs Color of cells with absence
#' @param ... other graphical paramenters to be passed to the function \code{image()}.
#' @details see  more graphical parameters in \code{help(image)}.
#' @seealso \code{image()}
#' @return a plot of the given matrix
#' @author Mauricio Cantor
#' @export
plotmatbin <- function(x, colpres, colabs, ...){
  par(mar = c(1,1,1,1))

  reverse <- nrow(x) : 1
  x <- x[reverse,]
  
  if(nrow(x)==ncol(x)){ # if the matrix is symmetrical
  image(1:nrow(x), 1:ncol(x), z=t(x) , col=c(colabs, colpres), xlab="", ylab="", axes=FALSE, zlim=c(min(x),max(x)))
  } else {              # if asymmetrical
    image(1:ncol(x), 1:nrow(x), z=t(x) , col=c(colabs, colpres), xlab="", ylab="", axes=FALSE, zlim=c(min(x),max(x)))
  }
}


#' @title Toy data to test Nestedness function
#' @description Creates square, symmetric binary matrices with diagonal=0 to test the function \code{unodf()}
#' @param size Number of row and columns
#' @return A square, symmetric binary matrices with diagonal=0
#' @examples perfnestgen(10)
#' @author Mauricio Cantor
#' @export
perfnestgen = function(size){
  A <- matrix(NA, size, size)
  for(i in 1:size){
    aux1  <- c(rep(1, size))
    if(i!=1) aux1[(size-(i-2)):size] = 0
    A[i,] <- aux1
  }
  diag(A) <- 0
  A
}



#' @title Ordering matrices
#' @description This functions takes a binary or weighted matrix and re-order its rows and columns according to partial totals, i.e. first rows and columns are the ones with higher sums.
#' @param mat Matrix to be ordered according to row and column sums
#' @param symmetry FALSE for asymetric matrices (different number of rows and columns) depicting two-mode networks; TRUE for symmetric matrices depicting one-mode networks
#' @return A matrix ordered by row and column sums
#' @examples 
#'# A weighted asymetric matrix
#'matwa =  rbind(c(0.1,0.9,0,0,0,0.3,1),
#'               c(0,0,0,0.6,0.7,0,0.1),
#'               c(0,0.2,0,0.4,0.4,0.4,0.2),
#'               c(0.6,0,0.7,0,0.5,0,0.3),
#'               c(0,0.2,0.6,0,0,0.8,0))
#'
#'# A binary asymetric matrix
#'matba = rbind(c(1,1,0,0,0,1,1),
#'              c(0,0,0,1,1,0,1),
#'              c(0,1,0,1,1,1,1),
#'              c(1,0,1,0,1,0,1),
#'              c(0,1,1,0,0,1,0))
#'              
#'# A weighted symmetric matrix
#'matws =  rbind(c(0,0.1,0.2,0,0.3),
#'               c(0.1,0,0.8,0.5,0),
#'               c(0.2,0.8,0,0,0.9),
#'               c(0,0.5,0,0,0.4),
#'               c(0.3,0,0.9,0.4,0))
#'
#'#A binary symmetric matrix
#'matbs =  rbind(c(0,1,1,0,1),
#'               c(1,0,1,1,0),
#'               c(1,1,0,0,1),
#'               c(0,1,0,0,1),
#'               c(1,0,1,1,0))
#' 
#' # Reordering matrices
#' ordermat(matwa, symmetry=FALSE)
#' ordermat(matba, symmetry=FALSE)
#' ordermat(matws, symmetry=TRUE)
#' ordermat(matbs, symmetry=TRUE)
#' 
#' @author Mauricio Cantor
#' @export
ordermat <- function(mat, symmetry=FALSE){
  if(symmetry==FALSE){
    mat <- mat[order(rowSums(mat), decreasing=TRUE), order(colSums(mat), decreasing=TRUE)]    
  } else {
    mat <- mat[order((rowSums(mat) + colSums(mat)), decreasing = TRUE), order(rowSums(mat) + colSums(mat), decreasing = TRUE)]  
  }
  mat
}