#deletion from sigma to gamma

###function 1: first step
deletion_gamma <- function(deletion, gamma, value){
  for (i in 1:nrow(deletion)){
    for (j in 1:ncol(deletion)){
      if (deletion[i,j] == value){
        gamma[[j]][,i] <- value
      }
    }
  }
  return(gamma)
}

#sum of the matrices within a list
matrix_sum <- function(A){ #object A is a list of matrices
  B <- A[[1]]
  for (i in 2:length(A)){
    B <- B + A[[i]]
  }
  return(B) #returns a matrix that is the sum of all the matrices within the list
}

#flattenlist (awesome!)
flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}