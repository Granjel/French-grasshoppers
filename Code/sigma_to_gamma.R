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

#conectance (specific for this study); to broaden its horizons, remove the part between "###"
conect <- function(M){
  r <- nrow(as.data.frame(M))
  c <- ncol(as.data.frame(M))
  total <- r*c
  actual <- total
  for (i in 1:r){
    for (j in 1:c){
      if (M[i,j] == 0){
        actual <- c(actual - 1)
      }
    }
  }
  con <- actual/total
  ###
  if(con >= 0.19 && con <= 0.21){
    con <- 0.2
  } else {
    if (con >= 0.39 && con <= 0.41){
      con <- 0.4
    } else {
      if (con >= 0.59 && con <= 0.61){
        con <- 0.6
      } else {
        if (con >= 0.79 && con <= 0.81){
          con <- 0.8
        }
      }
    }
  }
  ###
  return(con)
}