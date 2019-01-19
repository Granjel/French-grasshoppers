#INTERACTIONS BETWEEN THE DIFFERENT NETWORK METRICS AND THE VARIABILITY OF THE INTERACTION STRENGTH




#example
try <- list(list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
            list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
            list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
            list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())))

for (i in 1:4){
  for (j in 1:4){
    for (k in 1:4){
      for (q in 1:10){
        try[[i]][[j]][[k]][[q]] <- runif(1, 0, 1)
      }
    }
  }
}
try[[1]][[1]][[1]][[5]]












































