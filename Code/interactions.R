#INTERACTIONS BETWEEN THE DIFFERENT NETWORK METRICS AND THE VARIABILITY OF THE INTERACTION STRENGTH


#
matrices <- list(list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
                 list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
                 list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
                 list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())))

conectance_values <- seq(from = 0.2, to = 0.8, by = 0.2)
reps <- 10 #number of replicates
nsigma <- as.numeric(nrow(sigma)*ncol(sigma)) #number of entries in the matrix, to combine with 'percentage'
value <- 0 #value to write when making the deletions

for (i in 1:length(conectance_values)){ #conectance values
  
  for (j in (1 + (4 * (i - 1))):((1 + (4 * (i - 1))) + 3)){ #nestedness values
    
    for (k in (1 + (4 * (i - 1))):((1 + (4 * (i - 1))) + 3)){ #modularity values
      
      for (q in 1:reps){ #replicates for each combination
        
        condition <- FALSE #condition for the while() loop to meet
        
        while(isFALSE(condition)){ #don't stop until the specific matrix is found
          
          #reset matrix to sigma
          modify_sigma <- sigma
          
          #random, variable percentage of deletions
          modify_sigma[sample(nsigma, ((1 - conectance_values[i]) * nsigma))] <- value
          
          #calculate nestedness and modularity for the random matrix
          nest <- nested(modify_sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE) #nestedness
          modu <- computeModules(modify_sigma, method = "Beckett")@likelihood #modularity
          
          n_upper <- n_values[j] + n_errors[i]
          n_lower <- n_values[j] - n_errors[i]
          m_upper <- m_values[k] + m_errors[i]
          m_lower <- m_values[k] - m_errors[i]
          
          #see if it matches the requirements and, if so, save
          if ( (nest >= n_lower && nest <= n_upper) && (modu >= m_lower && modu <= m_upper) ){ #nestedness and modularity +- error
            matrices[[i]][[j]][[k]][[q]] <- modify_sigma #save to break the while() loop
            condition <- TRUE #condition met
            
          } #end if()
        } #end while()
        print("q")
        print(q)
      } #end q
      print("k")
      print(k)
    } #end k
    print("j")
    print(j)
  } #end j
  print("i")
  print(i)
} #end i


nested(matrices[[1]][[1]][[1]][[1]], method = "weighted NODF", rescale = FALSE, normalised = FALSE)
computeModules(matrices[[1]][[1]][[1]][[1]], method = "Beckett")@likelihood

nested(matrices[[1]][[1]][[1]][[2]], method = "weighted NODF", rescale = FALSE, normalised = FALSE)
computeModules(matrices[[1]][[1]][[1]][[2]], method = "Beckett")@likelihood

nested(matrices[[1]][[1]][[1]][[3]], method = "weighted NODF", rescale = FALSE, normalised = FALSE)
computeModules(matrices[[1]][[1]][[1]][[3]], method = "Beckett")@likelihood













