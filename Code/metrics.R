#CALCULATE THE NETWORK METRICS AND FIND MATRICES FOR COMPUTATION ACCORDING TO THOSE METRIC VALUES

RANDOMIZATION <- FALSE #WARNING - if TRUE, the script computes again the random matrices and (slightly) changes the results! Keep FALSE!

#load empirical matrix 'sigma'
sigma <- as.matrix(read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory

#load bipartite for network metrics: conectance, nestedness and moudlarity
library(bipartite)

#RANDOMIZATION COMPUTATION
if (isTRUE(RANDOMIZATION)){
  
  #check the possible range and distribution of nestedness
  percentage <- seq(from = 0.2, to = 0.8, by = 0.2) #sequence of percentage values (% of deletions)
  reps <- 9999 #number of repetitions per percentage value
  nsigma <- as.numeric(nrow(sigma)*ncol(sigma)) #number of entries in the matrix, to combine with 'percentage'
  value <- 0 #value to write when making the deletions
  
  conectance <- NULL #to save conectance values
  nest <- NULL #to save nestedness values
  modul <- NULL #to save modularity values
  
  #calculate the random matrices (n = 'reps') for each conectance value (1 - 'percentage')
  for (i in 1:length(percentage)){
    for (j in 1:reps){ #as many cicles as 'reps' for each 'percentage' value
      #reset matrix to sigma
      modify_sigma <- sigma
      #random, variable percentage of deletions
      modify_sigma[sample(nsigma, (percentage[i] * nsigma))] <- value
      #calculate the network metrics and save the values
      conectance <- c(conectance, (1 - percentage[i])) #conectance
      nest <- c(nest, nested(modify_sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE)) #nestedness
      modul <- c(modul,  computeModules(modify_sigma, method = "Beckett")@likelihood) #modularity
    } #end j (percentage value)
  } #end i
  
  random_metrics <- data.frame(conectance, nest, modul) #dataframe to play with
  
  #boxplot nestedness - conectance, to see the nestedness range for each conectance value
  boxplot(random_metrics$nest ~ random_metrics$conectance, xlab = "Conectance (realised links / possible links)", ylab = "Nestedness (weighted MODF)", outline = FALSE)
  abline(nested(sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE), 0, lty = "dotted") #empirical sigma's nestedness value
  
  #boxplot nestedness - conectance, to see the nestedness range for each conectance value
  boxplot(random_metrics$modul ~ random_metrics$conectance, xlab = "Conectance (realised links / possible links)", ylab = "Modularity (Beckett's algorithm)", outline = FALSE)
  abline(computeModules(sigma, method = "Beckett")@likelihood, 0, lty = "dotted") #empirical sigma's modularity value
  
  #boxplot nestedness - modularity, to have an idea of the relationship between them
  plot(random_metrics$modul, random_metrics$nest, xlab = "Modularity (Beckett's algorithm)", ylab = "Nestedness (weighted MODF)", ylim = c(0, 30), xlim = c(-25, 60))
  abline(nested(sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE), 0, lty = "dotted") #empirical sigma's nestedness value
  abline(v = computeModules(sigma, method = "Beckett")@likelihood, lty = "dotted") #empirical sigma's modularity value
  
  
  #subset 'random_metrics' depending on the conectance value
  C2 <- random_metrics[which(random_metrics$conectance == "0.2"),] #conectance 0.2
  row.names(C2) <- NULL #reset row names
  C4 <- random_metrics[which(random_metrics$conectance == "0.4"),] #conectance 0.4
  row.names(C4) <- NULL #reset row names
  C6 <- random_metrics[which(random_metrics$conectance == "0.6"),] #conectance 0.6
  row.names(C6) <- NULL #reset row names
  C8 <- random_metrics[which(random_metrics$conectance == "0.8"),] #conectance 0.8
  row.names(C8) <- NULL #reset row names
  
  percentage_error <- 0.05 #percentage of the variability that is allowed as error
  
  #Calculate the NESTEDNESS values, for each conectance level, to be selected for computation, based on the natural data distribution
  #conectance 0.2
  boxplot(C2$nest, ylab = "Nestedness (weighted MODF)", main = "Nestedness; conectance = 0.2", outline = FALSE)
  n1_1 <- boxplot.stats(C2$nest)$stats[1] + (boxplot.stats(C2$nest)$stats[2] - boxplot.stats(C2$nest)$stats[1]) * (2/3)
  n1_2 <- boxplot.stats(C2$nest)$stats[2] + (boxplot.stats(C2$nest)$stats[3] - boxplot.stats(C2$nest)$stats[2]) / 2
  n1_3 <- boxplot.stats(C2$nest)$stats[3] + (boxplot.stats(C2$nest)$stats[4] - boxplot.stats(C2$nest)$stats[3]) / 2
  n1_4 <- boxplot.stats(C2$nest)$stats[4] + (boxplot.stats(C2$nest)$stats[5] - boxplot.stats(C2$nest)$stats[4]) * (1/3)
  n_e1 <- (boxplot.stats(C2$nest)$stats[5] - boxplot.stats(C2$nest)$stats[1]) * percentage_error
  abline(n1_1, 0)
  abline(n1_1 + n_e1, 0, lty = "dotted")
  abline(n1_1 - n_e1, 0, lty = "dotted")
  abline(n1_2, 0)
  abline(n1_2 + n_e1, 0, lty = "dotted")
  abline(n1_2 - n_e1, 0, lty = "dotted")
  abline(n1_3, 0)
  abline(n1_3 + n_e1, 0, lty = "dotted")
  abline(n1_3 - n_e1, 0, lty = "dotted")
  abline(n1_4, 0)
  abline(n1_4 + n_e1, 0, lty = "dotted")
  abline(n1_4 - n_e1, 0, lty = "dotted")
  
  #conectance 0.4
  boxplot(C4$nest, ylab = "Nestedness (weighted MODF)", main = "Nestedness; conectance = 0.4", outline = FALSE)
  n2_1 <- boxplot.stats(C4$nest)$stats[1] + (boxplot.stats(C4$nest)$stats[2] - boxplot.stats(C4$nest)$stats[1]) * (2/3)
  n2_2 <- boxplot.stats(C4$nest)$stats[2] + (boxplot.stats(C4$nest)$stats[3] - boxplot.stats(C4$nest)$stats[2]) / 2
  n2_3 <- boxplot.stats(C4$nest)$stats[3] + (boxplot.stats(C4$nest)$stats[4] - boxplot.stats(C4$nest)$stats[3]) / 2
  n2_4 <- boxplot.stats(C4$nest)$stats[4] + (boxplot.stats(C4$nest)$stats[5] - boxplot.stats(C4$nest)$stats[4]) * (1/3)
  n_e2 <- (boxplot.stats(C4$nest)$stats[5] - boxplot.stats(C4$nest)$stats[1]) * percentage_error
  abline(n2_1, 0)
  abline(n2_1 + n_e2, 0, lty = "dotted")
  abline(n2_1 - n_e2, 0, lty = "dotted")
  abline(n2_2, 0)
  abline(n2_2 + n_e2, 0, lty = "dotted")
  abline(n2_2 - n_e2, 0, lty = "dotted")
  abline(n2_3, 0)
  abline(n2_3 + n_e2, 0, lty = "dotted")
  abline(n2_3 - n_e2, 0, lty = "dotted")
  abline(n2_4, 0)
  abline(n2_4 + n_e2, 0, lty = "dotted")
  abline(n2_4 - n_e2, 0, lty = "dotted")
  
  #conectance 0.6
  boxplot(C6$nest, ylab = "Nestedness (weighted MODF)", main = "Nestedness; conectance = 0.6", outline = FALSE)
  n3_1 <- boxplot.stats(C6$nest)$stats[1] + (boxplot.stats(C6$nest)$stats[2] - boxplot.stats(C6$nest)$stats[1]) * (2/3)
  n3_2 <- boxplot.stats(C6$nest)$stats[2] + (boxplot.stats(C6$nest)$stats[3] - boxplot.stats(C6$nest)$stats[2]) / 2
  n3_3 <- boxplot.stats(C6$nest)$stats[3] + (boxplot.stats(C6$nest)$stats[4] - boxplot.stats(C6$nest)$stats[3]) / 2
  n3_4 <- boxplot.stats(C6$nest)$stats[4] + (boxplot.stats(C6$nest)$stats[5] - boxplot.stats(C6$nest)$stats[4]) * (1/3)
  n_e3 <- (boxplot.stats(C6$nest)$stats[5] - boxplot.stats(C6$nest)$stats[1]) * percentage_error
  abline(n3_1, 0)
  abline(n3_1 + n_e3, 0, lty = "dotted")
  abline(n3_1 - n_e3, 0, lty = "dotted")
  abline(n3_2, 0)
  abline(n3_2 + n_e3, 0, lty = "dotted")
  abline(n3_2 - n_e3, 0, lty = "dotted")
  abline(n3_3, 0)
  abline(n3_3 + n_e3, 0, lty = "dotted")
  abline(n3_3 - n_e3, 0, lty = "dotted")
  abline(n3_4, 0)
  abline(n3_4 + n_e3, 0, lty = "dotted")
  abline(n3_4 - n_e3, 0, lty = "dotted")
  
  #conectance 0.8
  boxplot(C8$nest, ylab = "Nestedness (weighted MODF)", main = "Nestedness; conectance = 0.8", outline = FALSE)
  n4_1 <- boxplot.stats(C8$nest)$stats[1] + (boxplot.stats(C8$nest)$stats[2] - boxplot.stats(C8$nest)$stats[1]) * (2/3)
  n4_2 <- boxplot.stats(C8$nest)$stats[2] + (boxplot.stats(C8$nest)$stats[3] - boxplot.stats(C8$nest)$stats[2]) / 2
  n4_3 <- boxplot.stats(C8$nest)$stats[3] + (boxplot.stats(C8$nest)$stats[4] - boxplot.stats(C8$nest)$stats[3]) / 2
  n4_4 <- boxplot.stats(C8$nest)$stats[4] + (boxplot.stats(C8$nest)$stats[5] - boxplot.stats(C8$nest)$stats[4]) * (1/3)
  n_e4 <- (boxplot.stats(C8$nest)$stats[5] - boxplot.stats(C8$nest)$stats[1]) * percentage_error
  abline(n4_1, 0)
  abline(n4_1 + n_e4, 0, lty = "dotted")
  abline(n4_1 - n_e4, 0, lty = "dotted")
  abline(n4_2, 0)
  abline(n4_2 + n_e4, 0, lty = "dotted")
  abline(n4_2 - n_e4, 0, lty = "dotted")
  abline(n4_3, 0)
  abline(n4_3 + n_e4, 0, lty = "dotted")
  abline(n4_3 - n_e4, 0, lty = "dotted")
  abline(n4_4, 0)
  abline(n4_4 + n_e4, 0, lty = "dotted")
  abline(n4_4 - n_e4, 0, lty = "dotted")
  
  #save the selected nestedness values in vectors
  n_values <- c(n1_1, n1_2, n1_3, n1_4, n2_1, n2_2, n2_3, n2_4, n3_1, n3_2, n3_3, n3_4, n4_1, n4_2, n4_3, n4_4)
  n_errors <- c(n_e1, n_e2, n_e3, n_e4)
  
  #Calculate the MODULARITY values, for each conectance level, to be selected for computation, based on the natural data distribution
  #conectance 0.2
  boxplot(C2$modul, ylab = "Modularity (Beckett's algorithm)", main = "Modularity; conectance = 0.2", outline = FALSE)
  m1_1 <- boxplot.stats(C2$modul)$stats[1] + (boxplot.stats(C2$modul)$stats[2] - boxplot.stats(C2$modul)$stats[1]) * (2/3)
  m1_2 <- boxplot.stats(C2$modul)$stats[2] + (boxplot.stats(C2$modul)$stats[3] - boxplot.stats(C2$modul)$stats[2]) / 2
  m1_3 <- boxplot.stats(C2$modul)$stats[3] + (boxplot.stats(C2$modul)$stats[4] - boxplot.stats(C2$modul)$stats[3]) / 2
  m1_4 <- boxplot.stats(C2$modul)$stats[4] + (boxplot.stats(C2$modul)$stats[5] - boxplot.stats(C2$modul)$stats[4]) * (1/3)
  m_e1 <- (boxplot.stats(C2$modul)$stats[5] - boxplot.stats(C2$modul)$stats[1]) * percentage_error
  abline(m1_1, 0)
  abline(m1_1 + m_e1, 0, lty = "dotted")
  abline(m1_1 - m_e1, 0, lty = "dotted")
  abline(m1_2, 0)
  abline(m1_2 + m_e1, 0, lty = "dotted")
  abline(m1_2 - m_e1, 0, lty = "dotted")
  abline(m1_3, 0)
  abline(m1_3 + m_e1, 0, lty = "dotted")
  abline(m1_3 - m_e1, 0, lty = "dotted")
  abline(m1_4, 0)
  abline(m1_4 + m_e1, 0, lty = "dotted")
  abline(m1_4 - m_e1, 0, lty = "dotted")
  
  #conectance 0.4
  boxplot(C4$modul, ylab = "Modularity (Beckett's algorithm)", main = "Modularity; conectance = 0.4", outline = FALSE)
  m2_1 <- boxplot.stats(C4$modul)$stats[1] + (boxplot.stats(C4$modul)$stats[2] - boxplot.stats(C4$modul)$stats[1]) * (2/3)
  m2_2 <- boxplot.stats(C4$modul)$stats[2] + (boxplot.stats(C4$modul)$stats[3] - boxplot.stats(C4$modul)$stats[2]) / 2
  m2_3 <- boxplot.stats(C4$modul)$stats[3] + (boxplot.stats(C4$modul)$stats[4] - boxplot.stats(C4$modul)$stats[3]) / 2
  m2_4 <- boxplot.stats(C4$modul)$stats[4] + (boxplot.stats(C4$modul)$stats[5] - boxplot.stats(C4$modul)$stats[4]) * (1/3)
  m_e2 <- (boxplot.stats(C4$modul)$stats[5] - boxplot.stats(C4$modul)$stats[1]) * percentage_error
  abline(m2_1, 0)
  abline(m2_1 + m_e2, 0, lty = "dotted")
  abline(m2_1 - m_e2, 0, lty = "dotted")
  abline(m2_2, 0)
  abline(m2_2 + m_e2, 0, lty = "dotted")
  abline(m2_2 - m_e2, 0, lty = "dotted")
  abline(m2_3, 0)
  abline(m2_3 + m_e2, 0, lty = "dotted")
  abline(m2_3 - m_e2, 0, lty = "dotted")
  abline(m2_4, 0)
  abline(m2_4 + m_e2, 0, lty = "dotted")
  abline(m2_4 - m_e2, 0, lty = "dotted")
  
  #conectance 0.6
  boxplot(C6$modul, ylab = "Modularity (Beckett's algorithm)", main = "Modularity; conectance = 0.6", outline = FALSE)
  m3_1 <- boxplot.stats(C6$modul)$stats[1] + (boxplot.stats(C6$modul)$stats[2] - boxplot.stats(C6$modul)$stats[1]) * (2/3)
  m3_2 <- boxplot.stats(C6$modul)$stats[2] + (boxplot.stats(C6$modul)$stats[3] - boxplot.stats(C6$modul)$stats[2]) / 2
  m3_3 <- boxplot.stats(C6$modul)$stats[3] + (boxplot.stats(C6$modul)$stats[4] - boxplot.stats(C6$modul)$stats[3]) / 2
  m3_4 <- boxplot.stats(C6$modul)$stats[4] + (boxplot.stats(C6$modul)$stats[5] - boxplot.stats(C6$modul)$stats[4]) * (1/3)
  m_e3 <- (boxplot.stats(C6$modul)$stats[5] - boxplot.stats(C6$modul)$stats[1]) * percentage_error
  abline(m3_1, 0)
  abline(m3_1 + m_e3, 0, lty = "dotted")
  abline(m3_1 - m_e3, 0, lty = "dotted")
  abline(m3_2, 0)
  abline(m3_2 + m_e3, 0, lty = "dotted")
  abline(m3_2 - m_e3, 0, lty = "dotted")
  abline(m3_3, 0)
  abline(m3_3 + m_e3, 0, lty = "dotted")
  abline(m3_3 - m_e3, 0, lty = "dotted")
  abline(m3_4, 0)
  abline(m3_4 + m_e3, 0, lty = "dotted")
  abline(m3_4 - m_e3, 0, lty = "dotted")
  
  #conectance 0.8
  boxplot(C8$modul, ylab = "Modularity (Beckett's algorithm)", main = "Modularity; conectance = 0.8", outline = FALSE)
  m4_1 <- boxplot.stats(C8$modul)$stats[1] + (boxplot.stats(C8$modul)$stats[2] - boxplot.stats(C8$modul)$stats[1]) * (2/3)
  m4_2 <- boxplot.stats(C8$modul)$stats[2] + (boxplot.stats(C8$modul)$stats[3] - boxplot.stats(C8$modul)$stats[2]) / 2
  m4_3 <- boxplot.stats(C8$modul)$stats[3] + (boxplot.stats(C8$modul)$stats[4] - boxplot.stats(C8$modul)$stats[3]) / 2
  m4_4 <- boxplot.stats(C8$modul)$stats[4] + (boxplot.stats(C8$modul)$stats[5] - boxplot.stats(C8$modul)$stats[4]) * (1/3)
  m_e4 <- (boxplot.stats(C8$modul)$stats[5] - boxplot.stats(C8$modul)$stats[1]) * percentage_error
  abline(m4_1, 0)
  abline(m4_1 + m_e4, 0, lty = "dotted")
  abline(m4_1 - m_e4, 0, lty = "dotted")
  abline(m4_2, 0)
  abline(m4_2 + m_e4, 0, lty = "dotted")
  abline(m4_2 - m_e4, 0, lty = "dotted")
  abline(m4_3, 0)
  abline(m4_3 + m_e4, 0, lty = "dotted")
  abline(m4_3 - m_e4, 0, lty = "dotted")
  abline(m4_4, 0)
  abline(m4_4 + m_e4, 0, lty = "dotted")
  abline(m4_4 - m_e4, 0, lty = "dotted")
  
  #save the selected modularity values in vectors
  m_values <- c(m1_1, m1_2, m1_3, m1_4, m2_1, m2_2, m2_3, m2_4, m3_1, m3_2, m3_3, m3_4, m4_1, m4_2, m4_3, m4_4)
  m_errors <- c(m_e1, m_e2, m_e3, m_e4)
  
  #clean the enviroment
  rm(i, j, percentage, conectance, modul, nest, nsigma, reps, value, modify_sigma, percentage_error, C2, C4, C6, C8,
     n1_1, n1_2, n1_3, n1_4, n2_1, n2_2, n2_3, n2_4, n3_1, n3_2, n3_3, n3_4, n4_1, n4_2, n4_3, n4_4, n_e1, n_e2, n_e3, n_e4,
     m1_1, m1_2, m1_3, m1_4, m2_1, m2_2, m2_3, m2_4, m3_1, m3_2, m3_3, m3_4, m4_1, m4_2, m4_3, m4_4, m_e1, m_e2, m_e3, m_e4)
  
  #write results
  write.table(random_metrics, file = "Results/random_metrics.txt", sep = "\t", row.names = FALSE) #random metrics
  write.table(n_values, file = "Results/n_values.txt", sep = "\t", row.names = FALSE) #nestedness values
  write.table(n_errors, file = "Results/n_errors.txt", sep = "\t", row.names = FALSE) #nestedness allowed errors
  write.table(m_values, file = "Results/m_values.txt", sep = "\t", row.names = FALSE) #modularity values
  write.table(m_errors, file = "Results/m_errors.txt", sep = "\t", row.names = FALSE) #modularity allowed errors
  
} else {
  
  random_metrics <- read.table("Results/random_metrics.txt", header = TRUE, sep = "\t") #random metrics
  n_values <- as.vector(as.matrix(read.table("Results/n_values.txt", header = TRUE, sep = "\t"))) #load: nestedness values
  n_errors <- as.vector(as.matrix(read.table("Results/n_errors.txt", header = TRUE, sep = "\t"))) #load: nestedness allowed errors
  m_values <- as.vector(as.matrix(read.table("Results/m_values.txt", header = TRUE, sep = "\t"))) #load: modularity values
  m_errors <- as.vector(as.matrix(read.table("Results/m_errors.txt", header = TRUE, sep = "\t"))) #load: modularity allowed errors
  
} #end if(isTRUE(RANDOMIZATION)

rm(RANDOMIZATION) #not needed anymore


#LET'S FIND THE MATRICES FOR COMPUTATION!
#crazy list to save the matrices
matrices <- list(list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
                 list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
                 list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
                 list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())))

#settings:
conectance_values <- seq(from = 0.2, to = 0.8, by = 0.2) #conectance values
reps <- 4 #number of replicates
nsigma <- as.numeric(nrow(sigma)*ncol(sigma)) #number of entries in the matrix, to combine with 'percentage'
value <- 0 #value to write when making the deletions

#super loop of loops to find and save the matrices!
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
          n_upper <- n_values[j] + n_errors[i] #upper local bound for nestedness
          n_lower <- n_values[j] - n_errors[i] #lower local bound for nestedness
          m_upper <- m_values[k] + m_errors[i] #upper local bound for modularity
          m_lower <- m_values[k] - m_errors[i] #lower local bound for modularity
          #see if it matches the requirements and, if so, save
          if ( (nest >= n_lower && nest <= n_upper) && (modu >= m_lower && modu <= m_upper) ){ #nestedness and modularity +- error
            matrices[[i]][[(j - (4 * (i - 1)))]][[(k - (4 * (i - 1)))]][[q]] <- modify_sigma #save to break the while() loop
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

#clean the environment a bit
rm(i, j, k, q, modify_sigma, nest, modu, n_upper, n_lower, m_upper, m_lower, conectance_values,
   n_values, n_errors, m_values, m_errors, nsigma, random_metrics, condition)


#load empirical matrices
alpha <- as.matrix(read.table("Results/output_glinternet/alpha.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: competition
lambda <- as.matrix(read.table("Results/output_glinternet/lambda.txt", header = TRUE, sep = "\t"))
gamma1 <- as.matrix(read.table("Results/output_glinternet/gamma1.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma2 <- as.matrix(read.table("Results/output_glinternet/gamma2.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma3 <- as.matrix(read.table("Results/output_glinternet/gamma3.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma4 <- as.matrix(read.table("Results/output_glinternet/gamma4.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma5 <- as.matrix(read.table("Results/output_glinternet/gamma5.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma6 <- as.matrix(read.table("Results/output_glinternet/gamma6.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory



#gammas modified (hoes: high order effects)
hoes <- list(list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
             list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
             list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())),
             list(list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list()), list(list(), list(), list(), list())))

#from sigma to gamma
source("Code/sigma_to_gamma.R") #two functions to help here

#List of gamma matrices (from 1 to 6), to use below
gamma <- list(gamma1, gamma2, gamma3, gamma4, gamma5, gamma6)

#loop to asign the correspondent deletions to the gamma matrices
for (i in 1:4){ #conectance
  for (j in 1:4){ #nestedness
    for (k in 1:4){ #modularity
      for (q in 1:reps){ #replicates
        hoes[[i]][[j]][[k]][[q]] <- deletion_gamma(matrices[[i]][[j]][[k]][[q]], gamma, value) #deletions from sigma to gamma(s)
      } #end q
    } #end k
  } #end j
} #end i
rm(i, j, k, q, gamma, deletion_gamma)



#Population dynamics ---



#Intrinsic
intrinsic <- matrices #to save the intrinsic growth rate

for (i in 1:4){ #conectance
  for (j in 1:4){ #nestedness
    for (k in 1:4){ #modularity
      for (q in 1:reps){ #replicates
        intrinsic[[i]][[j]][[k]][[q]] <- as.vector(lambda) + apply(matrices[[i]][[j]][[k]][[q]], 1, sum)
      } #end q
    } #end k
  } #end j
} #end i
rm(i, j, k, q)

##intrinsic
#for(i in 1:reps){
#  deletions25[[i]] <- as.vector(lambda) + apply(deletions25[[i]], 1, sum)
#  deletions50[[i]] <- as.vector(lambda) + apply(deletions50[[i]], 1, sum)
#  deletions75[[i]] <- as.vector(lambda) + apply(deletions75[[i]], 1, sum)
#}

#alphas
alphas <- hoes #new big list to save the sums

for (i in 1:4){ #conectance
  for (j in 1:4){ #nestedness
    for (k in 1:4){ #modularity
      for (q in 1:reps){ #replicates
        #sum the matrices + alpha
        A <- hoes[[i]][[j]][[k]][[q]]
        B <- A[[1]]
        for (w in 2:6){
          B <- B + A[[w]]
        }
        alphas[[i]][[j]][[k]][[q]] <- B + alpha #sum the empirical alpha + the sumatory of modified gammas
      } #end q
    } #end k
  } #end j
} #end i
rm(i, j, k, q, w, A, B, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, lambda, alpha, sigma)



# OUTPUT WITH PARALLEL COMPUTING
library(foreach)
library(doParallel)
library(parallel)
library(MASS)
source("Code/functions_structural_coex_outputs.R")
#source("Code/toolbox_coexistence.R") #included in the previous source

#flatten the lists
alphas <- flattenlist(alphas)
intrinsic <- flattenlist(intrinsic)

#iniciate the implicit cluster
numCores <- detectCores() #number of cores in your computer
registerDoParallel(numCores) #important: define parallel computation to numCores


#metrics alone... RUN!
start_time <- Sys.time()
#foreach() to distribute the loops between all the cores of the computer: faster
metrics <- foreach (i = 1:reps, .combine = rbind, .packages='mvtnorm') %dopar% { #you need to say which packages you are using
  structural_coex_3spp(alpha = as.data.frame(alphas[[i]]), intrinsic = as.matrix(intrinsic[[i]]))
}
end_time <- Sys.time()
end_time <- start_time #time expent running code

stopImplicitCluster() #end parallel computation

write.table(metrics, file = "Results/structural_coex_metrics.txt", sep = "\t", row.names = FALSE)






#NEXT: WRITE THE METRICS WITHIN THE DATAFRAME FOR A BETTER HANDLING

##replicate number
#r <- NULL
#for (i in 1:reps){
#  r <- c(r, rep(i, ncol(combn(22, 3))))
#}
#metrics25$replicate <- r
#metrics50$replicate <- r
#metrics75$replicate <- r
#
#write.table(metrics, file = "Results/structural_coex_metrics.txt", sep = "\t", row.names = FALSE)



