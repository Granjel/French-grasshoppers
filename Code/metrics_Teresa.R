#TERESA

#function to bring back the alphas to a list
alphas_4reps <- as.data.frame(read.table("Results/metrics_4rep_alphas.txt", header = TRUE, sep = "\t")) ########################## CHANGE

list_from_df <- function(df, nspecies){
  
}

alphas_4reps <- list_from_df(alphas_4reps, 22)


#working
intr_4reps <- as.data.frame(read.table("Results/metrics_4rep_intr.txt", header = TRUE, sep = "\t")) ########################## CHANGE
steps <- seq(from = 1, to = nrow(intr_4reps), by = nspecies)
intrinsic_4reps <- list()
for (i in 1:length(steps)){
  pos <- steps[i]
  intrinsic_4reps[[i]] <- intr_4reps[pos:(pos+(nspecies-1)),1]
  names(intrinsic_4reps[[i]]) <- names(alphas_4reps[[1]])
}


###TOOLBOX_COEXISTENCE ##########################################################################################################

#install.packages("mvtnorm")
require(mvtnorm)

#input parameters:
#alpha = competition strenght matrix 
#r = vector of intrinsic growth rates

#structural niche difference (output on a log scale)
Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}

#vector defining the centroid of the feasibility domain
r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}


#structural fitness difference (in degree)
theta <- function(alpha,r){
  r_c <- r_centroid(alpha)
  out <- acos(sum(r_c*r)/(sqrt(sum(r^2))*sqrt(sum(r_c^2))))*180/pi
  return(out)
}


#test if a system (alpha and r) is feasible (output 1 = feasible, 0 = not feasible)
test_feasibility <- function(alpha,r){
  out <- prod(solve(alpha,r)>0)
  return(out)
}


#test which pairs in a system (alpha and r) are feasible (output 1 = feasible, 0 = not feasible)
test_feasibility_pairs <- function(alpha,r){
  n <- length(r)
  c <- combn(n,2)
  nc <- dim(c)[2]
  f <- rep(NA,nc)
  for (i in 1:nc){
    f[i] <- prod(solve(alpha[c[,i],c[,i]],r[c[,i]])>0)
  }
  out <- list(pairs = c, feasibility = f)
  return(out)
}


#compute the feasiblity domain, the feasibility domain of all pairs, and their overlap (Nrand = number of randomization)
compute_overlap <- function(alpha,Nrand){
  
  n <- dim(alpha)[1]
  
  counter_f <- 0
  counter_overlap <- 0
  counter_all <- 0
  
  for (i in 1:Nrand){
    
    r_rand <- abs(rnorm(n))  
    r_rand <- r_rand/sqrt(sum(r_rand^2))
    
    f1 <- test_feasibility(alpha,r_rand)  
    f2 <- test_feasibility_pairs(alpha,r_rand)$feasibility  
    
    counter_f <- counter_f + f1
    counter_all <- counter_all + prod(f2)
    counter_overlap <- counter_overlap + f1*prod(f2)
    
  }
  
  Omega <- counter_f/Nrand
  Omega_all <- counter_all/Nrand
  overlap <- counter_overlap/Nrand
  
  out <- list(Omega = Omega, Omega_all = Omega_all, overlap = overlap)
  return(out)
  
}

#########################################################################################################################################

### function --- structural coexistence for combinations of **3 species** ###############################################################
structural_coex_3spp <- function(alpha, intrinsic){
  
  n <- 3 #spp
  combn(n, 2)
  pair_names <- c("coex_1_2", "coex_1_3", "coex_2_3")
  
  combos3 <- t(combn(rownames(alpha), n))
  results_combos3 <- matrix(nrow=dim(combos3)[1], ncol=(5 + ncol(combn(n, 2)) + 1))
  row.names(results_combos3) <- apply(combos3, 1, paste, collapse = ".") 
  colnames(results_combos3) <- c("Omega", "theta", "differential", "overlap", "feasibility", pair_names, "coex_rate")
  
  for(i in 1:nrow(combos3)){
    alpha2 <- as.matrix(alpha[combos3[i,], combos3[i,]])
    intrinsic2 <- as.matrix(subset(intrinsic, rownames(intrinsic) %in% combos3[i,]))
    #omega
    results_combos3[i,1] <- 10^Omega(alpha2)
    #theta
    results_combos3[i,2] <- theta(alpha2, intrinsic2)
    #differential and overlap
    y <- compute_overlap(alpha2,10000)
    results_combos3[i,3] <- y$Omega - y$Omega_all
    results_combos3[i,4] <- y$overlap
    #feasibility (all)
    results_combos3[i,5] <- test_feasibility(alpha2, intrinsic2)
    x <- test_feasibility_pairs(alpha2, intrinsic2)
    #feasibility (pair by pair)
    results_combos3[i,6:(5+ncol(combn(n, 2)))] <- x$feasibility
    #coexistence rate
    results_combos3[i,ncol(results_combos3)] <- sum(x$feasibility) / ncol(combn(n, 2))
  }
  results_combos3 <- as.data.frame(results_combos3)
  
  return(results_combos3)
}

#########################################################################################################################################


# OUTPUT WITH PARALLEL COMPUTING
library(foreach)
library(doParallel)
library(parallel)
library(MASS)
source("Code/functions_structural_coex_outputs.R")
#source("Code/toolbox_coexistence.R") #included in the previous source

#flatten the lists
alphas <- alphas_4reps
intrinsic <- intrinsic_4reps

#iniciate the implicit cluster
numCores <- detectCores() #number of cores in your computer
registerDoParallel(numCores) #important: define parallel computation to numCores

#metrics alone... RUN!
start_time <- Sys.time()
#foreach() to distribute the loops between all the cores of the computer: faster
metrics_Teresa <- foreach (i = 5:length(alphas), .combine = rbind, .packages='mvtnorm') %dopar% { #you need to say which packages you are using
  structural_coex_3spp(alpha = as.data.frame(alphas[[i]]), intrinsic = as.matrix(intrinsic[[i]]))
}
end_time <- Sys.time()
end_time - start_time #time expent running code

stopImplicitCluster() #end parallel computation

write.table(metrics_Teresa, file = "Results/structural_coex_metrics_Teresa.txt", sep = "\t", row.names = FALSE)





























