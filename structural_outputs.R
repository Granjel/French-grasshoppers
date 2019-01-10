source("Code/toolbox_coexistence.R")

#only requires alpha (called 'alpha') and lambda (called 'intrinsic')
alpha <- read.table("Results/output_glinternet/alpha.txt", header = TRUE, sep = "\t") * (-1)
lambda <- read.table("Results/output_glinternet/lambda.txt", header = TRUE, sep = "\t")
intrinsic <- lambda
rm(lambda)
rownames(intrinsic) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                         "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

#################
### 3 species ########################################################################################################################
#################
start_time <- Sys.time()
n <- 3 #spp
combn(n, 2) #help for the next line
pair_names <- c("coex_1_2", "coex_1_3", "coex_2_3") #modify with the number of species

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

end_time <- Sys.time()
end_time - start_time #Time expent running code

rm(alpha2, combos3, intrinsic2, start_time, end_time, i, n, pair_names, x, y)


#################
### 4 species ########################################################################################################################
#################
start_time <- Sys.time()
n <- 4 #spp
combn(n, 2) #help for the next line
pair_names <- c("coex_1_2", "coex_1_3", "coex_1_4", "coex_2_3", "coex_2_4", "coex_3_4") #modify with the number of species

combos4 <- t(combn(rownames(alpha), n))
results_combos4 <- matrix(nrow=dim(combos4)[1], ncol=(5 + ncol(combn(n, 2)) + 1))
row.names(results_combos4) <- apply(combos4, 1, paste, collapse = ".") 
colnames(results_combos4) <- c("Omega", "theta", "differential", "overlap", "feasibility", pair_names, "coex_rate")

for(i in 1:nrow(combos4)){
  alpha2 <- as.matrix(alpha[combos4[i,], combos4[i,]])
  intrinsic2 <- as.matrix(subset(intrinsic, rownames(intrinsic) %in% combos4[i,]))
  #omega
  results_combos4[i,1] <- 10^Omega(alpha2)
  #theta
  results_combos4[i,2] <- theta(alpha2, intrinsic2)
  #differential and overlap
  y <- compute_overlap(alpha2,10000)
  results_combos4[i,3] <- y$Omega - y$Omega_all
  results_combos4[i,4] <- y$overlap
  #feasibility (all)
  results_combos4[i,5] <- test_feasibility(alpha2, intrinsic2)
  x <- test_feasibility_pairs(alpha2, intrinsic2)
  #feasibility (pair by pair)
  results_combos4[i,6:(5+ncol(combn(n, 2)))] <- x$feasibility
  #coexistence rate
  results_combos4[i,ncol(results_combos4)] <- sum(x$feasibility) / ncol(combn(n, 2))
}
results_combos4 <- as.data.frame(results_combos4)

end_time <- Sys.time()
end_time - start_time #Time expent running code

rm(alpha2, combos4, intrinsic2, start_time, end_time, i, n, pair_names, x, y)









#function --- structural coexistence for combinations of **3 species**
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

#function --- structural coexistence for combinations of **4 species**
structural_coex_4spp <- function(alpha, intrinsic){
  
  n <- 4 #spp
  combn(n, 2) #help for the next line
  pair_names <- c("coex_1_2", "coex_1_3", "coex_1_4", "coex_2_3", "coex_2_4", "coex_3_4") #modify with the number of species
  
  combos4 <- t(combn(rownames(alpha), n))
  results_combos4 <- matrix(nrow=dim(combos4)[1], ncol=(5 + ncol(combn(n, 2)) + 1))
  row.names(results_combos4) <- apply(combos4, 1, paste, collapse = ".") 
  colnames(results_combos4) <- c("Omega", "theta", "differential", "overlap", "feasibility", pair_names, "coex_rate")
  
  for(i in 1:nrow(combos4)){
    alpha2 <- as.matrix(alpha[combos4[i,], combos4[i,]])
    intrinsic2 <- as.matrix(subset(intrinsic, rownames(intrinsic) %in% combos4[i,]))
    #omega
    results_combos4[i,1] <- 10^Omega(alpha2)
    #theta
    results_combos4[i,2] <- theta(alpha2, intrinsic2)
    #differential and overlap
    y <- compute_overlap(alpha2,10000)
    results_combos4[i,3] <- y$Omega - y$Omega_all
    results_combos4[i,4] <- y$overlap
    #feasibility (all)
    results_combos4[i,5] <- test_feasibility(alpha2, intrinsic2)
    x <- test_feasibility_pairs(alpha2, intrinsic2)
    #feasibility (pair by pair)
    results_combos4[i,6:(5+ncol(combn(n, 2)))] <- x$feasibility
    #coexistence rate
    results_combos4[i,ncol(results_combos4)] <- sum(x$feasibility) / ncol(combn(n, 2))
  }
  results_combos4 <- as.data.frame(results_combos4)
  
  return(results_combos4)
}



a <- alpha[1:3, 1:3]
b <- intrinsic[1:3]

structural_coex_3spp(a, b)




























