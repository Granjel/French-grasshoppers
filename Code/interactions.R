#INTERACTIONS BETWEEN THE DIFFERENT NETWORK METRICS AND THE VARIABILITY OF THE INTERACTION STRENGTH

#from sigma to gamma
source("Code/sigma_to_gamma.R") #two functions to help here

matrices_4reps <- as.data.frame(read.table("Results/metrics_4rep_matrices.txt", header = TRUE, sep = "\t"))
matrices <- recover_sigma(matrices_4reps, 22, 6)
rm(matrices_4reps)

vis <- c(0.33, 3, 9, 27)

matrices1 <- list()
matrices2 <- list()
matrices3 <- list()
matrices4 <- list()
for(i in 1:length(matrices)){
  matrices1[[i]] <- matrices[[i]] * vis[1]
  matrices2[[i]] <- matrices[[i]] * vis[2]
  matrices3[[i]] <- matrices[[i]] * vis[3]
  matrices4[[i]] <- matrices[[i]] * vis[4]
}

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
hoes1 <- list()
hoes2 <- list()
hoes3 <- list()
hoes4 <- list()

#List of gamma matrices (from 1 to 6), to use below
gamma <- list(gamma1, gamma2, gamma3, gamma4, gamma5, gamma6)
value <- 0

#loop to asign the correspondent deletions to the gamma matrices
for (i in 1:length(matrices)){
  hoes1[[i]] <- deletion_gamma(matrices1[[i]], gamma, value)
  hoes2[[i]] <- deletion_gamma(matrices2[[i]], gamma, value)
  hoes3[[i]] <- deletion_gamma(matrices3[[i]], gamma, value)
  hoes4[[i]] <- deletion_gamma(matrices4[[i]], gamma, value)
}

rm(i, gamma, deletion_gamma)



#Population dynamics ---

#Intrinsic
intrinsic1 <- list() #to save the intrinsic growth rate
intrinsic2 <- list() #to save the intrinsic growth rate
intrinsic3 <- list() #to save the intrinsic growth rate
intrinsic4 <- list() #to save the intrinsic growth rate

for (i in 1:length(matrices)){
  intrinsic1[[i]] <- as.vector(lambda) + apply(matrices1[[i]], 1, sum)
  intrinsic2[[i]] <- as.vector(lambda) + apply(matrices2[[i]], 1, sum)
  intrinsic3[[i]] <- as.vector(lambda) + apply(matrices3[[i]], 1, sum)
  intrinsic4[[i]] <- as.vector(lambda) + apply(matrices4[[i]], 1, sum)
}

rm(i)


#alphas
alphas1 <- list() #new big list to save the sums
alphas2 <- list() #new big list to save the sums
alphas3 <- list() #new big list to save the sums
alphas4 <- list() #new big list to save the sums

for (i in 1:length(matrices)){
  A1 <- hoes1[[i]]
  B1 <- A1[[1]]
  for (j in 2:6){
    B1 <- B1 + A1[[j]]
  }
  alphas1[[i]] <- B1 + alpha #sum the empirical alpha + the sumatory of modified gammas
  
  A2 <- hoes2[[i]]
  B2 <- A2[[1]]
  for (j in 2:6){
    B2 <- B2 + A2[[j]]
  }
  alphas2[[i]] <- B2 + alpha #sum the empirical alpha + the sumatory of modified gammas
  
  A3 <- hoes3[[i]]
  B3 <- A3[[1]]
  for (j in 2:6){
    B3 <- B3 + A3[[j]]
  }
  alphas3[[i]] <- B3 + alpha #sum the empirical alpha + the sumatory of modified gammas
  
  A4 <- hoes4[[i]]
  B4 <- A4[[1]]
  for (j in 2:6){
    B4 <- B4 + A4[[j]]
  }
  alphas4[[i]] <- B4 + alpha #sum the empirical alpha + the sumatory of modified gammas
}

rm(i, j, A1, A2, A3, A4, B1, B2, B3, B4, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, lambda, alpha, matrices, value, vis)


#write!




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
metrics2 <- foreach (i = 1:length(metrics1), .combine = rbind, .packages='mvtnorm') %dopar% { #you need to say which packages you are using
  structural_coex_3spp(alpha = as.data.frame(alphas[[i]]), intrinsic = as.matrix(intrinsic[[i]]))
}
end_time <- Sys.time()
end_time <- start_time #time expent running code

stopImplicitCluster() #end parallel computation

write.table(metrics, file = "Results/structural_coex_metrics.txt", sep = "\t", row.names = FALSE)




#SAVE A LIST AS A DATA.FRAME

###1###

#alphas1
alphas1_df <- data.frame()
for (i in 1:length(alphas1)){
  alphas1_df <- rbind(alphas1_df, as.data.frame(alphas1[[i]]))
}
write.table(alphas1_df, file = "Results/interactions_alphas1.txt", sep = "\t", row.names = FALSE)

#intrinsic1
intrinsic1_df <- as.data.frame(unlist(intrinsic1))
write.table(intrinsic1_df, file = "Results/interactions_intrinsic1.txt", sep = "\t", row.names = FALSE)

#original deletion matrices 1
matrices1_df <- data.frame()
for (i in 1:length(matrices1)){
  matrices1_df <- rbind(matrices1_df, as.data.frame(matrices1[[i]]))
}
write.table(matrices1_df, file = "Results/interactions_matrices1.txt", sep = "\t", row.names = FALSE)


###2###

#alphas2
alphas2_df <- data.frame()
for (i in 1:length(alphas2)){
  alphas2_df <- rbind(alphas2_df, as.data.frame(alphas2[[i]]))
}
write.table(alphas2_df, file = "Results/interactions_alphas2.txt", sep = "\t", row.names = FALSE)

#intrinsic2
intrinsic2_df <- as.data.frame(unlist(intrinsic2))
write.table(intrinsic2_df, file = "Results/interactions_intrinsic2.txt", sep = "\t", row.names = FALSE)

#original deletion matrices 2
matrices2_df <- data.frame()
for (i in 1:length(matrices2)){
  matrices2_df <- rbind(matrices2_df, as.data.frame(matrices2[[i]]))
}
write.table(matrices2_df, file = "Results/interactions_matrices2.txt", sep = "\t", row.names = FALSE)


###3###

#alphas3
alphas3_df <- data.frame()
for (i in 1:length(alphas3)){
  alphas3_df <- rbind(alphas3_df, as.data.frame(alphas3[[i]]))
}
write.table(alphas3_df, file = "Results/interactions_alphas3.txt", sep = "\t", row.names = FALSE)

#intrinsic3
intrinsic3_df <- as.data.frame(unlist(intrinsic3))
write.table(intrinsic3_df, file = "Results/interactions_intrinsic3.txt", sep = "\t", row.names = FALSE)

#original deletion matrices 3
matrices3_df <- data.frame()
for (i in 1:length(matrices3)){
  matrices3_df <- rbind(matrices3_df, as.data.frame(matrices3[[i]]))
}
write.table(matrices3_df, file = "Results/interactions_matrices3.txt", sep = "\t", row.names = FALSE)


###4###

#alphas4
alphas4_df <- data.frame()
for (i in 1:length(alphas4)){
  alphas4_df <- rbind(alphas4_df, as.data.frame(alphas4[[i]]))
}
write.table(alphas4_df, file = "Results/interactions_alphas4.txt", sep = "\t", row.names = FALSE)

#intrinsic4
intrinsic4_df <- as.data.frame(unlist(intrinsic4))
write.table(intrinsic4_df, file = "Results/interactions_intrinsic4.txt", sep = "\t", row.names = FALSE)

#original deletion matrices 4
matrices4_df <- data.frame()
for (i in 1:length(matrices4)){
  matrices4_df <- rbind(matrices4_df, as.data.frame(matrices4[[i]]))
}
write.table(matrices4_df, file = "Results/interactions_matrices4.txt", sep = "\t", row.names = FALSE)











