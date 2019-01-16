alpha <- as.matrix(read.table("Results/output_glinternet/alpha.txt", header = TRUE, sep = "\t") * (-1))
sigma <- as.matrix(read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1))
lambda <- as.matrix(read.table("Results/output_glinternet/lambda.txt", header = TRUE, sep = "\t"))
gamma1 <- as.matrix(read.table("Results/output_glinternet/gamma1.txt", header = TRUE, sep = "\t") * (-1))
gamma2 <- as.matrix(read.table("Results/output_glinternet/gamma2.txt", header = TRUE, sep = "\t") * (-1))
gamma3 <- as.matrix(read.table("Results/output_glinternet/gamma3.txt", header = TRUE, sep = "\t") * (-1))
gamma4 <- as.matrix(read.table("Results/output_glinternet/gamma4.txt", header = TRUE, sep = "\t") * (-1))
gamma5 <- as.matrix(read.table("Results/output_glinternet/gamma5.txt", header = TRUE, sep = "\t") * (-1))
gamma6 <- as.matrix(read.table("Results/output_glinternet/gamma6.txt", header = TRUE, sep = "\t") * (-1))

value <- 0 #zero or NA?
nsigma <- as.numeric(nrow(sigma)*ncol(sigma))
nmiss <- nsigma * 0.25
reps <- 20

deletions25 <- list()
deletions50 <- list()
deletions75 <- list()
for (i in 1:reps){
  sigma25 <- sigma
  sigma50 <- sigma
  sigma75 <- sigma
  
  sigma25[sample(nsigma, nmiss)] <- value
  sigma50[sample(nsigma, nmiss*2)] <- value
  sigma75[sample(nsigma, nmiss*3)] <- value
  
  deletions25[[i]] <- sigma25
  deletions50[[i]] <- sigma50
  deletions75[[i]] <- sigma75
}

rm(i, nsigma, nmiss, sigma25, sigma50, sigma75)


#from sigma to gamma
source("Code/sigma_to_gamma.R")

gamma <- list(gamma1, gamma2, gamma3, gamma4, gamma5, gamma6)

matrices25 <- list()
matrices50 <- list()
matrices75 <- list()
for (i in 1:reps){
  matrices25[[i]] <- deletion_gamma(deletions25[[i]], gamma, value)
  matrices50[[i]] <- deletion_gamma(deletions50[[i]], gamma, value)
  matrices75[[i]] <- deletion_gamma(deletions75[[i]], gamma, value)
}
rm(i, gamma)



#intrinsic
for(i in 1:reps){
  deletions25[[i]] <- as.vector(lambda) + apply(deletions25[[i]], 1, sum)
  deletions50[[i]] <- as.vector(lambda) + apply(deletions50[[i]], 1, sum)
  deletions75[[i]] <- as.vector(lambda) + apply(deletions75[[i]], 1, sum)
}

#alpha
matrices25 <- lapply(matrices25, matrix_sum)
matrices50 <- lapply(matrices50, matrix_sum)
matrices75 <- lapply(matrices75, matrix_sum)

for(i in 1:reps){
  matrices25[[i]] <- matrices25[[i]] + alpha
  matrices50[[i]] <- matrices50[[i]] + alpha
  matrices75[[i]] <- matrices75[[i]] + alpha
}


# OUTPUT WITH PARALLEL COMPUTING
library(foreach)
library(doParallel)
library(parallel)
library(MASS)
source("Code/functions_structural_coex_outputs.R")
#source("Code/toolbox_coexistence.R") #included in the previous source


numCores <- detectCores() #number of cores in your computer
registerDoParallel(numCores) #important: define parallel computation to numCores

##25%
start_time25 <- Sys.time()
metrics25 <- foreach (i = 1:reps, .combine = rbind, .packages='mvtnorm') %dopar% {
  structural_coex_3spp(alpha = as.data.frame(matrices25[[i]]), intrinsic = as.matrix(deletions25[[i]]))
}
end_time25 <- Sys.time()
end_time25 - start_time25 #time expent running code

##50%
start_time50 <- Sys.time()
metrics50 <- foreach (i = 1:reps, .combine = rbind, .packages='mvtnorm') %dopar% {
  structural_coex_3spp(alpha = as.data.frame(matrices50[[i]]), intrinsic = as.matrix(deletions50[[i]]))
}
end_time50 <- Sys.time()
end_time50 - start_time50 #time expent running code

##75%
start_time75 <- Sys.time()
metrics75 <- foreach (i = 1:reps, .combine = rbind, .packages='mvtnorm') %dopar% {
  structural_coex_3spp(alpha = as.data.frame(matrices75[[i]]), intrinsic = as.matrix(deletions75[[i]]))
}
end_time75 <- Sys.time()
end_time75 - start_time75 #time expent running code

stopImplicitCluster() #end parallel computation


#replicate number
r <- NULL
for (i in 1:reps){
  r <- c(r, rep(i, ncol(combn(22, 3))))
}
metrics25$replicate <- r
metrics50$replicate <- r
metrics75$replicate <- r

write.table(metrics25, file = "Results/structural_coex_metrics25.txt", sep = "\t", row.names = FALSE)
write.table(metrics50, file = "Results/structural_coex_metrics50.txt", sep = "\t", row.names = FALSE)
write.table(metrics75, file = "Results/structural_coex_metrics75.txt", sep = "\t", row.names = FALSE)



#boxplot(t(sigma), main = "Overall direct effect of grasshoppers on plants",
#        ylab = "Herbivory coefficient (sum of the 6 grasshopper spp.)", las = 2)
#abline(0,0)

#next step -> remove interactions on the gamma matrices following the deletions over sigma 


##check deviation
#media <- NULL
#deviation <- NULL
#for (j in 1:length(deletions)){
#  media <- c(media, mean(deletions[[j]]))
#  deviation <- c(deviation, sd(deletions[[j]]))
#}
#mean(media)
#sd(media)
#mean(deviation)
#sd(deviation)

