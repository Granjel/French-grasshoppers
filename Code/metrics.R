sigma <- as.matrix(read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1))
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

rm(i, value, nsigma, nmiss, reps, sigma25, sigma50, sigma75)

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

