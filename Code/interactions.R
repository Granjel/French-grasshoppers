#script to calculate the coexistence outputs for the interactions between VIM and network metrics

library(bipartite)

#nested(sigma)
#mean(as.numeric(lapply(deletions25, nested)))
#sd(as.numeric(lapply(deletions25, nested)))
#mean(as.numeric(lapply(deletions50, nested)))
#sd(as.numeric(lapply(deletions50, nested)))
#mean(as.numeric(lapply(deletions75, nested)))
#sd(as.numeric(lapply(deletions75, nested)))

a <- nested(sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE)

n <- NULL

for (i in 1:reps){
  n <- c(n, nested(deletions25[[i]], method = "weighted NODF", rescale = FALSE, normalised = FALSE),
         nested(deletions50[[i]], method = "weighted NODF", rescale = FALSE, normalised = FALSE),
         nested(deletions75[[i]], method = "weighted NODF", rescale = FALSE, normalised = FALSE))
}

deletions_percentage <- c(1, rep(25, 20), rep(50, 20), rep(75, 20))
nestedness <- c(a, n)

df <- data.frame(deletions_percentage, nestedness)

plot(df)

mean(df$nestedness[df$deletions_percentage == 25])
mean(df$nestedness[df$deletions_percentage == 50])
mean(df$nestedness[df$deletions_percentage == 75])




networklevel(t(sigma))

wine(sigma, nreps = 10)$'win'
plot(wine(sigma, nreps = 1000))

wine(deletions25[[4]], nreps = 10)$'win'
wine(deletions25[[4]], nreps = 10)$wine


V.ratio(sigma) #Calculates the variance-ratio as suggested by Schluter (1984)
computeModules(sigma)
nested
degreedistr(sigma)


###################################################################################
###################################################################################


#seleccionar matrices variando la conectancia

sigma <- as.matrix(read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1))

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

##
sigma <- as.matrix(read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1))

conectance_matrices <- list()
modifs_conectance <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) #conectance values --- modify as needed
reps <- 10 #reps for each conectance value --- modify as needed
nsigma <- as.numeric(nrow(sigma)*ncol(sigma))
value <- 0 #zero or NA?

for (i in 1:length(modifs_conectance)){
  conectance_matrices[[i]] <- list()
  ndeletions <- 1 - modifs_conectance[i]
  for (j in 1:reps){
    modify_sigma <- sigma
    modify_sigma[sample(nsigma, ndeletions)] <- value
    conectance_matrices[[i]][[j]] <- modify_sigma
  }
}


conectance_matrices[[1]][[1]]







sigma <- as.matrix(read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1))

conectance_values <- c(0.2, 0.4, 0.6, 0.8) #modif if needed
conectance_matrices <- list() #selected matrices
order_matrices <- list() #order in which these matrices showed up and were selected
conectance_residuals <- list() #matrices not selected
order_residuals <- list() ##order in which these matrices showed up and weren't selected

reps <- 10 #reps for each conectance value --- modify as needed
nsigma <- as.numeric(nrow(sigma)*ncol(sigma))
value <- 0 #zero or NA?

for (i in 1:conectance_values){
  
  conectance_matrices[[i]] <- list() 
  order_matrices[[i]] <- list()
  conectance_residuals[[i]] <- list()
  order_residuals[[i]] <- list()
  
  ndeletions <- 1 - modifs_conectance[i]
  
  modify_sigma <- sigma
  modify_sigma[sample(nsigma, ndeletions)] <- value
  conectance_matrices[[i]][[j]] <- modify_sigma  
  
  
  
  
}





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









sigma <- as.matrix(read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1))

nsigma <- as.numeric(nrow(sigma)*ncol(sigma))
value <- 0 #zero or NA?
ndeletions <- nsigma * 0.8

modify_sigma <- sigma
modify_sigma[sample(nsigma, ndeletions)] <- value
 

library(bipartite)

nested(sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE)
nested(modify_sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE)



#check the possible range and distribution of nestedness
sigma <- as.matrix(read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1))

percentage <- seq(from = 0.2, to = 0.8, by = 0.2)
reps <- 1000
nsigma <- as.numeric(nrow(sigma)*ncol(sigma))
value <- 0 #zero or NA?

value <- 0 #zero or NA?

conectance <- NULL
nest <- NULL
modul <- NULL

for (i in 1:length(percentage)){
  for (j in 1:reps){
    
    modify_sigma <- sigma
    modify_sigma[sample(nsigma, (percentage[i] * nsigma))] <- value
    
    conectance <- c(conectance, (1 - percentage[i]))
    nest <- c(nest, nested(modify_sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE))
    modul <- c(modul,  computeModules(modify_sigma, method = "Beckett")@likelihood)
    
  }
}

df <- data.frame(conectance, nest, modul) #play

boxplot(df$nest ~ df$conectance, xlab = "Conectance (realised links / possible links)", ylab = "Nestedness (weighted MODF)", outline = FALSE)
abline(nested(sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE), 0, lty = "dotted")

boxplot(df$modul ~ df$conectance, xlab = "Conectance (realised links / possible links)", ylab = "Modularity (Beckett's algorithm)", outline = FALSE)
abline(computeModules(sigma, method = "Beckett")@likelihood, 0, lty = "dotted")

plot(df$modul, df$nest, xlab = "Modularity (Beckett's algorithm)", ylab = "Nestedness (weighted MODF)", ylim = c(0, 30), xlim = c(-25, 60))
abline(nested(sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE), 0, lty = "dotted")
abline(v = computeModules(sigma, method = "Beckett")@likelihood, lty = "dotted")



is.numeric(df$conectance)

C2 <- df[which(df$conectance == "0.2"),]
row.names(C2) <- NULL
C4 <- df[which(df$conectance == "0.4"),]
row.names(C4) <- NULL
C6 <- df[which(df$conectance == "0.6"),]
row.names(C6) <- NULL
C8 <- df[which(df$conectance == "0.8"),]
row.names(C8) <- NULL

#conectance 0.2
boxplot(C2$nest, outline = FALSE)
n1_1 <- boxplot.stats(C2$nest)$stats[1] + (boxplot.stats(C2$nest)$stats[2] - boxplot.stats(C2$nest)$stats[1]) * (2/3)
n1_2 <- boxplot.stats(C2$nest)$stats[2] + (boxplot.stats(C2$nest)$stats[3] - boxplot.stats(C2$nest)$stats[2]) / 2
n1_3 <- boxplot.stats(C2$nest)$stats[3] + (boxplot.stats(C2$nest)$stats[4] - boxplot.stats(C2$nest)$stats[3]) / 2
n1_4 <- boxplot.stats(C2$nest)$stats[4] + (boxplot.stats(C2$nest)$stats[5] - boxplot.stats(C2$nest)$stats[4]) * (1/3)
n_e1 <- (boxplot.stats(C2$nest)$stats[5] - boxplot.stats(C2$nest)$stats[1]) * 0.02
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
boxplot(C4$nest, outline = FALSE)
n2_1 <- boxplot.stats(C4$nest)$stats[1] + (boxplot.stats(C4$nest)$stats[2] - boxplot.stats(C4$nest)$stats[1]) * (2/3)
n2_2 <- boxplot.stats(C4$nest)$stats[2] + (boxplot.stats(C4$nest)$stats[3] - boxplot.stats(C4$nest)$stats[2]) / 2
n2_3 <- boxplot.stats(C4$nest)$stats[3] + (boxplot.stats(C4$nest)$stats[4] - boxplot.stats(C4$nest)$stats[3]) / 2
n2_4 <- boxplot.stats(C4$nest)$stats[4] + (boxplot.stats(C4$nest)$stats[5] - boxplot.stats(C4$nest)$stats[4]) * (1/3)
n_e2 <- (boxplot.stats(C4$nest)$stats[5] - boxplot.stats(C4$nest)$stats[1]) * 0.02
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
boxplot(C6$nest, outline = FALSE)
n3_1 <- boxplot.stats(C6$nest)$stats[1] + (boxplot.stats(C6$nest)$stats[2] - boxplot.stats(C6$nest)$stats[1]) * (2/3)
n3_2 <- boxplot.stats(C6$nest)$stats[2] + (boxplot.stats(C6$nest)$stats[3] - boxplot.stats(C6$nest)$stats[2]) / 2
n3_3 <- boxplot.stats(C6$nest)$stats[3] + (boxplot.stats(C6$nest)$stats[4] - boxplot.stats(C6$nest)$stats[3]) / 2
n3_4 <- boxplot.stats(C6$nest)$stats[4] + (boxplot.stats(C6$nest)$stats[5] - boxplot.stats(C6$nest)$stats[4]) * (1/3)
n_e3 <- (boxplot.stats(C6$nest)$stats[5] - boxplot.stats(C6$nest)$stats[1]) * 0.02
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
boxplot(C8$nest, outline = FALSE)
n4_1 <- boxplot.stats(C8$nest)$stats[1] + (boxplot.stats(C8$nest)$stats[2] - boxplot.stats(C8$nest)$stats[1]) * (2/3)
n4_2 <- boxplot.stats(C8$nest)$stats[2] + (boxplot.stats(C8$nest)$stats[3] - boxplot.stats(C8$nest)$stats[2]) / 2
n4_3 <- boxplot.stats(C8$nest)$stats[3] + (boxplot.stats(C8$nest)$stats[4] - boxplot.stats(C8$nest)$stats[3]) / 2
n4_4 <- boxplot.stats(C8$nest)$stats[4] + (boxplot.stats(C8$nest)$stats[5] - boxplot.stats(C8$nest)$stats[4]) * (1/3)
n_e4 <- (boxplot.stats(C8$nest)$stats[5] - boxplot.stats(C8$nest)$stats[1]) * 0.02
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


##################################


#conectance 0.2
boxplot(C2$modul, outline = FALSE)
m1_1 <- boxplot.stats(C2$modul)$stats[1] + (boxplot.stats(C2$modul)$stats[2] - boxplot.stats(C2$modul)$stats[1]) * (2/3)
m1_2 <- boxplot.stats(C2$modul)$stats[2] + (boxplot.stats(C2$modul)$stats[3] - boxplot.stats(C2$modul)$stats[2]) / 2
m1_3 <- boxplot.stats(C2$modul)$stats[3] + (boxplot.stats(C2$modul)$stats[4] - boxplot.stats(C2$modul)$stats[3]) / 2
m1_4 <- boxplot.stats(C2$modul)$stats[4] + (boxplot.stats(C2$modul)$stats[5] - boxplot.stats(C2$modul)$stats[4]) * (1/3)
m_e1 <- (boxplot.stats(C2$modul)$stats[5] - boxplot.stats(C2$modul)$stats[1]) * 0.02
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
boxplot(C4$modul, outline = FALSE)
m2_1 <- boxplot.stats(C4$modul)$stats[1] + (boxplot.stats(C4$modul)$stats[2] - boxplot.stats(C4$modul)$stats[1]) * (2/3)
m2_2 <- boxplot.stats(C4$modul)$stats[2] + (boxplot.stats(C4$modul)$stats[3] - boxplot.stats(C4$modul)$stats[2]) / 2
m2_3 <- boxplot.stats(C4$modul)$stats[3] + (boxplot.stats(C4$modul)$stats[4] - boxplot.stats(C4$modul)$stats[3]) / 2
m2_4 <- boxplot.stats(C4$modul)$stats[4] + (boxplot.stats(C4$modul)$stats[5] - boxplot.stats(C4$modul)$stats[4]) * (1/3)
m_e2 <- (boxplot.stats(C4$modul)$stats[5] - boxplot.stats(C4$modul)$stats[1]) * 0.02
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
boxplot(C6$modul, outline = FALSE)
m3_1 <- boxplot.stats(C6$modul)$stats[1] + (boxplot.stats(C6$modul)$stats[2] - boxplot.stats(C6$modul)$stats[1]) * (2/3)
m3_2 <- boxplot.stats(C6$modul)$stats[2] + (boxplot.stats(C6$modul)$stats[3] - boxplot.stats(C6$modul)$stats[2]) / 2
m3_3 <- boxplot.stats(C6$modul)$stats[3] + (boxplot.stats(C6$modul)$stats[4] - boxplot.stats(C6$modul)$stats[3]) / 2
m3_4 <- boxplot.stats(C6$modul)$stats[4] + (boxplot.stats(C6$modul)$stats[5] - boxplot.stats(C6$modul)$stats[4]) * (1/3)
m_e3 <- (boxplot.stats(C6$modul)$stats[5] - boxplot.stats(C6$modul)$stats[1]) * 0.02
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
boxplot(C8$modul, outline = FALSE)
m4_1 <- boxplot.stats(C8$modul)$stats[1] + (boxplot.stats(C8$modul)$stats[2] - boxplot.stats(C8$modul)$stats[1]) * (2/3)
m4_2 <- boxplot.stats(C8$modul)$stats[2] + (boxplot.stats(C8$modul)$stats[3] - boxplot.stats(C8$modul)$stats[2]) / 2
m4_3 <- boxplot.stats(C8$modul)$stats[3] + (boxplot.stats(C8$modul)$stats[4] - boxplot.stats(C8$modul)$stats[3]) / 2
m4_4 <- boxplot.stats(C8$modul)$stats[4] + (boxplot.stats(C8$modul)$stats[5] - boxplot.stats(C8$modul)$stats[4]) * (1/3)
m_e4 <- (boxplot.stats(C8$modul)$stats[5] - boxplot.stats(C8$modul)$stats[1]) * 0.02
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




GHGHGHGHGHGHG
  JHJGJHJHJG
    KHKHKGKFKF
GJMDKJKFSNFSNHS
  fmnñf
















