alpha <- read.table("Results/output_glinternet/alpha.txt", header = TRUE, sep = "\t") * (-1)
sigma <- read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1)
lambda <- read.table("Results/output_glinternet/lambda.txt", header = TRUE, sep = "\t")
gamma1 <- read.table("Results/output_glinternet/gamma1.txt", header = TRUE, sep = "\t") * (-1)
gamma2 <- read.table("Results/output_glinternet/gamma2.txt", header = TRUE, sep = "\t") * (-1)
gamma3 <- read.table("Results/output_glinternet/gamma3.txt", header = TRUE, sep = "\t") * (-1)
gamma4 <- read.table("Results/output_glinternet/gamma4.txt", header = TRUE, sep = "\t") * (-1)
gamma5 <- read.table("Results/output_glinternet/gamma5.txt", header = TRUE, sep = "\t") * (-1)
gamma6 <- read.table("Results/output_glinternet/gamma6.txt", header = TRUE, sep = "\t") * (-1)


#STRENGTH MODIFICATIONS
strength_seq <- seq(from = 0.01, to = 1, by = 0.01) #change as needed

### sigma -> lambda
intr_strength <- list()
xx <- matrix(nrow = 22, ncol = 1)
colnames(xx) <- c("intrinsic")
rownames(xx) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                  "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")
for (i in 1:length(strength_seq)){
  for (j in 1:nrow(xx)){
    xx[j, 1] <- lambda[j,] + rowSums(sigma[j,]) * strength_seq[i]
  }
  intr_strength[[i]] <- xx
}

### gamma(s) -> alpha
gamma <- gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6
comp_strength <- list()
for (i in 1:length(strength_seq)){
  yy <- alpha + gamma * strength_seq[i]
  comp_strength[[i]] <- yy
}

#rm
rm(i, j, xx, yy)


#OUTPUT
### --- strength modification
source("Code/toolbox_coexistence.R")
start_time <- Sys.time()
alpha <- comp_strength[[1]]
intrinsic <- intr_strength[[1]]
output_strength_3spp <- structural_coex_3spp(alpha, intrinsic)
for (i in 2:length(comp_strength)){
  alpha <- comp_strength[[i]]
  intrinsic <- intr_strength[[i]]
  output_strength_3spp <- rbind(output_strength_3spp, structural_coex_3spp(alpha, intrinsic))
}
end_time <- Sys.time()
end_time - start_time #Time expent running code









#VARIABILITY MODIFICATIONS (sd)
variability_seq <- seq(from = 1, to = 100, by = 1) #change as needed

### sigma -> lambda
intr_variability <- list()
xx <- matrix(nrow = 22, ncol = 1)
colnames(xx) <- c("intrinsic")
rownames(xx) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                  "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")
for (i in 1:length(variability_seq)){
  for (j in 1:nrow(xx)){
    xx[j, 1] <- lambda[j,] + rowSums(sigma[j,]) * variability_seq[i]
  }
  intr_variability[[i]] <- xx
}

### gamma(s) -> alpha
gamma <- gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6
comp_variability <- list()
for (i in 1:length(variability_seq)){
  yy <- alpha + gamma * variability_seq[i]
  comp_variability[[i]] <- yy
}

#rm
rm(i, j, xx, yy)


















#VARIABILITY (sd) --- graphics

### first, check how we're doing the modifications and their effect
# new data.frame with density modifications of sigma
n <- nrow(sigma) * ncol(sigma)
densities <- data.frame("sigmas" = as.vector(as.matrix(sigma)), "modif_sd" = as.factor(rep(1, n)))
modif_sd <- c(seq(from = 25, to = 100, by = 25)) #change as needed
for(i in modif_sd){
  data <- data.frame("sigmas" = as.vector(as.matrix(sigma * i)), "modif_sd" = as.factor(rep(i, n)))
  densities <- rbind(densities, data)
}

# requires 'sm' package
if(!require(sm)){
  install.packages("sm")
}
library(sm)

# create value labels 
sd.f <- factor(densities$sigmas, levels = levels(densities$modif_sd), labels = as.character(levels(densities$modif_sd)))

# plot densities 
sm.density.compare(densities$sigmas, densities$modif_sd, xlab = "Sigma values")
title(main = "Sigma values by sd modifications")

# add legend via mouse click
colfill <- c(2:(2 + length(levels(sd.f)))) 
legend(locator(1), levels(sd.f), fill = colfill, title = "sd", cex = 1) #click on the plot to place the legend

#plot(densities$sigmas, as.vector(densities$modif_sd), xlim = c(-250, 250), ylim = c(-5, 110),
#     xlab = "Sigma", ylab = "Standard deviation factor", main = "Variability of sigma values") #a different point of view

# rm
rm(n, densities, modif_sd, i, data, sd.f, colfill)





