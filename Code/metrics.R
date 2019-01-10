alpha <- read.table("Results/output_glinternet/alpha.txt", header = TRUE, sep = "\t") * (-1)
sigma <- read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1)
lambda <- read.table("Results/output_glinternet/lambda.txt", header = TRUE, sep = "\t")
gamma1 <- read.table("Results/output_glinternet/gamma1.txt", header = TRUE, sep = "\t") * (-1)
gamma2 <- read.table("Results/output_glinternet/gamma2.txt", header = TRUE, sep = "\t") * (-1)
gamma3 <- read.table("Results/output_glinternet/gamma3.txt", header = TRUE, sep = "\t") * (-1)
gamma4 <- read.table("Results/output_glinternet/gamma4.txt", header = TRUE, sep = "\t") * (-1)
gamma5 <- read.table("Results/output_glinternet/gamma5.txt", header = TRUE, sep = "\t") * (-1)
gamma6 <- read.table("Results/output_glinternet/gamma6.txt", header = TRUE, sep = "\t") * (-1)

#intrinsic --- from 0.05 to 1 by 0.05 (strength)
intr_strength_seq <- seq(from = 0.25, to = 1, by = 0.25) #modif
intr_strength <- list()
xx <- matrix(nrow = 22, ncol = 1)
colnames(xx) <- c("intr")
rownames(xx) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                  "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")
for (i in 1:length(intr_strength_seq)){
  for (j in 1:nrow(xx)){
    xx[j, 1] <- lambda[j,] + rowSums(sigma[j,]) * intr_strength_seq[i]
  }
  intr_strength[[i]] <- xx
}

#alphas --- from 0.01 to 1 by 0.01 (strength)
gamma <- gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6
comp_strength_seq <- seq(from = 0.25, to = 1, by = 0.25) #modif
comp_strength <- list()
for (i in 1:length(comp_strength_seq)){
  yy <- alpha + gamma * comp_strength_seq[i]
  comp_strength[[i]] <- yy
}

#rm
rm(i, j, xx, yy)


#############

### first try --- strength modification
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










#play before leaving
start_time <- Sys.time()
alpha <- comp_strength[[1]]
intrinsic <- intr_strength[[1]]
intento <- structural_coex_3spp(alpha, intrinsic)
for (i in 2:length(comp_strength)){
  alpha <- comp_strength[[i]]
  intrinsic <- intr_strength[[i]]
  intento <- rbind(intento, structural_coex_3spp(alpha, intrinsic))
}
end_time <- Sys.time()
end_time - start_time #Time expent running code








