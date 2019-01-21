alpha <- read.table("Results/output_glinternet/alpha.txt", header = TRUE, sep = "\t") * (-1)
sigma <- read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1)
lambda <- read.table("Results/output_glinternet/lambda.txt", header = TRUE, sep = "\t")
gamma1 <- read.table("Results/output_glinternet/gamma1.txt", header = TRUE, sep = "\t") * (-1)
gamma2 <- read.table("Results/output_glinternet/gamma2.txt", header = TRUE, sep = "\t") * (-1)
gamma3 <- read.table("Results/output_glinternet/gamma3.txt", header = TRUE, sep = "\t") * (-1)
gamma4 <- read.table("Results/output_glinternet/gamma4.txt", header = TRUE, sep = "\t") * (-1)
gamma5 <- read.table("Results/output_glinternet/gamma5.txt", header = TRUE, sep = "\t") * (-1)
gamma6 <- read.table("Results/output_glinternet/gamma6.txt", header = TRUE, sep = "\t") * (-1)


#VARIABILITY OF INTERACTION STRENGTH'S MODIFICATIONS
strength_seq <- c((1/3), 3, 9, 27) #change as needed

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


####################################################################################

# OUTPUT WITH PARALLEL COMPUTING
library(foreach)
library(doParallel)
library(parallel)
library(MASS)
source("Code/functions_structural_coex_outputs.R")
#source("Code/toolbox_coexistence.R") #included in the previous source


numCores <- detectCores() #number of cores in your computer
registerDoParallel(numCores) #important: define parallel computation to numCores

start_time <- Sys.time()
vis <- foreach (i = 1:length(comp_strength), .combine = rbind, .packages='mvtnorm') %dopar% {
  structural_coex_3spp(alpha = comp_strength[[i]], intrinsic = intr_strength[[i]])
}
end_time <- Sys.time()
end_time - start_time #time expent running code

stopImplicitCluster() #end parallel computation

####################################################################################

for (i in 1:length(strength_seq)){
  elem <- rep(strength_seq[i], (nrow(vis)/length(strength_seq)))
  if (i == 1){
    strength_values <- elem
  } else {
    strength_values <- c(strength_values, elem)
  }
}
#output_strength_3spp$strength <- strength_values
vis$VIS <- strength_values
write.table(vis, file = "Results/structural_coex_VIS_2.txt", sep = "\t", row.names = FALSE)

####################################################################################

interruptor <- FALSE #switch to TRUE if you need visualisation

if(isTRUE(interruptor)){

### Quick plots and summaries --- VIS
plot(vis$Omega, vis$theta, xlab = "SND (Omega)", ylab = "SFD (theta)", main = "Coexistence determinants, 3 spp. combinations")

#Omega ~ VIS
plot(vis$VIS, vis$Omega, xlab = "Variation of the Interaction Strength (VIS)", ylab = "SND (Omega)", main = "Omega ~ VIS")
abline(lm(vis$Omega ~ vis$VIS))
summary(lm(vis$Omega ~ vis$VIS))

plot(vis$VIS, vis$Omega, xlab = "Variation of the Interaction Strength (VIS)", ylab = "Omega", main = "Omega ~ VIS")
g1 = glm(Omega ~ VIS, family = gaussian(link = "identity"), data = vis) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
curve(predict(g1, data.frame(VIS = x), type = "resp"), add = TRUE) # draws a curve based on prediction from logistic regression model
points(vis$VIS, fitted(g1), pch=20) # optional: you could skip this draws an invisible set of points of body size survival based on a 'fit' to glm model. pch= changes type of dots.

#theta ~ VIS
plot(vis$VIS, vis$theta, xlab = "Variation of the Interaction Strength (VIS)", ylab = "SFD (theta)", main = "theta ~ VIS")
abline(lm(vis$theta ~ vis$VIS))
summary(lm(vis$theta ~ vis$VIS))

plot(vis$VIS, vis$theta, xlab = "Variation of the Interaction Strength (VIS)", ylab = "theta", main = "theta ~ VIS")
g1 = glm(theta ~ VIS, family = gaussian(link = "identity"), data = vis) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
curve(predict(g1, data.frame(VIS = x), type = "resp"), add = TRUE) # draws a curve based on prediction from logistic regression model
points(vis$VIS, fitted(g1), pch=20) # optional: you could skip this draws an invisible set of points of body size survival based on a 'fit' to glm model. pch= changes type of dots.

#differential ~ VIS
plot(vis$VIS, vis$differential, xlab = "Variation of the Interaction Strength (VIS)", ylab = "differential", main = "differential ~ VIS")
abline(lm(vis$differential ~ vis$VIS))
summary(lm(vis$differential ~ vis$VIS))

plot(vis$VIS, vis$differential, xlab = "Variation of the Interaction Strength (VIS)", ylab = "differential", main = "differential ~ VIS")
g1 = glm(differential ~ VIS, family = gaussian(link = "identity"), data = vis) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
curve(predict(g1, data.frame(VIS = x), type = "resp"), add = TRUE) # draws a curve based on prediction from logistic regression model
points(vis$VIS, fitted(g1), pch=20) # optional: you could skip this draws an invisible set of points of body size survival based on a 'fit' to glm model. pch= changes type of dots.

#overlap ~ VIS
plot(vis$VIS, vis$overlap, xlab = "Variation of the Interaction Strength (VIS)", ylab = "overlap", main = "overlap ~ VIS")
abline(lm(vis$overlap ~ vis$VIS))
summary(lm(vis$overlap ~ vis$VIS))

plot(vis$VIS, vis$overlap, xlab = "Variation of the Interaction Strength (VIS)", ylab = "overlap", main = "overlap ~ VIS")
g1 = glm(overlap ~ VIS, family = gaussian(link = "identity"), data = vis) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
curve(predict(g1, data.frame(VIS = x), type = "resp"), add = TRUE) # draws a curve based on prediction from logistic regression model
points(vis$VIS, fitted(g1), pch=20) # optional: you could skip this draws an invisible set of points of body size survival based on a 'fit' to glm model. pch= changes type of dots.

#feasibility logistic regression:
plot(vis$VIS, vis$feasibility, xlab = "Variation of the Interaction Strength (VIS)", ylab = "feasibility", main = "feasibility ~ VIS")
g1 = glm(feasibility ~ VIS, family = binomial(link='logit'), data = vis) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
curve(predict(g1, data.frame(VIS = x), type = "resp"), add = TRUE) # draws a curve based on prediction from logistic regression model
points(vis$VIS, fitted(g1), pch=20) # optional: you could skip this draws an invisible set of points of body size survival based on a 'fit' to glm model. pch= changes type of dots.

#coex_rate logistic regression:
plot(vis$VIS, vis$coex_rate, xlab = "Variation of the Interaction Strength (VIS)", ylab = "coexistence rate", main = "coexistence rate ~ VIS")
g1 = glm(coex_rate ~ VIS, family = gaussian(link = "identity"), data = vis) # run a logistic regression model (in this case, generalized linear model with logit link). see ?glm
curve(predict(g1, data.frame(VIS = x), type = "resp"), add = TRUE) # draws a curve based on prediction from logistic regression model
points(vis$VIS, fitted(g1), pch=20) # optional: you could skip this draws an invisible set of points of body size survival based on a 'fit' to glm model. pch= changes type of dots.

} #end interruptor



##VARIABILITY (sd) --- graphics
#
#### first, check how we're doing the modifications and their effect
## new data.frame with density modifications of sigma
#n <- nrow(sigma) * ncol(sigma)
#densities <- data.frame("sigmas" = as.vector(as.matrix(sigma)), "modif_sd" = as.factor(rep(1, n)))
##modif_sd <- c(seq(from = 0.2, to = 1, by = 0.2), seq(from = 2, to = 5, by = 1), 20) #change as needed
#modif_sd <- c(0.25, 0.5, 1, 2, 5, 10, 20, 40, 60, 80, 100) #change as needed
#for(i in modif_sd){
#  data <- data.frame("sigmas" = as.vector(as.matrix(sigma * i)), "modif_sd" = as.factor(rep(i, n)))
#  densities <- rbind(densities, data)
#}
#
## requires 'sm' package
#if(!require(sm)){
#  install.packages("sm")
#}
#library(sm)
#
## create value labels 
#sd.f <- factor(densities$sigmas, levels = levels(densities$modif_sd), labels = as.character(levels(densities$modif_sd)))
#
## plot densities 
#sm.density.compare(densities$sigmas, densities$modif_sd, xlab = "Sigma values")
#title(main = "Sigma values by sd modifications")
#
## add legend via mouse click
#colfill <- c(2:(2 + length(levels(sd.f)))) 
#legend(locator(1), levels(sd.f), fill = colfill, title = "sd", cex = 0.75) #click on the plot to place the legend
#
##plot(densities$sigmas, as.vector(densities$modif_sd), xlim = c(-250, 250), ylim = c(-5, 110),
##     xlab = "Sigma", ylab = "Standard deviation factor", main = "Variability of sigma values") #a different point of view
#
## rm
#rm(n, densities, modif_sd, i, data, sd.f, colfill)
#modif_sd

