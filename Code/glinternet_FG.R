#####################################################
### 'glinternet' analysis --- French grasshoppers ###
#####################################################
### Rodrigo R. Granjel --- December 2018 ############
#####################################################

### packages required:
#install.packages("glinternet")
require(glinternet)
require(beepr)
#library(glinternet)

### load dataset:
d <- read.table("Data_Fg/FG.txt", header = TRUE, sep = "\t")

### subsample for different dates
#d1 <- d[d$time == "1",] #Jun 2012
d2 <- d[d$time == "2",] #Sep 2012
#d3 <- d[d$time == "3",] #May 2013
#d4 <- d[d$time == "4",] #Jul 2013
d5 <- d[d$time == "5",] #Oct 2013
#d6 <- d[d$time == "6",] #May 2014
dh <- rbind(d2, d5) #bind of times after grasshopper grazing
data <- dh #WORKING DATASET (change as needed)
rm(d, d2, d5, dh)

### number of plant species, grashopper species and lambdas
plants <- 36
grasshoppers <- 6
n <- 204 #for the 'glinternet' function, optimised for the species

### define outputs (lambda, alpha, sigma, gamma1, gamma2, gamma3, gamma4, gamma5, and gamma6)
lambda <- NULL #vector of intercepts
alpha <- as.data.frame(matrix(NA, nrow = 0, ncol = plants)) #plant-plant matrix (no grasshoppers)
sigma <- as.data.frame(matrix(NA, nrow = 0, ncol = grasshoppers)) #plant-grasshopper matrix
gamma1 <- as.data.frame(matrix(NA, nrow = 0, ncol = plants)) #plant-plant matrix (grasshopper 1)
gamma2 <- as.data.frame(matrix(NA, nrow = 0, ncol = plants)) #plant-plant matrix (grasshopper 2)
gamma3 <- as.data.frame(matrix(NA, nrow = 0, ncol = plants)) #plant-plant matrix (grasshopper 3)
gamma4 <- as.data.frame(matrix(NA, nrow = 0, ncol = plants)) #plant-plant matrix (grasshopper 4)
gamma5 <- as.data.frame(matrix(NA, nrow = 0, ncol = plants)) #plant-plant matrix (grasshopper 5)
gamma6 <- as.data.frame(matrix(NA, nrow = 0, ncol = plants)) #plant-plant matrix (grasshopper 6)

#APPLY 'glinternet' TO DIFF PLANTS:

### ACHMIL #######################################################################################################################

### specific dataset
d_ACHMIL <- data[data$Focal == "ACHMIL",]

### X_PLANT and Y_PLANT matrices
X_ACHMIL <- as.data.frame(cbind(as.matrix(d_ACHMIL[, seq(15, 91, by = 2)]), as.matrix(d_ACHMIL[6:11]))) #explanative variables
Y_ACHMIL <- d_ACHMIL$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_ACHMIL <- glinternet(X_ACHMIL, Y_ACHMIL, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_ACHMIL #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(1e-15, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_ACHMIL, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_ACHMIL, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_ACHMIL, MODEL)



### ARRELA #######################################################################################################################

### specific dataset
d_ARRELA <- data[data$Focal == "ARRELA",]

### X_PLANT and Y_PLANT matrices
X_ARRELA <- as.data.frame(cbind(as.matrix(d_ARRELA[, seq(15, 91, by = 2)]), as.matrix(d_ARRELA[6:11]))) #explanative variables
Y_ARRELA <- d_ARRELA$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_ARRELA <- glinternet(X_ARRELA, Y_ARRELA, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_ARRELA #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(1e-15, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_ARRELA, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_ARRELA, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_ARRELA, MODEL)



### BROERE #######################################################################################################################

### specific dataset
d_BROERE <- data[data$Focal == "BROERE",]

### X_PLANT and Y_PLANT matrices
X_BROERE <- as.data.frame(cbind(as.matrix(d_BROERE[, seq(15, 91, by = 2)]), as.matrix(d_BROERE[6:11]))) #explanative variables
Y_BROERE <- d_BROERE$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_BROERE <- glinternet(X_BROERE, Y_BROERE, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_BROERE #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(1e-15, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_BROERE, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_BROERE, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_BROERE, MODEL)



### CENJAC #######################################################################################################################

### specific dataset
d_CENJAC <- data[data$Focal == "CENJAC",]

### X_PLANT and Y_PLANT matrices
X_CENJAC <- as.data.frame(cbind(as.matrix(d_CENJAC[, seq(15, 91, by = 2)]), as.matrix(d_CENJAC[6:11]))) #explanative variables
Y_CENJAC <- d_CENJAC$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_CENJAC <- glinternet(X_CENJAC, Y_CENJAC, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_CENJAC #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(1e-15, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_CENJAC, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_CENJAC, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_CENJAC, MODEL)



### CONARV #######################################################################################################################

### specific dataset
d_CONARV <- data[data$Focal == "CONARV",]

### X_PLANT and Y_PLANT matrices
X_CONARV <- as.data.frame(cbind(as.matrix(d_CONARV[, seq(15, 91, by = 2)]), as.matrix(d_CONARV[6:11]))) #explanative variables
Y_CONARV <- d_CONARV$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_CONARV <- glinternet(X_CONARV, Y_CONARV, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_CONARV #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(1e-15, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_CONARV, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_CONARV, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_CONARV, MODEL)



### CREPIS #######################################################################################################################

### specific dataset
d_CREPIS <- data[data$Focal == "CREPIS",]

### X_PLANT and Y_PLANT matrices
X_CREPIS <- as.data.frame(cbind(as.matrix(d_CREPIS[, seq(15, 91, by = 2)]), as.matrix(d_CREPIS[6:11]))) #explanative variables
Y_CREPIS <- d_CREPIS$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_CREPIS <- glinternet(X_CREPIS, Y_CREPIS, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_CREPIS #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(1e-15, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_CREPIS, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_CREPIS, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_CREPIS, MODEL)


### DACGLO #######################################################################################################################

### specific dataset
d_DACGLO <- data[data$Focal == "DACGLO",]

### X_PLANT and Y_PLANT matrices
X_DACGLO <- as.data.frame(cbind(as.matrix(d_DACGLO[, seq(15, 91, by = 2)]), as.matrix(d_DACGLO[6:11]))) #explanative variables
Y_DACGLO <- d_DACGLO$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_DACGLO <- glinternet(X_DACGLO, Y_DACGLO, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_DACGLO #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_DACGLO, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_DACGLO, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_DACGLO, MODEL)



### DAUCAR #######################################################################################################################

### specific dataset
d_DAUCAR <- data[data$Focal == "DAUCAR",]

### X_PLANT and Y_PLANT matrices
X_DAUCAR <- as.data.frame(cbind(as.matrix(d_DAUCAR[, seq(15, 91, by = 2)]), as.matrix(d_DAUCAR[6:11]))) #explanative variables
Y_DAUCAR <- d_DAUCAR$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_DAUCAR <- glinternet(X_DAUCAR, Y_DAUCAR, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_DAUCAR #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_DAUCAR, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_DAUCAR, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_DAUCAR, MODEL)



### ELYREP #######################################################################################################################

### specific dataset
d_ELYREP <- data[data$Focal == "ELYREP",]

### X_PLANT and Y_PLANT matrices
X_ELYREP <- as.data.frame(cbind(as.matrix(d_ELYREP[, seq(15, 91, by = 2)]), as.matrix(d_ELYREP[6:11]))) #explanative variables
Y_ELYREP <- d_ELYREP$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_ELYREP <- glinternet(X_ELYREP, Y_ELYREP, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_ELYREP #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_ELYREP, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_ELYREP, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_ELYREP, MODEL)



### ERYNGE #######################################################################################################################

### specific dataset
d_ERYNGE <- data[data$Focal == "ERYNGE",]

### X_PLANT and Y_PLANT matrices
X_ERYNGE <- as.data.frame(cbind(as.matrix(d_ERYNGE[, seq(15, 91, by = 2)]), as.matrix(d_ERYNGE[6:11]))) #explanative variables
Y_ERYNGE <- d_ERYNGE$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_ERYNGE <- glinternet(X_ERYNGE, Y_ERYNGE, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_ERYNGE #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_ERYNGE, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_ERYNGE, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_ERYNGE, MODEL)



### FESRUB #######################################################################################################################

### specific dataset
d_FESRUB <- data[data$Focal == "FESRUB",]

### X_PLANT and Y_PLANT matrices
X_FESRUB <- as.data.frame(cbind(as.matrix(d_FESRUB[, seq(15, 91, by = 2)]), as.matrix(d_FESRUB[6:11]))) #explanative variables
Y_FESRUB <- d_FESRUB$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_FESRUB <- glinternet(X_FESRUB, Y_FESRUB, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_FESRUB #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_FESRUB, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_FESRUB, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_FESRUB, MODEL)



### GALVER #######################################################################################################################

### specific dataset
d_GALVER <- data[data$Focal == "GALVER",]

### X_PLANT and Y_PLANT matrices
X_GALVER <- as.data.frame(cbind(as.matrix(d_GALVER[, seq(15, 91, by = 2)]), as.matrix(d_GALVER[6:11]))) #explanative variables
Y_GALVER <- d_GALVER$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_GALVER <- glinternet(X_GALVER, Y_GALVER, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_GALVER #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_GALVER, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_GALVER, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_GALVER, MODEL)



### LEUVUL #######################################################################################################################

### specific dataset
d_LEUVUL <- data[data$Focal == "LEUVUL",]

### X_PLANT and Y_PLANT matrices
X_LEUVUL <- as.data.frame(cbind(as.matrix(d_LEUVUL[, seq(15, 91, by = 2)]), as.matrix(d_LEUVUL[6:11]))) #explanative variables
Y_LEUVUL <- d_LEUVUL$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_LEUVUL <- glinternet(X_LEUVUL, Y_LEUVUL, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_LEUVUL #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_LEUVUL, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_LEUVUL, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_LEUVUL, MODEL)



### ONOREP #######################################################################################################################

### specific dataset
d_ONOREP <- data[data$Focal == "ONOREP",]

### X_PLANT and Y_PLANT matrices
X_ONOREP <- as.data.frame(cbind(as.matrix(d_ONOREP[, seq(15, 91, by = 2)]), as.matrix(d_ONOREP[6:11]))) #explanative variables
Y_ONOREP <- d_ONOREP$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_ONOREP <- glinternet(X_ONOREP, Y_ONOREP, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_ONOREP #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_ONOREP, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_ONOREP, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_ONOREP, MODEL)



### PICECH #######################################################################################################################

### specific dataset
d_PICECH <- data[data$Focal == "PICECH",]

### X_PLANT and Y_PLANT matrices
X_PICECH <- as.data.frame(cbind(as.matrix(d_PICECH[, seq(15, 91, by = 2)]), as.matrix(d_PICECH[6:11]))) #explanative variables
Y_PICECH <- d_PICECH$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_PICECH <- glinternet(X_PICECH, Y_PICECH, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_PICECH #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_PICECH, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_PICECH, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_PICECH, MODEL)



### PLALAN #######################################################################################################################

### specific dataset
d_PLALAN <- data[data$Focal == "PLALAN",]

### X_PLANT and Y_PLANT matrices
X_PLALAN <- as.data.frame(cbind(as.matrix(d_PLALAN[, seq(15, 91, by = 2)]), as.matrix(d_PLALAN[6:11]))) #explanative variables
Y_PLALAN <- d_PLALAN$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_PLALAN <- glinternet(X_PLALAN, Y_PLALAN, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_PLALAN #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_PLALAN, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_PLALAN, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_PLALAN, MODEL)



### POAANG #######################################################################################################################

### specific dataset
d_POAANG <- data[data$Focal == "POAANG",]

### X_PLANT and Y_PLANT matrices
X_POAANG <- as.data.frame(cbind(as.matrix(d_POAANG[, seq(15, 91, by = 2)]), as.matrix(d_POAANG[6:11]))) #explanative variables
Y_POAANG <- d_POAANG$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_POAANG <- glinternet(X_POAANG, Y_POAANG, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_POAANG #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_POAANG, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_POAANG, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_POAANG, MODEL)



### RANACR #######################################################################################################################

### specific dataset
d_RANACR <- data[data$Focal == "RANACR",]

### X_PLANT and Y_PLANT matrices
X_RANACR <- as.data.frame(cbind(as.matrix(d_RANACR[, seq(15, 91, by = 2)]), as.matrix(d_RANACR[6:11]))) #explanative variables
Y_RANACR <- d_RANACR$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_RANACR <- glinternet(X_RANACR, Y_RANACR, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_RANACR #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_RANACR, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_RANACR, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_RANACR, MODEL)



### RUMACE #######################################################################################################################

### specific dataset
d_RUMACE <- data[data$Focal == "RUMACE",]

### X_PLANT and Y_PLANT matrices
X_RUMACE <- as.data.frame(cbind(as.matrix(d_RUMACE[, seq(15, 91, by = 2)]), as.matrix(d_RUMACE[6:11]))) #explanative variables
Y_RUMACE <- d_RUMACE$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_RUMACE <- glinternet(X_RUMACE, Y_RUMACE, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_RUMACE #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_RUMACE, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_RUMACE, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_RUMACE, MODEL)



### SALPRA #######################################################################################################################

### specific dataset
d_SALPRA <- data[data$Focal == "SALPRA",]

### X_PLANT and Y_PLANT matrices
X_SALPRA <- as.data.frame(cbind(as.matrix(d_SALPRA[, seq(15, 91, by = 2)]), as.matrix(d_SALPRA[6:11]))) #explanative variables
Y_SALPRA <- d_SALPRA$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_SALPRA <- glinternet(X_SALPRA, Y_SALPRA, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_SALPRA #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_SALPRA, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_SALPRA, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_SALPRA, MODEL)



### TAROFF #######################################################################################################################

### specific dataset
d_TAROFF <- data[data$Focal == "TAROFF",]

### X_PLANT and Y_PLANT matrices
X_TAROFF <- as.data.frame(cbind(as.matrix(d_TAROFF[, seq(15, 91, by = 2)]), as.matrix(d_TAROFF[6:11]))) #explanative variables
Y_TAROFF <- d_TAROFF$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_TAROFF <- glinternet(X_TAROFF, Y_TAROFF, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_TAROFF #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_TAROFF, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_TAROFF, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_TAROFF, MODEL)



### TRIPRA #######################################################################################################################

### specific dataset
d_TRIPRA <- data[data$Focal == "TRIPRA",]

### X_PLANT and Y_PLANT matrices
X_TRIPRA <- as.data.frame(cbind(as.matrix(d_TRIPRA[, seq(15, 91, by = 2)]), as.matrix(d_TRIPRA[6:11]))) #explanative variables
Y_TRIPRA <- d_TRIPRA$Cover #cover; dependent variable

### 'glinternet' MODEL
gli_TRIPRA <- glinternet(X_TRIPRA, Y_TRIPRA, numLevels = rep(1, 45), lambda = NULL, nLambda = n, lambdaMinRatio = 0.01,
                         interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian",
                         tol = 1e-05, maxIter = 5000, verbose = TRUE, numCores = 1)
MODEL <- gli_TRIPRA #for future proceedings

### which lambda should be chosen? phi
bet <- as.data.frame(matrix(NA, nrow = n, ncol = 2))
colnames(bet) <- c("n_var", "n_lambda")
for(i in 1:n){
  bet[i, 1] <- length(MODEL$betahat[[i]])
  bet[i, 2] <- i
}
plot(bet$n_lambda, bet$n_var)
phi <- min(which(bet == max(bet$n_var))) #lambda value selected (max number of variables)

### save fix and interaction coefficients
fix <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$`mainEffects`$cont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$mainEffectsCoef$`cont`))
colnames(fix) <- c("pos", "coef")
fix$pos <- as.numeric(fix$pos)
fix$coef <- as.numeric(fix$coef)

int <- as.data.frame(cbind(coef(MODEL, lambdaIndex = phi)[[1]]$interactions$contcont,
                           coef(MODEL, lambdaIndex = phi)[[1]]$interactionsCoef$`contcont`))
colnames(int) <- c("pos_i", "pos_j", "coef")
int$pos_i <- as.numeric(int$pos_i)
int$pos_j <- as.numeric(int$pos_j)
int$coef <- as.numeric(int$coef)
int <- int[int$pos_i <= 36,]

### order fix dataframe
fix <- fix[order(fix$pos),]
rownames(fix) <- seq(length = nrow(fix))

### plant-plant without grasshopper effects
pp_fix <- fix[fix$pos <= 36,]
pg_fix <- fix[fix$pos >= 40,]

pp_spp <- rep(1e-15, plants)
for (i in 1:nrow(pp_fix)){
  pp_spp[pp_fix$pos[i]] <- pp_fix$coef[i]
}

### plant-grasshopper
pg_spp <- rep(NA, grasshoppers)
for (i in 1:nrow(pg_fix)){
  pg_spp[pg_fix$pos[i]-39] <- pg_fix$coef[i] #trick: removing 39 positions
}

### plant-plant with grasshopper effects
gh1_int <- int[int$pos_j == 40,] #subset only with gh1 (pos. 40)
ppg1_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh1_int)){
  ppg1_spp[gh1_int$pos_i[i]] <- gh1_int$coef[i]
}

gh2_int <- int[int$pos_j == 41,] #subset only with gh2 (pos. 41)
ppg2_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh2_int)){
  ppg2_spp[gh2_int$pos_i[i]] <- gh2_int$coef[i]
}

gh3_int <- int[int$pos_j == 42,] #subset only with gh3 (pos. 42)
ppg3_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh3_int)){
  ppg3_spp[gh3_int$pos_i[i]] <- gh3_int$coef[i]
}

gh4_int <- int[int$pos_j == 43,] #subset only with gh4 (pos. 43)
ppg4_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh4_int)){
  ppg4_spp[gh4_int$pos_i[i]] <- gh4_int$coef[i]
}

gh5_int <- int[int$pos_j == 44,] #subset only with gh5 (pos. 44)
ppg5_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh5_int)){
  ppg5_spp[gh5_int$pos_i[i]] <- gh5_int$coef[i]
}

gh6_int <- int[int$pos_j == 45,] #subset only with gh6 (pos. 45)
ppg6_spp <- rep(1e-15, plants)
for (i in 1:nrow(gh6_int)){
  ppg6_spp[gh6_int$pos_i[i]] <- gh6_int$coef[i]
}

### saving outputs
lambda <- c(lambda, MODEL$betahat[[phi]][1])
alpha <- rbind(alpha, pp_spp)
sigma <- rbind(sigma, pg_spp)
gamma1 <- rbind(gamma1, ppg1_spp)
gamma2 <- rbind(gamma2, ppg2_spp)
gamma3 <- rbind(gamma3, ppg3_spp)
gamma4 <- rbind(gamma4, ppg4_spp)
gamma5 <- rbind(gamma5, ppg5_spp)
gamma6 <- rbind(gamma6, ppg6_spp)

### remove stuff
rm(bet, d_TRIPRA, fix, gh1_int, gh2_int, gh3_int, gh4_int, gh5_int, gh6_int, int, pg_fix, pp_fix, X_TRIPRA, i,
   pg_spp, phi, pp_spp, ppg1_spp, ppg2_spp, ppg3_spp, ppg4_spp, ppg5_spp, ppg6_spp, Y_TRIPRA, MODEL)



##########################################
####### NAME AND SAVE THE MATRICES #######
##########################################

### rownames
rownames(alpha) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                     "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

rownames(sigma) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                     "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

rownames(gamma1) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

rownames(gamma2) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

rownames(gamma3) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

rownames(gamma4) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

rownames(gamma5) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

rownames(gamma6) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

#colnames
alpha <- alpha[,c(-2, -12, -15, -16, -18, -19, -20, -23, -26, -27, -31, -33, -35, -36)] #remove not desired plant species
colnames(alpha) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                     "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

colnames(sigma) <- c("Cb", "Cd", "Ci", "Ee", "Pg", "Pp")

gamma1 <- gamma1[,c(-2, -12, -15, -16, -18, -19, -20, -23, -26, -27, -31, -33, -35, -36)] #remove not desired plant species
colnames(gamma1) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

gamma2 <- gamma2[,c(-2, -12, -15, -16, -18, -19, -20, -23, -26, -27, -31, -33, -35, -36)] #remove not desired plant species
colnames(gamma2) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

gamma3 <- gamma3[,c(-2, -12, -15, -16, -18, -19, -20, -23, -26, -27, -31, -33, -35, -36)] #remove not desired plant species
colnames(gamma3) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

gamma4 <- gamma4[,c(-2, -12, -15, -16, -18, -19, -20, -23, -26, -27, -31, -33, -35, -36)] #remove not desired plant species
colnames(gamma4) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

gamma5 <- gamma5[,c(-2, -12, -15, -16, -18, -19, -20, -23, -26, -27, -31, -33, -35, -36)] #remove not desired plant species
colnames(gamma5) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

gamma6 <- gamma6[,c(-2, -12, -15, -16, -18, -19, -20, -23, -26, -27, -31, -33, -35, -36)] #remove not desired plant species
colnames(gamma6) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                      "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

### names lambda vector
names(lambda) <- c("ACHMIL", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESRUB",
                   "GALVER", "LEUVUL", "ONOREP", "PICECH", "PLALAN", "POAANG", "RANACR", "RUMACE", "SALPRA", "TAROFF", "TRIPRA")

### write the matrices and vector
write.table(lambda, file = "Results/output_glinternet/lambda.txt", sep = "\t", row.names = TRUE)
write.table(alpha, file = "Results/output_glinternet/alpha.txt", sep = "\t", row.names = TRUE)
write.table(sigma, file = "Results/output_glinternet/sigma.txt", sep = "\t", row.names = TRUE)
write.table(gamma1, file = "Results/output_glinternet/gamma1.txt", sep = "\t", row.names = TRUE)
write.table(gamma2, file = "Results/output_glinternet/gamma2.txt", sep = "\t", row.names = TRUE)
write.table(gamma3, file = "Results/output_glinternet/gamma3.txt", sep = "\t", row.names = TRUE)
write.table(gamma4, file = "Results/output_glinternet/gamma4.txt", sep = "\t", row.names = TRUE)
write.table(gamma5, file = "Results/output_glinternet/gamma5.txt", sep = "\t", row.names = TRUE)
write.table(gamma6, file = "Results/output_glinternet/gamma6.txt", sep = "\t", row.names = TRUE)

#end
beep(4)
