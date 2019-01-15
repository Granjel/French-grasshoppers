################################
### OPTIM WITH GRASSHOPPERS: ###
################################

### predictive functions depending on the number of plant and grasshopper species:

f_ACHMIL_gh <- function(n, theta, Y, X_plant, X_gh){
  lambda <- theta[1] #max cover without competition
  alpha <- t(t(rep(theta[2], 39))) #as many alphas as plant species
  gamma1 <- rep(theta[3], n) #comp. coef. for plant-plant with gh1
  gamma2 <- rep(theta[4], n) #comp. coef. for plant-plant with gh2
  gamma3 <- rep(theta[5], n) #comp. coef. for plant-plant with gh3
  gamma4 <- rep(theta[6], n) #comp. coef. for plant-plant with gh4
  gamma5 <- rep(theta[7], n) #comp. coef. for plant-plant with gh5
  gamma6 <- rep(theta[8], n) #comp. coef. for plant-plant with gh6
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) -
    log(1 + ((X_gh[, 1] * gamma1) %*% t(X_plant))) - log(1 + ((X_gh[, 2] * gamma2) %*% t(X_plant))) - log(1 + ((X_gh[, 3] * gamma3) %*% t(X_plant))) -
    log(1 + ((X_gh[, 4] * gamma4) %*% t(X_plant))) - log(1 + ((X_gh[, 5] * gamma5) %*% t(X_plant))) - log(1 + ((X_gh[, 5] * gamma5) %*% t(X_plant)))
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}


### load dataset:

d <- read.table("Data_Fg/FG.txt", header = TRUE, sep = "\t")


### subsample for different dates

d2 <- d[d$time == "2",] #June 2012 ### changing this command changes everything !!!
summary(d2$Focal) #zero ANTODO, GERDIS, TRIFLA and VERPER


### one dataset for each focal species

d_ACHMIL_gh <- d2[d2$Focal == "ACHMIL",]
n_ACHMIL_gh <- nrow(d_ACHMIL_gh)

### X_plant, X_gh and Y matrices for each plant species:

X_ACHMIL_plant <- as.matrix(d_ACHMIL_gh[, seq(15, 91, by = 2)]) #competition
X_ACHMIL_gh <- as.matrix(d_ACHMIL_gh[6:11]) #depends on the number of links for every species (suppl. mat. Gross)
Y_ACHMIL_gh <- d_ACHMIL_gh$Cover #cover


### functions for optim, species-specific

f_ACHMIL_gh_optim <- function(theta){f_ACHMIL_gh(n_ACHMIL_gh, theta, Y_ACHMIL_gh, X_ACHMIL_plant, X_ACHMIL_gh)}


### optim wrap for the diff. species (one 'ini' and 'lower' for each spp.) ---:

ini_ACHMIL_gh <- c(1, rep(0.000001, 274))
lower_ACHMIL_gh <- rep(0.000001, 274)
out_ACHMIL_gh <- optim(ini_ACHMIL_gh, f_ACHMIL_gh_optim, lower = lower_ACHMIL_gh, method = 'L-BFGS-B', hessian = T)












