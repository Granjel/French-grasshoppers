##################################################
## FRENCH GRASSHOPPERS ------------------------ ##
## Optim Analyses of the 36 plant species ----- ##
## Rodrigo R. Granjel ---- 28th June 2018 ----- ##
##################################################

#load library -----
#install.packages("likelihood")
library(likelihood)

#load data -----
d <- read.table(file = "C:/Users/Granjel RR/Desktop/Nico Gross/FG.txt", header = TRUE, sep = "\t")
# d <- read.table(file = "/Users/oscargodoy/Dropbox/French grasshoppers/FG.txt", header = TRUE, sep = "\t")
total <- rep(NA, nrow(d))
for (i in 1:nrow(d)){
  if(d$datapoint[i] == 1 || d$datapoint[i] == 3 || d$datapoint[i] == 7 || d$datapoint[i] == 9){
    total[i] <- (sum(d[i, seq(15, 91,by = 2)])/3)
  } else {
    if(d$datapoint[i] == 2 || d$datapoint[i] == 4 || d$datapoint[i] == 6 || d$datapoint[i] == 8){
      total[i] <- (sum(d[i, seq(15, 91,by = 2)])/4)
    } else {
      total[i] <- (sum(d[i, seq(15, 91, by = 2)])/5)
    }
  }
  print((i/nrow(d))*100)
}
d <- cbind(d, "total" = total)

# subset of databases per plant -----
d_achmil <- subset(d, Focal == "ACHMIL")
d_antodo <- subset(d, Focal == "ANTODO")
d_arrela <- subset(d, Focal == "ARRELA")
d_broere <- subset(d, Focal == "BROERE")
d_cenjac <- subset(d, Focal == "CENJAC")
d_conarv <- subset(d, Focal == "CONARV")
d_crepis <- subset(d, Focal == "CREPIS")
d_dacglo <- subset(d, Focal == "DACGLO")
d_daucar <- subset(d, Focal == "DAUCAR")
d_elyrep <- subset(d, Focal == "ELYREP")
d_erynge <- subset(d, Focal == "ERYNGE")
d_fesaru <- subset(d, Focal == "FESARU")
d_fesrub <- subset(d, Focal == "FESRUB")
d_galver <- subset(d, Focal == "GALVER")
d_gerdis <- subset(d, Focal == "GERDIS")
d_gerrot <- subset(d, Focal == "GERROT")
d_leuvul <- subset(d, Focal == "LEUVUL")
d_lolper <- subset(d, Focal == "LOLPER")
d_lotcor <- subset(d, Focal == "LOTCOR")
d_medara <- subset(d, Focal == "MEDARA")
d_onorep <- subset(d, Focal == "ONOREP")
d_picech <- subset(d, Focal == "PICECH")
d_pichie <- subset(d, Focal == "PICHIE")
d_plalan <- subset(d, Focal == "PLALAN")
d_poaang <- subset(d, Focal == "POAANG")
d_poapra <- subset(d, Focal == "POAPRA")
d_poatri <- subset(d, Focal == "POATRI")
d_ranacr <- subset(d, Focal == "RANACR")
d_rumace <- subset(d, Focal == "RUMACE")
d_salpra <- subset(d, Focal == "SALPRA")
d_sonchu <- subset(d, Focal == "SONCHU")
d_taroff <- subset(d, Focal == "TAROFF")
d_trifla <- subset(d, Focal == "TRIFLA")
d_tripra <- subset(d, Focal == "TRIPRA")
d_verbof <- subset(d, Focal == "VERBOF")
d_verper <- subset(d, Focal == "VERPER")

## get a list of target species to work through sequentially:
splist <- as.vector(levels(d$Focal))

## objects to hold the final parameter estimates from model 3 and 5:
alpha_matrix <- matrix(0, nrow = length(splist), ncol = length(splist))
row.names(alpha_matrix) <- splist
colnames(alpha_matrix) <- splist
matrix<-matrix(0, nrow = length(splist), ncol = 15)
colnames(matrix) <- c("lambda", "theta", "eta", "zeta", "kappa", "gamma", "iota", "common_alpha", "omega", "psi", "epsilon", "tau", "omicron", "rho", "sigma")
row.names(matrix) <- splist

#ACHMIL -----
#model 1 - no effect of density (no competitive effects)
compmodel1 <- function(par){
  
  ## mu (=lambda later) and sigma parameters for the normal distribution (assuming lognormal error- cover data are logged)
  lambda <- par[1]
  sigma <- par[2]
  
  #this the predictive model- here is just fitting a horizontal line through the data:
  pred <- rep(lambda, times = length(log_cover))
  #these are the log likelihoods of the data given the model + parameters
  llik <- dnorm(log_cover, log(pred), sd = sigma, log = TRUE) 
  #return the sum of negative log likelihoods - what optim minimizes
  return(sum(-1 * llik))
}

#model 2 - competition, but no difference between species
compmodel2 <- function(par){
  lambda <- par[1] ## same as model 1
  alpha <- par[2]  ## new parameter introduced in model 2
  sigma <- par[3] ## same as model 1
  # predictive model:
  pred <- lambda/(1+alpha*(d_achmil$total))
  # log likelihoods of data given the model + parameters:
  llik <- dnorm(log_cover, log(pred), sd = sigma, log = TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1 * llik))
}

#model 3a- common effect of competition including salinity and pollinators. 
compmodel3a<-function(par){
  lambda <- par[1] ## same as model 1
  theta <- par[2]
  eta <- par[3]
  alpha <- par[4]  ## new parameter introduced in model 2
  omega <- par[5]
  psi <- par[6]
  sigma <- par[7] ## same as model 1
  # predictive model:
  pred <- lambda * (1 + theta * d_achmil$Cb + eta * d_achmil$Cd)/ (1 + (alpha + omega * d_achmil$Cb + psi * d_achmil$Cd) * d_achmil$total) 
  # log likelihoods of data given the model + parameters:
  llik <- dnorm(log_cover, log(pred), sd = sigma, log = TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1 * llik))
}

#model 3- common effect of competition including salinity and pollinators. 
compmodel3<-function(par){
  lambda <- par[1] ## same as model 1
  theta <- par[2]
  eta <- par[3]
  zeta <- par[4]
  kappa <- par[5]
  epsilon <- par[6]
  iota <- par[7]
  alpha <- par[8]  ## new parameter introduced in model 2
  omega <- par[9]
  psi <- par[10]
  epsilon <- par[11]
  tau <- par[12]
  omicron <- par[13]
  rho <- par[14]
  sigma <- par[15] ## same as model 1
  # predictive model:
  pred <- lambda * (1 + theta * d_achmil$Cb + eta * d_achmil$Cd + zeta * d_achmil$Ci + kappa * d_achmil$Ee + iota * d_achmil$Pp)/ (1 + (alpha + omega * d_achmil$Cb + psi * d_achmil$Cd + epsilon * d_achmil$Ci + tau * d_achmil$Ee + omicron * d_achmil$Pg + rho * d_achmil$Pg) * d_achmil$total) 
  # log likelihoods of data given the model + parameters:
  llik <- dnorm(log_cover, log(pred), sd = sigma, log = TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1 * llik))
}

#gamma * d_achmil$Pg
##### MODIFIED UP TO THIS POINT ---------- next model is not complicated ###Update Sept 5th this needs to be updated accordingly to include only competition between species, no grasshopper effect. 
#model 4 - all species have different competitive effects
compmodel4 <- function(par){
  lambda <- par[1] #same as model 2
  a_ACHMIL <- par[2]	## new parameters - use alpha estimate from model 2 as start value for fitting
  a_ANTODO <- par[3]
  a_ARRELA <- par[4]
  a_BROERE <- par[5]
  a_CENJAC <- par[6]
  a_CONARV <- par[7]
  a_CREPIS <- par[8]
  a_DACGLO <- par[9]
  a_DAUCAR <- par[10]
  a_ELYREP <- par[11]
  a_ERYNGE <- par[12]
  a_FESARU <- par[13]
  a_FESRUB <- par[14]
  a_GALVER <- par[15]
  a_GERDIS <- par[16]
  a_GERROT <- par[17]
  a_LEUVUL <- par[18]
  a_LOLPER <- par[19]
  a_LOTCOR <- par[20]
  a_MEDARA <- par[21]
  a_ONOREP <- par[22]
  a_PICECH <- par[23]
  a_PICHIE <- par[24]
  a_PLALAN <- par[25]
  a_POAANG <- par[26]
  a_POAPRA <- par[27]
  a_POATRI <- par[28]
  a_RANACR <- par[29]
  a_RUMACE <- par[30]
  a_SALPRA <- par[31]
  a_SONCHU <- par[32]
  a_TAROFF <- par[33]
  a_TRIFLA <- par[34]
  a_TRIPRA <- par[35]
  a_VERBOF <- par[36]
  a_VERPER <- par[37]
  a_legumes <- par[38]
  a_grasses <- par[39]
  a_other <- par[40]
  sigma <- par[41] ## same as model 2
  
  ##probably overkill to name all the alphas, but it helps keep things straight
  
  # same form as model 1, but super long given all the species:
  pred <- lambda / (1 + a_CHFU * d_chfu$CHFU + a_BEMA * d_chfu$BEMA + a_LEMA * d_chfu$LEMA + a_MEEL * d_chfu$MEEL + a_MESU * d_chfu$MESU + a_PUPA * d_chfu$PUPA + a_NON_FOCAL * d_chfu$NON_FOCAL)
  
  # likelihood as before:
  llik <- dnorm(log_cover, log(pred), sd = sigma, log = TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1 * llik)) #sum of negative log likelihoods
}


#model 5 - 
compmodel5 <- function(par){
  lambda <- par[1] 
  gamma <- par[2]
  theta <- par[3]
  a_CHFU <- par[4]
  a_BEMA <- par[5]
  a_LEMA <- par[6]
  a_MEEL <- par[7]
  a_MESU <- par[8]
  a_PUPA <- par[9]
  a_NON_FOCAL <- par[10]
  omega <- par[11]
  psi <- par[12]
  sigma <- par[13]
  
  pred <- lambda * (1 + theta * d_chfu$salinity + gamma * d_chfu$pol_sum) / 
    (1 + (a_CHFU + omega* d_chfu$pol_sum + psi* d_chfu$salinity)* d_chfu$CHFU + (a_BEMA + omega* d_chfu$pol_sum + psi* d_chfu$salinity) * d_chfu$BEMA 
     + (a_LEMA + omega* d_chfu$pol_sum + psi* d_chfu$salinity)* d_chfu$LEMA  + (a_MEEL + omega* d_chfu$pol_sum + psi* d_chfu$salinity)* d_chfu$MEEL 
     + (a_MESU + omega* d_chfu$pol_sum + psi* d_chfu$salinity)* d_chfu$MESU + (a_PUPA + omega* d_chfu$pol_sum + psi* d_chfu$salinity)* d_chfu$PUPA 
     + (a_NON_FOCAL + omega* d_chfu$pol_sum + psi* d_chfu$salinity) * d_achmil$total)
  
  
  # likelihood as before:
  llik <- dnorm(log_cover, log(pred), sd = sigma, log = TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1 * llik)) #sum of negative log likelihoods
}


#model 6 - 
compmodel6<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega_CHFU<-par[11]
  omega_BEMA<-par[12]
  omega_LEMA<-par[13]
  omega_MEEL<-par[14]
  omega_MESU<-par[15]
  omega_PUPA<-par[16]
  omega_NON_FOCAL<-par[17]
  psi_CHFU<-par[18]
  psi_BEMA<-par[19]
  psi_LEMA<-par[20]
  psi_MEEL<-par[21]
  psi_MESU<-par[22]
  psi_PUPA<-par[23]
  psi_NON_FOCAL<-par[24]
  sigma<-par[25]
  
  pred<- lambda*(1 + theta* d_chfu$salinity + gamma* d_chfu$pol_sum)/ 
    
    (1+ (a_CHFU + omega_CHFU* d_chfu$pol_sum + psi_CHFU* d_chfu$salinity)* d_chfu$CHFU + (a_BEMA + omega_BEMA* d_chfu$pol_sum + psi_BEMA* d_chfu$salinity) * d_chfu$BEMA 
     + (a_LEMA + omega_LEMA* d_chfu$pol_sum + psi_LEMA* d_chfu$salinity)* d_chfu$LEMA  + (a_MEEL + omega_MEEL* d_chfu$pol_sum + psi_MEEL* d_chfu$salinity)* d_chfu$MEEL 
     + (a_MESU + omega_CHFU* d_chfu$pol_sum + psi_CHFU* d_chfu$salinity)* d_chfu$MESU + (a_PUPA + omega_PUPA* d_chfu$pol_sum + psi_PUPA* d_chfu$salinity)* d_chfu$PUPA 
     + (a_NON_FOCAL + omega_NON_FOCAL* d_chfu$pol_sum + psi_NON_FOCAL* d_chfu$salinity)* d_chfu$NON_FOCAL)
  
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

###########################################################
#models fitting using optim and earlier likelihood functions
###########################################################

## model 1, no competition ----
###############################

log_cover <- jitter(log(d_achmil$Cover +1),factor = 4)

###recall parameters are lambda and sigma- initialize these with estimates from the data:
par1 <- c(mean(log_cover), sd(log_cover))

##repeat optimization until we get convergence (or we try 25 times)
for(k in 1:25){
  result_achmil1 <- optim(par1, compmodel1, method = "L-BFGS-B", lower = c(1, 0.0000000001) , control = list(maxit = 1000, parscale = c(100,0.1), trace = T, REPORT = 100))
  ##update start parameters to final estimate to use in next run in case of nonconvergence
  par1 <- result_achmil1$par
  if(result_achmil1$convergence == 0){
    print(paste("ACHMIL", "model 1 converged on rep", k, sep = " "))
    break
  }
}

## model 2, one alpha ----
#############################

## pars here are lambda, alpha and sigma- use lambda and sigma from model 1 as starting esimtates

par2 <- c(result_achmil1$par[1], 0.0001, result_achmil1$par[2])

##as before:
for(k in 1:25){
  ##now using a specific method that permits constrained optimization so that alpha has to be nonzero- this is an issue in some of the fits, especially in model 3. lower parameter has same order as par2
  result_achmil2<-optim(par2, compmodel2, method = "L-BFGS-B", lower = c(1,-1,0.0000000001), control = list(maxit=1000, parscale = c(100,0.0001,0.1), trace = T, REPORT = 100))
  par2<-result_achmil2$par
  if(result_achmil2$convergence == 0){
    print(paste("ACHMIL", "model 2 converged on rep", k, sep = " "))
    break
  }
}

## model 3a, one alpha, one grasshopper ----
###############################################

par3 <- c(result_achmil2$par[1], 0, 0, result_achmil2$par[2],0 ,0, result_achmil2$par[3])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_achmil3a <- optim(par3, compmodel3a, method = "L-BFGS-B", lower = c(0.1, 0, 0, -0.002052232, 0, 0, 0.0001), 
                           upper=c(100, 0.001, 0.001, 2, 0.001, 0.001, 0.70), control = list(maxit=1000)) #parscale = c(1, 0.01, 0.01, 0.01), trace = T, REPORT = 100))
  par3 <- result_achmil3a$par
  if(result_achmil3a$convergence == 0){
    print(paste("ACHMIL", "model 3 converged on rep", k, sep = " "))
    break
  }
}


## model 3, one alpha, salt and polinators ----
###############################################

par3 <- c(result_achmil2$par[1], rep(0.1, times = 6), result_achmil2$par[2], rep(0.1, times = 5), result_achmil2$par[3])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_achmil3 <- optim(par3, compmodel3, method = "L-BFGS-B", lower = c(-1, rep(-0.5, times = 6), 0, rep(-0.2, times = 5), 0.0000000001), control = list(maxit=1000, parscale = c(10, rep(0.01, times = 6), 0.1, rep(0.01, times = 5), 0.01), trace = T, REPORT = 100))
  par3 <- result_achmil3$par
  if(result_achmil3$convergence == 0){
    print(paste("ACHMIL", "model 3 converged on rep", k, sep = " "))
    break
  }
}

###  MODIFIED UP TO THIS POINT ------------------------
## model 4, several alphas, no salt and polinators ----
#######################################################

par4 <- c(result_chfu3$par[1], rep(result_chfu3$par[4], times = 7), result_chfu3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_chfu4 <- optim(par4, compmodel4, method="L-BFGS-B", lower = c(1, rep(0, times = 7), 0.0000000001), control = list(maxit = 1000, parscale = c(100, rep(0.1, times = 8)), trace = T, REPORT = 100))
  par4 <- result_chfu4$par
  if(result_chfu4$convergence == 0){
    print(paste("CHFU", "model 4 converged on rep", k, sep = " "))
    break
  }
}


#############################
## model 5, several alphas, salt and polinators ----
#############################

par5<-c(result_chfu3$par[1], result_chfu3$par[2], result_chfu3$par[3], result_chfu4$par[2],
        result_chfu4$par[3], result_chfu4$par[4], result_chfu4$par[5], result_chfu4$par[6], result_chfu4$par[7],
        result_chfu4$par[8], result_chfu3$par[5], result_chfu3$par[6], result_chfu3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_chfu5<-optim(par5,compmodel5, method="L-BFGS-B", hessian=TRUE, lower=c(1,-5,-5, rep(0.0001, times=7), -5, -5, 0.0000000001), control=list(maxit=1000, parscale=c(100, 0.1, 10, 0.01, 0.01, 0.1, 0.001, 0.1, 0.1, 0.1, 0.01, 0.01, 0.1), trace= T, REPORT= 100))
  par5<-result_chfu5$par  
  if(result_chfu5$convergence==0){
    print(paste("CHFU", "model 5 converged on rep", k, sep=" "))
    break
    
  }}

inverse <- solve(result_chfu5$hessian)
errors <- sqrt(diag(inverse))  
upper <- result_chfu5$par+1.96*errors
lower <- result_chfu5$par-1.96*errors

det(result_chfu5$hessian)

#############################
## model 6, several alphas, salt and polinators ----
#############################

par6<-c(result_chfu3$par[1], result_chfu3$par[2], result_chfu3$par[3], result_chfu4$par[2],
        result_chfu4$par[3], result_chfu4$par[4], result_chfu4$par[5], result_chfu4$par[6], result_chfu4$par[7],
        result_chfu4$par[8], rep(0.01, times=7), rep(0.01, times=7), result_chfu3$par[7])
        
        
##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_chfu6<-optim(par6,compmodel6, method="L-BFGS-B", lower=c(1,-5,-5, rep(0.0001, times=7), rep(-5, times=7), rep(-5, times=7), 0.0000000001), control=list(maxit=1000, parscale=c(100, 0.1, 10, 0.01, 0.01, 0.1, 0.001, 0.1, 0.1, 0.1, rep(0.01, times=7), rep(0.01, times=7), 0.1), trace= T, REPORT= 100))
  par6<-result_chfu6$par  
  if(result_chfu6$convergence==0){
    print(paste("CHFU", "model 6 converged on rep", k, sep=" "))
    break
    
  }}


###Save results
alpha_matrix[1, 1:7]<-par5[4:10]
matrix[1, 1:7]<-par3[1:7]


saveRDS(result_chfu1, file = "results/chfu_results/result_chfu1.rds")
saveRDS(result_chfu2, file = "results/chfu_results/result_chfu2.rds")
saveRDS(result_chfu3, file = "results/chfu_results/result_chfu3.rds")
saveRDS(result_chfu4, file = "results/chfu_results/result_chfu4.rds")
saveRDS(result_chfu5, file = "results/chfu_results/result_chfu5.rds")
saveRDS(result_chfu6, file = "results/chfu_results/result_chfu6.rds")






#################
################# end of one unit (1 spp)
#################






#BEMA----
#model 1 - no effect of density (no competitive effects)
compmodel1<-function(par){
  
  ## mu (=lambda later) and sigma parameters for the normal distribution (assuming lognormal error- seed data are logged)
  lambda<-par[1]
  sigma<-par[2]
  
  #this the predictive model- here is just fitting a horizontal line through the data:
  pred<-rep(lambda, times=length(log_seeds)) 
  #these are the log likelihoods of the data given the model + parameters
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE) 
  #return the sum of negative log likelihoods - what optim minimizes
  return(sum(-1*llik)) 
}

#model 2 - competition, but no difference between species
compmodel2<-function(par){
  lambda <- par[1] ## same as model 1
  alpha <- par[2]  ## new parameter introduced in model 2
  sigma <- par[3] ## same as model 1
  # predictive model:
  pred <- lambda/(1+alpha*(d_bema$n_nei)) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}

#model 3- common effect of competition including salinity and pollinators. 
compmodel3<-function(par){
  lambda <- par[1] ## same as model 1
  gamma<-par[2]
  theta<-par[3]
  alpha <- par[4]  ## new parameter introduced in model 2
  omega<-par[5]
  psi<-par[6]
  sigma <- par[7] ## same as model 1
  # predictive model:
  pred <- lambda*(1 + theta* d_bema$salinity + gamma* d_bema$pol_sum)/(1+ (alpha + omega* d_bema$pol_sum + psi* d_bema$salinity)* d_bema$n_nei) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}


#model 4 - all species have different competitive effects
compmodel4<-function(par){
  lambda<-par[1] #same as model 2
  a_CHFU<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
  a_BEMA<-par[3]
  a_LEMA<-par[4]
  a_MEEL<-par[5]
  a_MESU<-par[6]
  a_PUPA<-par[7]
  a_NON_FOCAL<-par[8]
  sigma<-par[9] ## same as model 2
  
  ##probably overkill to name all the alphas, but it helps keep things straight
  
  # same form as model 1, but super long given all the species:
  pred<- lambda/ (1+ a_CHFU* d_bema$CHFU + a_BEMA * d_bema$BEMA + a_LEMA* d_bema$LEMA + a_MEEL* d_bema$MEEL + a_MESU* d_bema$MESU + a_PUPA* d_bema$PUPA + a_NON_FOCAL* d_bema$NON_FOCAL)
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


#model 5 - 
compmodel5<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega<-par[11]
  psi<-par[12]
  sigma<-par[13]
  
  pred<- lambda*(1 + theta* d_bema$salinity + gamma* d_bema$pol_sum)/ 
    (1+ (a_CHFU + omega* d_bema$pol_sum + psi* d_bema$salinity)* d_bema$CHFU + (a_BEMA + omega* d_bema$pol_sum + psi* d_bema$salinity) * d_bema$BEMA 
     + (a_LEMA + omega* d_bema$pol_sum + psi* d_bema$salinity)* d_bema$LEMA  + (a_MEEL + omega* d_bema$pol_sum + psi* d_bema$salinity)* d_bema$MEEL 
     + (a_MESU + omega* d_bema$pol_sum + psi* d_bema$salinity)* d_bema$MESU + (a_PUPA + omega* d_bema$pol_sum + psi* d_bema$salinity)* d_bema$PUPA 
     + (a_NON_FOCAL + omega* d_bema$pol_sum + psi* d_bema$salinity)* d_bema$NON_FOCAL)
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


#model6


#model 6 - 
compmodel6<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega_CHFU<-par[11]
  omega_BEMA<-par[12]
  omega_LEMA<-par[13]
  omega_MEEL<-par[14]
  omega_MESU<-par[15]
  omega_PUPA<-par[16]
  omega_NON_FOCAL<-par[17]
  psi_CHFU<-par[18]
  psi_BEMA<-par[19]
  psi_LEMA<-par[20]
  psi_MEEL<-par[21]
  psi_MESU<-par[22]
  psi_PUPA<-par[23]
  psi_NON_FOCAL<-par[24]
  sigma<-par[25]
  
  pred<- lambda*(1 + theta* d_bema$salinity + gamma* d_bema$pol_sum)/ 
    
    (1+ (a_CHFU + omega_CHFU* d_bema$pol_sum + psi_CHFU* d_bema$salinity)* d_bema$CHFU + (a_BEMA + omega_BEMA* d_bema$pol_sum + psi_BEMA* d_bema$salinity) * d_bema$BEMA 
     + (a_LEMA + omega_LEMA* d_bema$pol_sum + psi_LEMA* d_bema$salinity)* d_bema$LEMA  + (a_MEEL + omega_MEEL* d_bema$pol_sum + psi_MEEL* d_bema$salinity)* d_bema$MEEL 
     + (a_MESU + omega_CHFU* d_bema$pol_sum + psi_CHFU* d_bema$salinity)* d_bema$MESU + (a_PUPA + omega_PUPA* d_bema$pol_sum + psi_PUPA* d_bema$salinity)* d_bema$PUPA 
     + (a_NON_FOCAL + omega_NON_FOCAL* d_bema$pol_sum + psi_NON_FOCAL* d_bema$salinity)* d_bema$NON_FOCAL)
  
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


#las funciones, son identicas, verdad? Si es asi, no las cargues cada vez. Multiplica la posibilidad de arrastrat errores si modificas una, y no aporta nada.
log_seeds<-log(d_bema$seed_number)

#model fitting using optim and earlier likelihood functions

#############################
## model 1, no competition ----
#############################

###recall parameters are lambda and sigma- initialize these with estimates from the data:
par1<-c(mean(log_seeds), sd(log_seeds))

##repeat optimization until we get convergence (or we try 25 times)
for(k in 1:25){
  result_bema1<-optim(par1,compmodel1, method="L-BFGS-B", lower=c(1, 0.0000000001), control=list(maxit=1000,parscale=c(10,0.1), trace= T, REPORT= 100))
  ##update start parameters to final estimate to use in next run in case of nonconvergence
  par1<-result_bema1$par
  if(result_bema1$convergence==0){
    print(paste("BEMA", "model 1 converged on rep", k, sep=" "))
    break
    
  }}



#############################
## model 2, one alpha ----
#############################

## pars here are lambda, alpha and sigma- use lambda and sigma from model 1 as starting esimtates

par2<-c(result_bema1$par[1], 0.0001, result_bema1$par[2])

##as before:
for(k in 1:25){
  ##now using a specific method that permits constrained optimization so that alpha has to be nonzero- this is an issue in some of the fits, especially in model 3. lower parameter has same order as par2
  result_bema2<-optim(par2,compmodel2, method="L-BFGS-B", lower=c(1,0, 0.0000000001), control=list(maxit=1000,parscale=c(10,0.0001,0.1), trace= T, REPORT= 100))
  par2<-result_bema2$par
  if(result_bema2$convergence==0){
    print(paste("BEMA", "model 2 converged on rep", k, sep=" "))
    break
    
  }}


#############################
## model 3, one alpha, salt and polinators ----
#############################


par3<-c(result_bema2$par[1], 0.0001, 0.0001, result_bema2$par[2], 0.0001, 0.0001, result_bema2$par[3]) 

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_bema3<-optim(par3,compmodel3, method="L-BFGS-B", lower=c(1, -5, -5, 0, -5, -5, 0.0000000001), control=list(maxit=1000, parscale=c(10,0.0001,0.0001, 0.0001, 0.0001, 0.0001, 0.1), trace= T, REPORT= 100))
  par3<-result_bema3$par
  if(result_bema3$convergence==0){
    print(paste("BEMA", "model 3 converged on rep", k, sep=" "))
    break
    
  }}

#############################
## model 4, several alphas, no salt and polinators ----
#############################

par4<-c(result_bema3$par[1], rep(result_bema3$par[4], times=7),result_bema3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_bema4<-optim(par4,compmodel4, method="L-BFGS-B", hessian=TRUE, lower=c(1, rep(0.00001, times=7), 0.0000000001), control=list(maxit=1000, parscale=c(10, rep(0.0001, times=7), 0.1), trace= T, REPORT= 100))
  par4<-result_bema4$par
  if(result_bema4$convergence==0){
    print(paste("BEMA", "model 4 converged on rep", k, sep=" "))
    break
    
  }}
#voy a modificar la matrix para ver si consigo hacerla invertible

result_bema4$hessian[result_bema4$hessian > 10] <- 10
result_bema4$hessian[result_bema4$hessian == 0] <- 0.001


inverse <- solve(result_bema4$hessian)
errors <- sqrt(diag(inverse))  
upper <- result_bema4$par+1.96*errors
lower <- result_bema4$par-1.96*errors


#############################
## model 5, several alphas, salt and polinators ----
#############################

par5<-c(result_bema3$par[1], result_bema3$par[2], result_bema3$par[3], result_bema4$par[2],
        result_bema4$par[3], result_bema4$par[4], result_bema4$par[5], result_bema4$par[6], result_bema4$par[7],
        result_bema4$par[8], result_bema3$par[5], result_bema3$par[6], result_bema3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_bema5<-optim(par5,compmodel5, method="L-BFGS-B", hessian=T, lower=c(1, -5,-5, rep(0.00001, times=7), -5, -5, 0.0000000001),  control=list(maxit=1000, parscale=c(10, rep(0.0001, times=11), 0.1), trace= T, REPORT= 100))
  par5<-result_bema5$par
  if(result_bema5$convergence==0){
    print(paste("BEMA", "model 5 converged on rep", k, sep=" "))
    break
    
  }}

inverse <- solve(result_bema4$hessian)
errors <- sqrt(diag(inverse))  
upper <- result_bema4$par+1.96*errors
lower <- result_bema4$par-1.96*errors

#############################
## model 6
#############################

par6<-c(result_bema3$par[1], result_bema3$par[2], result_bema3$par[3], result_bema4$par[2],
        result_bema4$par[3], result_bema4$par[4], result_bema4$par[5], result_bema4$par[6], result_bema4$par[7],
        result_bema4$par[8], rep(0.0001, times=7), rep(0.0001, times=7), result_bema3$par[7])


##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_bema6<-optim(par6,compmodel6, method="L-BFGS-B", lower=c(1,-5,-5, rep(0.0001, times=7), rep(-5, times=7), rep(-5, times=7), 0.0000000001), control=list(maxit=1000, parscale=c(100, 0.001, 0.0001, 0.00001, 0.01, 0.00001, 0.0001, 0.0001, 0.00001, 0.001, rep(0.0001, times=7), rep(0.0001, times=7), 0.1), trace= T, REPORT= 100))
  par6<-result_bema6$par  
  if(result_bema6$convergence==0){
    print(paste("BEMA", "model 6 converged on rep", k, sep=" "))
    break
    
  }}



###Save results
#Poner par4 en vez de par5

alpha_matrix[2, 1:7]<-par4[2:8]
matrix[2, 1:7]<-par3[1:7]


saveRDS(result_bema1, file = "results/bema_results/result_bema1.rds")
saveRDS(result_bema2, file = "results/bema_results/result_bema2.rds")
saveRDS(result_bema3, file = "results/bema_results/result_bema3.rds")
saveRDS(result_bema4, file = "results/bema_results/result_bema4.rds")
saveRDS(result_bema5, file = "results/bema_results/result_bema5.rds")
saveRDS(result_bema6, file = "results/bema_results/result_bema6.rds")




















#LEMA----
#model 1 - no effect of density (no competitive effects)
compmodel1<-function(par){
  
  ## mu (=lambda later) and sigma parameters for the normal distribution (assuming lognormal error- seed data are logged)
  lambda<-par[1]
  sigma<-par[2]
  
  #this the predictive model- here is just fitting a horizontal line through the data:
  pred<-rep(lambda, times=length(log_seeds)) 
  #these are the log likelihoods of the data given the model + parameters
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE) 
  #return the sum of negative log likelihoods - what optim minimizes
  return(sum(-1*llik)) 
}

#model 2 - competition, but no difference between species
compmodel2<-function(par){
  lambda <- par[1] ## same as model 1
  alpha <- par[2]  ## new parameter introduced in model 2
  sigma <- par[3] ## same as model 1
  # predictive model:
  pred <- lambda/(1+alpha*(d_lema$n_nei)) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}

#model 3- common effect of competition including salinity and pollinators. 
compmodel3<-function(par){
  lambda <- par[1] ## same as model 1
  gamma<-par[2]
  theta<-par[3]
  alpha <- par[4]  ## new parameter introduced in model 2
  omega<-par[5]
  psi<-par[6]
  sigma <- par[7] ## same as model 1
  # predictive model:
  pred <- lambda*(1 + theta* d_lema$salinity + gamma* d_lema$pol_sum)/(1+ (alpha + omega* d_lema$pol_sum + psi* d_lema$salinity)* d_lema$n_nei) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}


#model 4 - all species have different competitive effects
compmodel4<-function(par){
  lambda<-par[1] #same as model 2
  a_CHFU<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
  a_BEMA<-par[3]
  a_LEMA<-par[4]
  a_MEEL<-par[5]
  a_MESU<-par[6]
  a_PUPA<-par[7]
  a_NON_FOCAL<-par[8]
  sigma<-par[9] ## same as model 2
  
  ##probably overkill to name all the alphas, but it helps keep things straight
  
  # same form as model 1, but super long given all the species:
  pred<- lambda/ (1+ a_CHFU* d_lema$CHFU + a_BEMA * d_lema$BEMA + a_LEMA* d_lema$LEMA + a_MEEL* d_lema$MEEL + a_MESU* d_lema$MESU + a_PUPA* d_lema$PUPA + a_NON_FOCAL* d_lema$NON_FOCAL)
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


#model 5 - 
compmodel5<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega<-par[11]
  psi<-par[12]
  sigma<-par[13]
  
  pred<- lambda*(1 + theta* d_lema$salinity + gamma* d_lema$pol_sum)/ 
    (1+ (a_CHFU + omega* d_lema$pol_sum + psi* d_lema$salinity)* d_lema$CHFU + (a_BEMA + omega* d_lema$pol_sum + psi* d_lema$salinity) * d_lema$BEMA 
     + (a_LEMA + omega* d_lema$pol_sum + psi* d_lema$salinity)* d_lema$LEMA  + (a_MEEL + omega* d_lema$pol_sum + psi* d_lema$salinity)* d_lema$MEEL 
     + (a_MESU + omega* d_lema$pol_sum + psi* d_lema$salinity)* d_lema$MESU + (a_PUPA + omega* d_lema$pol_sum + psi* d_lema$salinity)* d_lema$PUPA 
     + (a_NON_FOCAL + omega* d_lema$pol_sum + psi* d_lema$salinity)* d_lema$NON_FOCAL)
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

#model6


#model 6 - 
compmodel6<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega_CHFU<-par[11]
  omega_BEMA<-par[12]
  omega_LEMA<-par[13]
  omega_MEEL<-par[14]
  omega_MESU<-par[15]
  omega_PUPA<-par[16]
  omega_NON_FOCAL<-par[17]
  psi_CHFU<-par[18]
  psi_BEMA<-par[19]
  psi_LEMA<-par[20]
  psi_MEEL<-par[21]
  psi_MESU<-par[22]
  psi_PUPA<-par[23]
  psi_NON_FOCAL<-par[24]
  sigma<-par[25]
  
  pred<- lambda*(1 + theta* d_lema$salinity + gamma* d_lema$pol_sum)/ 
    
    (1+ (a_CHFU + omega_CHFU* d_lema$pol_sum + psi_CHFU* d_lema$salinity)* d_lema$CHFU + (a_BEMA + omega_BEMA* d_lema$pol_sum + psi_BEMA* d_lema$salinity) * d_lema$BEMA 
     + (a_LEMA + omega_LEMA* d_lema$pol_sum + psi_LEMA* d_lema$salinity)* d_lema$LEMA  + (a_MEEL + omega_MEEL* d_lema$pol_sum + psi_MEEL* d_lema$salinity)* d_lema$MEEL 
     + (a_MESU + omega_CHFU* d_lema$pol_sum + psi_CHFU* d_lema$salinity)* d_lema$MESU + (a_PUPA + omega_PUPA* d_lema$pol_sum + psi_PUPA* d_lema$salinity)* d_lema$PUPA 
     + (a_NON_FOCAL + omega_NON_FOCAL* d_lema$pol_sum + psi_NON_FOCAL* d_lema$salinity)* d_lema$NON_FOCAL)
  
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

log_seeds<-log(d_lema$seed_number)

#model fitting using optim and earlier likelihood functions

#############################
## model 1, no competition ----
#############################

###recall parameters are lambda and sigma- initialize these with estimates from the data:

par1<-c(mean(log_seeds), sd(log_seeds))

##repeat optimization until we get convergence (or we try 25 times)
for(k in 1:25){
  result_lema1<-optim(par1,compmodel1, method="L-BFGS-B", lower=c(1, 0.0000000001), control=list(maxit=1000,parscale=c(10,0.1), trace= T, REPORT= 100))
  ##update start parameters to final estimate to use in next run in case of nonconvergence
  par1<-result_lema1$par
  if(result_lema1$convergence==0){
    print(paste("LEMA", "model 1 converged on rep", k, sep=" "))
    break
    
  }}

#############################
## model 2, one alpha ----
#############################

## pars here are lambda, alpha and sigma- use lambda and sigma from model 1 as starting esimtates

par2<-c(result_lema1$par[1],0.001,result_lema1$par[2])

##as before:
for(k in 1:25){
  ##now using a specific method that permits constrained optimization so that alpha has to be nonzero- this is an issue in some of the fits, especially in model 3. lower parameter has same order as par2
  result_lema2<-optim(par2,compmodel2, method="L-BFGS-B", lower=c(1,0,0.0000000001), control=list(maxit=1000,parscale=c(100,0.001, 0.1), trace= T, REPORT= 100))
  par2<-result_lema2$par
  if(result_lema2$convergence==0){
    print(paste("LEMA", "model 2 converged on rep", k, sep=" "))
    break
    
  }}


#############################
## model 3, one alpha, salt and polinators ----
#############################

par3<-c(result_lema2$par[1], 0.0001,0.1, result_lema2$par[2], 0.0001,0.1,result_lema2$par[3])

##as before

for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_lema3<-optim(par3,compmodel3, method="L-BFGS-B", lower=c(1, -5,-5, 0, -5, -5, 0.0000000001), control=list(maxit=1000,parscale=c(100, 0.0001, 0.1, 0.1, 0.0001, 0.1, 0.1), trace= T, REPORT= 100))
  par3<-result_lema3$par
  if(result_lema3$convergence==0){
    print(paste("LEMA", "model 3 converged on rep", k, sep=" "))
    break
    
  }}

#############################
## model 4, several alphas, no salt and polinators ----
#############################

par4<-c(result_lema2$par[1], rep(result_lema3$par[4], times=7), result_lema3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_lema4<-optim(par4,compmodel4, method="L-BFGS-B", lower=c(1, rep(0.00001, times=7), 0.0000000001), control=list(maxit=1000,parscale=c(100, rep(0.01, times=7), 0.1), trace= T, REPORT= 100))
  par4<-result_lema4$par
  if(result_lema4$convergence==0){
    print(paste("LEMA", "model 4 converged on rep", k, sep=" "))
    break
    
  }}


#############################
## model 5, several alphas, salt and polinators ----
#############################

par5<-c(result_lema2$par[1], result_lema3$par[2], result_lema3$par[3], result_lema4$par[2],
        result_lema4$par[3], result_lema4$par[4], result_lema4$par[5], result_lema4$par[6], result_lema4$par[7],
        result_lema4$par[8], result_lema3$par[5], result_lema3$par[6], result_lema3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_lema5<-optim(par5,compmodel5, method="L-BFGS-B", hessian=TRUE, lower=c(1, -5,-5, rep(0.00001, times=7), -5, -5, 0.0000000001), control=list(maxit=1000,parscale=c(100, 0.0001, 1, 0.01, 0.00001, 0.01, 0.00001, 0.01, 0.01, 0.01, 0.001, 0.01, 0.1), trace= T, REPORT= 100))
  par5<-result_lema5$par
  if(result_lema5$convergence==0){
    print(paste("LEMA", "model 5 converged on rep", k, sep=" "))
    break
    
  }}

inverse <- solve(result_lema5$hessian)
errors <- sqrt(diag(inverse))  
upper <- result_lema5$par+1.96*errors
lower <- result_lema5$par-1.96*errors


#############################
## model 6
#############################

par6<-c(result_lema3$par[1], result_lema3$par[2], result_lema3$par[3], result_lema4$par[2],
        result_lema4$par[3], result_lema4$par[4], result_lema4$par[5], result_lema4$par[6], result_lema4$par[7],
        result_lema4$par[8], rep(0.0001, times=7), rep(0.0001, times=7), result_lema3$par[7])


##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_lema6<-optim(par6,compmodel6, method="L-BFGS-B", lower=c(1,-5,-5, rep(0.0001, times=7), rep(-5, times=7), rep(-5, times=7), 0.0000000001), control=list(maxit=1000, parscale=c(100, 0.01, 0.01, 0.01, 0.00001, 0.01, 0.01, 0.01, 0.00001, 0.00001, rep(0.00001, times=7), rep(0.00001, times=7), 0.1), trace= T, REPORT= 100))
  par6<-result_lema6$par  
  if(result_lema6$convergence==0){
    print(paste("LEMA", "model 6 converged on rep", k, sep=" "))
    break
    
  }}


###Save results
alpha_matrix[3, 1:7]<-par5[4:10]
matrix[3, 1:7]<-par3[1:7]



saveRDS(result_lema1, file = "results/lema_results/result_lema1.rds")
saveRDS(result_lema2, file = "results/lema_results/result_lema2.rds")
saveRDS(result_lema3, file = "results/lema_results/result_lema3.rds")
saveRDS(result_lema4, file = "results/lema_results/result_lema4.rds")
saveRDS(result_lema5, file = "results/lema_results/result_lema5.rds")
saveRDS(result_lema6, file = "results/lema_results/result_lema6.rds")






#MEEL----
#model 1 - no effect of density (no competitive effects)
compmodel1<-function(par){
  
  ## mu (=lambda later) and sigma parameters for the normal distribution (assuming lognormal error- seed data are logged)
  lambda<-par[1]
  sigma<-par[2]
  
  #this the predictive model- here is just fitting a horizontal line through the data:
  pred<-rep(lambda, times=length(log_seeds)) 
  #these are the log likelihoods of the data given the model + parameters
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE) 
  #return the sum of negative log likelihoods - what optim minimizes
  return(sum(-1*llik)) 
}

#model 2 - competition, but no difference between species
compmodel2<-function(par){
  lambda <- par[1] ## same as model 1
  alpha <- par[2]  ## new parameter introduced in model 2
  sigma <- par[3] ## same as model 1
  # predictive model:
  pred <- lambda/(1+alpha*(d_meel$n_nei)) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}

#model 3- common effect of competition including salinity and pollinators. 
compmodel3<-function(par){
  lambda <- par[1] ## same as model 1
  gamma<-par[2]
  theta<-par[3]
  alpha <- par[4]  ## new parameter introduced in model 2
  omega<-par[5]
  psi<-par[6]
  sigma <- par[7] ## same as model 1
  # predictive model:
  pred <- lambda*(1 + theta* d_meel$salinity + gamma* d_meel$pol_sum)/(1+ (alpha + omega* d_meel$pol_sum + psi* d_meel$salinity)* d_meel$n_nei) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}


#model 4 - all species have different competitive effects
compmodel4<-function(par){
  lambda<-par[1] #same as model 2
  a_CHFU<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
  a_BEMA<-par[3]
  a_LEMA<-par[4]
  a_MEEL<-par[5]
  a_MESU<-par[6]
  a_PUPA<-par[7]
  a_NON_FOCAL<-par[8]
  sigma<-par[9] ## same as model 2
  
  ##probably overkill to name all the alphas, but it helps keep things straight
  
  # same form as model 1, but super long given all the species:
  pred<- lambda/ (1+ a_CHFU* d_meel$CHFU + a_BEMA * d_meel$BEMA + a_LEMA* d_meel$LEMA + a_MEEL* d_meel$MEEL + a_MESU* d_meel$MESU + a_PUPA* d_meel$PUPA + a_NON_FOCAL* d_meel$NON_FOCAL)
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


#model 5 - 
compmodel5<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega<-par[11]
  psi<-par[12]
  sigma<-par[13]
  
  pred<- lambda*(1 + theta* d_meel$salinity + gamma* d_meel$pol_sum)/ 
    (1+ (a_CHFU + omega* d_meel$pol_sum + psi* d_meel$salinity)* d_meel$CHFU + (a_BEMA + omega* d_meel$pol_sum + psi* d_meel$salinity) * d_meel$BEMA 
     + (a_LEMA + omega* d_meel$pol_sum + psi* d_meel$salinity)* d_meel$LEMA  + (a_MEEL + omega* d_meel$pol_sum + psi* d_meel$salinity)* d_meel$MEEL 
     + (a_MESU + omega* d_meel$pol_sum + psi* d_meel$salinity)* d_meel$MESU + (a_PUPA + omega* d_meel$pol_sum + psi* d_meel$salinity)* d_meel$PUPA 
     + (a_NON_FOCAL + omega* d_meel$pol_sum + psi* d_meel$salinity)* d_meel$NON_FOCAL)
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

#model6


#model 6 - 
compmodel6<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega_CHFU<-par[11]
  omega_BEMA<-par[12]
  omega_LEMA<-par[13]
  omega_MEEL<-par[14]
  omega_MESU<-par[15]
  omega_PUPA<-par[16]
  omega_NON_FOCAL<-par[17]
  psi_CHFU<-par[18]
  psi_BEMA<-par[19]
  psi_LEMA<-par[20]
  psi_MEEL<-par[21]
  psi_MESU<-par[22]
  psi_PUPA<-par[23]
  psi_NON_FOCAL<-par[24]
  sigma<-par[25]
  
  pred<- lambda*(1 + theta* d_meel$salinity + gamma* d_meel$pol_sum)/ 
    
    (1+ (a_CHFU + omega_CHFU* d_meel$pol_sum + psi_CHFU* d_meel$salinity)* d_meel$CHFU + (a_BEMA + omega_BEMA* d_meel$pol_sum + psi_BEMA* d_meel$salinity) * d_meel$BEMA 
     + (a_LEMA + omega_LEMA* d_meel$pol_sum + psi_LEMA* d_meel$salinity)* d_meel$LEMA  + (a_MEEL + omega_MEEL* d_meel$pol_sum + psi_MEEL* d_meel$salinity)* d_meel$MEEL 
     + (a_MESU + omega_CHFU* d_meel$pol_sum + psi_CHFU* d_meel$salinity)* d_meel$MESU + (a_PUPA + omega_PUPA* d_meel$pol_sum + psi_PUPA* d_meel$salinity)* d_meel$PUPA 
     + (a_NON_FOCAL + omega_NON_FOCAL* d_meel$pol_sum + psi_NON_FOCAL* d_meel$salinity)* d_meel$NON_FOCAL)
  
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

log_seeds<-log(d_meel$seed_number)

#model fitting using optim and earlier likelihood functions

#############################
## model 1, no competition ----
#############################

###recall parameters are lambda and sigma- initialize these with estimates from the data:
par1<-c(mean(log_seeds), sd(log_seeds))

##repeat optimization until we get convergence (or we try 25 times)
for(k in 1:25){
  result_meel1<-optim(par1,compmodel1, method="L-BFGS-B", lower=c(1, 0.0000000001), control=list(maxit=1000,parscale=c(1,0.1), trace= T, REPORT= 100))
  ##update start parameters to final estimate to use in next run in case of nonconvergence
  par1<-result_meel1$par
  if(result_meel1$convergence==0){
    print(paste("MEEL", "model 1 converged on rep", k, sep=" "))
    break
  }}


#############################
## model 2, one alpha ----
#############################

## pars here are lambda, alpha and sigma- use lambda and sigma from model 1 as starting esimtates

par2<-c(result_meel1$par[1],0.0001, result_meel1$par[2])

##as before:
for(k in 1:25){
  ##now using a specific method that permits constrained optimization so that alpha has to be nonzero- this is an issue in some of the fits, especially in model 3. lower parameter has same order as par2
  result_meel2<-optim(par2,compmodel2, method="L-BFGS-B", lower=c(1,0,0.0000000001), control=list(maxit=1000,parscale=c(10,0.0001, 0.1), trace= T, REPORT= 100))
  par2<-result_meel2$par
  if(result_meel2$convergence==0){
    print(paste("MEEL", "model 2 converged on rep", k, sep=" "))
    break
  }}

#############################
## model 3, one alpha, salt and polinators ----
#############################

par3<-c(result_meel2$par[1], 0,0, result_meel2$par[2], 0,0,result_meel2$par[3])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_meel3<-optim(par3,compmodel3, method="L-BFGS-B", lower=c(1, -5,-5, 0, -5, -5, 0.0000000001), control=list(maxit=1000,parscale=c(10, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.1), trace= T, REPORT= 100))
  par3<-result_meel3$par
  if(result_meel3$convergence==0){
    print(paste("MEEL", "model 3 converged on rep", k, sep=" "))
    break
  }}


#############################
## model 4, several alphas, no salt and polinators ----
#############################

par4<-c(result_meel3$par[1], rep(result_meel3$par[4], times=7), result_meel3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_meel4<-optim(par4,compmodel4, method="L-BFGS-B", hessian=TRUE, lower=c(1, rep(0.00001, times=7), 0.0000000001), control=list(maxit=1000,parscale=c(10, rep(0.001, times=7), 0.1), trace= T, REPORT= 100))
  par4<-result_meel4$par
  if(result_meel4$convergence==0){
    print(paste("MEEL", "model 4 converged on rep", k, sep=" "))
    break
  }}

inverse <- solve(result_meel4$hessian)
errors <- sqrt(diag(inverse))  
upper <- result_meel4$par+1.96*errors
lower <- result_meel4$par-1.96*errors


#############################
## model 5, several alphas, salt and polinators ----
#############################

par5<-c(result_meel3$par[1], result_meel3$par[2], result_meel3$par[3], result_meel4$par[2],
        result_meel4$par[3], result_meel4$par[4], result_meel4$par[5], result_meel4$par[6], result_meel4$par[7],
        result_meel4$par[8], result_meel3$par[5], result_meel3$par[6], result_meel3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  ##parece dar lo mismo auqneu ponga un parscale especifico para cada valor
  result_meel5<-optim(par5,compmodel5, method="L-BFGS-B", lower=c(1, -5,-5, rep(0.00001, times=7), -5, -5, 0.0000000001), control=list(maxit=1000,parscale=c(10, 0.01, 0.001, 0.001, 0.00001, 0.00001, 0.00001, 0.1, 0.01, 0.001, 0.001, 0.1 , 0.1), trace= T, REPORT= 100))
  if(result_meel5$convergence==0){
    print(paste("MEEL", "model 5 converged on rep", k, sep=" "))
    break
  }}


#############################
## model 6
#############################

par6<-c(result_meel3$par[1], result_meel3$par[2], result_meel3$par[3], result_meel4$par[2],
        result_meel4$par[3], result_meel4$par[4], result_meel4$par[5], result_meel4$par[6], result_meel4$par[7],
        result_meel4$par[8], rep(0.0001, times=7), rep(0.0001, times=7), result_meel3$par[7])


##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_meel6<-optim(par6,compmodel6, method="L-BFGS-B", lower=c(1,-5,-5, rep(0.0001, times=7), rep(-5, times=7), rep(-5, times=7), 0.0000000001), control=list(maxit=1000, parscale=c(10, 0.01, 0.001, 0.001, 0.00001, 0.00001, 0.00001, 0.1, 0.01, 0.001, rep(0.0001, times=7), rep(0.0001, times=7), 0.1), trace= T, REPORT= 100))
  par6<-result_meel6$par  
  if(result_meel6$convergence==0){
    print(paste("MEEL", "model 6 converged on rep", k, sep=" "))
    break
    
  }}

###Save results
alpha_matrix[4, 1:7]<-par4[2:8]
matrix[4, 1:7]<-par3[1:7]



saveRDS(result_meel1, file = "results/meel_results/result_meel1.rds")
saveRDS(result_meel2, file = "results/meel_results/result_meel2.rds")
saveRDS(result_meel3, file = "results/meel_results/result_meel3.rds")
saveRDS(result_meel4, file = "results/meel_results/result_meel4.rds")
saveRDS(result_meel5, file = "results/meel_results/result_meel5.rds")
saveRDS(result_meel6, file = "results/meel_results/result_meel6.rds")





#MESU----
#model 1 - no effect of density (no competitive effects)
compmodel1<-function(par){
  
  ## mu (=lambda later) and sigma parameters for the normal distribution (assuming lognormal error- seed data are logged)
  lambda<-par[1]
  sigma<-par[2]
  
  #this the predictive model- here is just fitting a horizontal line through the data:
  pred<-rep(lambda, times=length(log_seeds)) 
  #these are the log likelihoods of the data given the model + parameters
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE) 
  #return the sum of negative log likelihoods - what optim minimizes
  return(sum(-1*llik)) 
}

#model 2 - competition, but no difference between species
compmodel2<-function(par){
  lambda <- par[1] ## same as model 1
  alpha <- par[2]  ## new parameter introduced in model 2
  sigma <- par[3] ## same as model 1
  # predictive model:
  pred <- lambda/(1+alpha*(d_pupa$n_nei)) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}

#model 3- common effect of competition including salinity and pollinators. 
compmodel3<-function(par){
  lambda <- par[1] ## same as model 1
  gamma<-par[2]
  theta<-par[3]
  alpha <- par[4]  ## new parameter introduced in model 2
  omega<-par[5]
  psi<-par[6]
  sigma <- par[7] ## same as model 1
  # predictive model:
  pred <- lambda*(1 + theta* d_mesu$salinity + gamma* d_mesu$pol_sum)/(1+ (alpha + omega* d_mesu$pol_sum + psi* d_mesu$salinity)* d_mesu$n_nei) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}


#model 4 - all species have different competitive effects
compmodel4<-function(par){
  lambda<-par[1] #same as model 2
  a_CHFU<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
  a_BEMA<-par[3]
  a_LEMA<-par[4]
  a_MEEL<-par[5]
  a_MESU<-par[6]
  a_PUPA<-par[7]
  a_NON_FOCAL<-par[8]
  sigma<-par[9] ## same as model 2
  
  ##probably overkill to name all the alphas, but it helps keep things straight
  
  # same form as model 1, but super long given all the species:
  pred<- lambda/ (1+ a_CHFU* d_mesu$CHFU + a_BEMA * d_mesu$BEMA + a_LEMA* d_mesu$LEMA + a_MEEL* d_mesu$MEEL + a_MESU* d_mesu$MESU + a_PUPA* d_mesu$PUPA + a_NON_FOCAL* d_mesu$NON_FOCAL)
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


#model 5 - 
compmodel5<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega<-par[11]
  psi<-par[12]
  sigma<-par[13]
  
  pred<- lambda*(1 + theta* d_mesu$salinity + gamma* d_mesu$pol_sum)/ 
    (1+ (a_CHFU + omega* d_mesu$pol_sum + psi* d_mesu$salinity)* d_mesu$CHFU + (a_BEMA + omega* d_mesu$pol_sum + psi* d_mesu$salinity) * d_mesu$BEMA 
     + (a_LEMA + omega* d_mesu$pol_sum + psi* d_mesu$salinity)* d_mesu$LEMA  + (a_MEEL + omega* d_mesu$pol_sum + psi* d_mesu$salinity)* d_mesu$MEEL 
     + (a_MESU + omega* d_mesu$pol_sum + psi* d_mesu$salinity)* d_mesu$MESU + (a_PUPA + omega* d_mesu$pol_sum + psi* d_mesu$salinity)* d_mesu$PUPA 
     + (a_NON_FOCAL + omega* d_mesu$pol_sum + psi* d_mesu$salinity)* d_mesu$NON_FOCAL)
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


#model 6 - 
compmodel6<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega_CHFU<-par[11]
  omega_BEMA<-par[12]
  omega_LEMA<-par[13]
  omega_MEEL<-par[14]
  omega_MESU<-par[15]
  omega_PUPA<-par[16]
  omega_NON_FOCAL<-par[17]
  psi_CHFU<-par[18]
  psi_BEMA<-par[19]
  psi_LEMA<-par[20]
  psi_MEEL<-par[21]
  psi_MESU<-par[22]
  psi_PUPA<-par[23]
  psi_NON_FOCAL<-par[24]
  sigma<-par[25]
  
  pred<- lambda*(1 + theta* d_mesu$salinity + gamma* d_mesu$pol_sum)/ 
    
    (1+ (a_CHFU + omega_CHFU* d_mesu$pol_sum + psi_CHFU* d_mesu$salinity)* d_mesu$CHFU + (a_BEMA + omega_BEMA* d_mesu$pol_sum + psi_BEMA* d_mesu$salinity) * d_mesu$BEMA 
     + (a_LEMA + omega_LEMA* d_mesu$pol_sum + psi_LEMA* d_mesu$salinity)* d_mesu$LEMA  + (a_MEEL + omega_MEEL* d_mesu$pol_sum + psi_MEEL* d_mesu$salinity)* d_mesu$MEEL 
     + (a_MESU + omega_CHFU* d_mesu$pol_sum + psi_CHFU* d_mesu$salinity)* d_mesu$MESU + (a_PUPA + omega_PUPA* d_mesu$pol_sum + psi_PUPA* d_mesu$salinity)* d_mesu$PUPA 
     + (a_NON_FOCAL + omega_NON_FOCAL* d_mesu$pol_sum + psi_NON_FOCAL* d_mesu$salinity)* d_mesu$NON_FOCAL)
  
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


log_seeds<-log(d_mesu$seed_number)

#model fitting using optim and earlier likelihood functions

#############################
## model 1, no competition ----
#############################

###recall parameters are lambda and sigma- initialize these with estimates from the data:
par1<-c(mean(log_seeds), sd(log_seeds))

##repeat optimization until we get convergence (or we try 25 times)
for(k in 1:25){
  result_mesu1<-optim(par1,compmodel1, method="L-BFGS-B", lower=c(1, 0.0000000001), control=list(maxit=1000,parscale=c(10,0.1), trace= T, REPORT= 100))
  ##update start parameters to final estimate to use in next run in case of nonconvergence
  par1<-result_mesu1$par
  if(result_mesu1$convergence==0){
    print(paste("MESU", "model 1 converged on rep", k, sep=" "))
    break
  }}



#############################
## model 2, one alpha ----
#############################

## pars here are lambda, alpha and sigma- use lambda and sigma from model 1 as starting esimtates

par2<-c(result_mesu1$par[1],0.001,result_mesu1$par[2])

##as before:
for(k in 1:25){
  ##now using a specific method that permits constrained optimization so that alpha has to be nonzero- this is an issue in some of the fits, especially in model 3. lower parameter has same order as par2
  result_mesu2<-optim(par2,compmodel2, method="L-BFGS-B", lower=c(1,0,0.0000000001), control=list(maxit=1000,parscale=c(10,0.0001, 0.1), trace= T, REPORT= 100))
  par2<-result_mesu2$par  
  if(result_mesu2$convergence==0){
    print(paste("MESU", "model 2 converged on rep", k, sep=" "))
    break
  }}


#############################
## model 3, one alpha, salt and polinators ----
#############################

par3<-c(result_mesu2$par[1], 0.01,0.01, result_mesu2$par[2], 0.01,0.01,result_mesu2$par[3])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_mesu3<-optim(par3,compmodel3, method="L-BFGS-B", lower=c(1, -5,-5, 0, -5, -5, 0.0000000001), control=list(maxit=1000,parscale=c(100, 0.001, 0.001, 0.0001, 0.001, 0.001, 0.1), trace= T, REPORT= 100))
  par3<-result_mesu3$par
  if(result_mesu3$convergence==0){
    print(paste("MESU", "model 2 converged on rep", k, sep=" "))
    break
  }}


#############################
## model 4, several alphas, no salt and polinators ----
#############################

par4<-c(result_mesu3$par[1], rep(result_mesu3$par[4], times=7), result_mesu3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_mesu4<-optim(par4,compmodel4, method="L-BFGS-B", hessian=TRUE, lower=c(1, 0.00001,0.05,0.00001,0.00001,0.00001,0.05,0.00001, 0.0000000001), control=list(maxit=1000,parscale=c(100, rep(0.001, times=7), 0.1), trace= T, REPORT= 100))
  par4<-result_mesu4$par
  if(result_mesu4$convergence==0){
    print(paste("MESU", "model 4 converged on rep", k, sep=" "))
    break
  }}

inverse <- solve(result_mesu4$hessian)
errors <- sqrt(diag(inverse))  
upper <- result_mesu4$par+1.96*errors
lower <- result_mesu4$par-1.96*errors


#############################
## model 5, several alphas, salt and polinators ----
#############################

par5<-c(result_mesu3$par[1], result_mesu3$par[2], result_mesu3$par[3], result_mesu4$par[2],
        result_mesu4$par[3],result_mesu4$par[4], result_mesu4$par[5], result_mesu4$par[6], result_mesu4$par[7],
        result_mesu4$par[8], result_mesu3$par[5], result_mesu3$par[6], result_mesu3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_mesu5<-optim(par5,compmodel5, method="L-BFGS-B", lower=c(1, -5,-5, rep(0.00001, times=7), -5, -5, 0.0000000001), control=list(maxit=1000,parscale=c(100, 0.0001, 0.0001, 0.01, 0.00001, 0.1, 0.01, 0.01, 0.00001, 0.01, 0.001, 0.0001 , 0.1), trace= T, REPORT= 100))
  par5<-result_mesu5$par
  if(result_mesu5$convergence==0){
    print(paste("MESU", "model 5 converged on rep", k, sep=" "))
    break
  }}


#############################
## model 6
#############################

par6<-c(result_mesu3$par[1], result_mesu3$par[2], result_mesu3$par[3], result_mesu4$par[2],
        result_mesu4$par[3], result_mesu4$par[4], result_mesu4$par[5], result_mesu4$par[6], result_mesu4$par[7],
        result_mesu4$par[8], rep(0.0001, times=7), rep(0.0001, times=7), result_mesu3$par[7])


##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_mesu6<-optim(par6,compmodel6, method="L-BFGS-B", lower=c(1,-5,-5, rep(0.0001, times=7), rep(-5, times=7), rep(-5, times=7), 0.0000000001), control=list(maxit=1000, parscale=c(100, 0.001, 0.0001, 0.01, 0.00001, 0.1, 0.01, 0.01, 0.00001, 0.01, rep(0.0001, times=7), rep(0.0001, times=7), 0.1), trace= T, REPORT= 100))
  par6<-result_mesu6$par  
  if(result_mesu6$convergence==0){
    print(paste("MESU", "model 6 converged on rep", k, sep=" "))
    break
    
  }}

###Save results
alpha_matrix[5, 1:7]<-par4[2:8]
matrix[5, 1:7]<-par3[1:7]


saveRDS(result_mesu1, file = "results/mesu_results/result_mesu1.rds")
saveRDS(result_mesu2, file = "results/mesu_results/result_mesu2.rds")
saveRDS(result_mesu3, file = "results/mesu_results/result_mesu3.rds")
saveRDS(result_mesu4, file = "results/mesu_results/result_mesu4.rds")
saveRDS(result_mesu5, file = "results/mesu_results/result_mesu5.rds")
saveRDS(result_mesu6, file = "results/mesu_results/result_mesu6.rds")



#PUPA----
#model 1 - no effect of density (no competitive effects)
compmodel1<-function(par){
  
  ## mu (=lambda later) and sigma parameters for the normal distribution (assuming lognormal error- seed data are logged)
  lambda<-par[1]
  sigma<-par[2]
  
  #this the predictive model- here is just fitting a horizontal line through the data:
  pred<-rep(lambda, times=length(log_seeds)) 
  #these are the log likelihoods of the data given the model + parameters
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE) 
  #return the sum of negative log likelihoods - what optim minimizes
  return(sum(-1*llik)) 
}

#model 2 - competition, but no difference between species
compmodel2<-function(par){
  lambda <- par[1] ## same as model 1
  alpha <- par[2]  ## new parameter introduced in model 2
  sigma <- par[3] ## same as model 1
  # predictive model:
  pred <- lambda/(1+alpha*(d_pupa$n_nei)) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}

#model 3- common effect of competition including salinity and pollinators. 
compmodel3<-function(par){
  lambda <- par[1] ## same as model 1
  gamma<-par[2]
  theta<-par[3]
  alpha <- par[4]  ## new parameter introduced in model 2
  omega<-par[5]
  psi<-par[6]
  sigma <- par[7] ## same as model 1
  # predictive model:
  pred <- lambda*(1 + theta* d_pupa$salinity + gamma* d_pupa$pol_sum)/(1+ (alpha + omega* d_pupa$pol_sum + psi* d_pupa$salinity)* d_pupa$n_nei) 
  # log likelihoods of data given the model + parameters:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  # return sum of negative log likelihoods:
  return(sum(-1*llik)) 
}


#model 4 - all species have different competitive effects
compmodel4<-function(par){
  lambda<-par[1] #same as model 2
  a_CHFU<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
  a_BEMA<-par[3]
  a_LEMA<-par[4]
  a_MEEL<-par[5]
  a_MESU<-par[6]
  a_PUPA<-par[7]
  a_NON_FOCAL<-par[8]
  sigma<-par[9] ## same as model 2
  
  ##probably overkill to name all the alphas, but it helps keep things straight
  
  # same form as model 1, but super long given all the species:
  pred<- lambda/ (1+ a_CHFU* d_pupa$CHFU + a_BEMA * d_pupa$BEMA + a_LEMA* d_pupa$LEMA + a_MEEL* d_pupa$MEEL + a_MESU* d_pupa$MESU + a_PUPA* d_pupa$PUPA + a_NON_FOCAL* d_pupa$NON_FOCAL)
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


#model 5 - 
compmodel5<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega<-par[11]
  psi<-par[12]
  sigma<-par[13]
  
  pred<- lambda*(1 + theta* d_pupa$salinity + gamma* d_pupa$pol_sum)/ 
    (1+ (a_CHFU + omega* d_pupa$pol_sum + psi* d_pupa$salinity)* d_pupa$CHFU + (a_BEMA + omega* d_pupa$pol_sum + psi* d_pupa$salinity) * d_pupa$BEMA 
     + (a_LEMA + omega* d_pupa$pol_sum + psi* d_pupa$salinity)* d_pupa$LEMA  + (a_MEEL + omega* d_pupa$pol_sum + psi* d_pupa$salinity)* d_pupa$MEEL 
     + (a_MESU + omega* d_pupa$pol_sum + psi* d_pupa$salinity)* d_pupa$MESU + (a_PUPA + omega* d_pupa$pol_sum + psi* d_pupa$salinity)* d_pupa$PUPA 
     + (a_NON_FOCAL + omega* d_pupa$pol_sum + psi* d_pupa$salinity)* d_pupa$NON_FOCAL)
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}


#model 6 - 
compmodel6<-function(par){
  lambda<-par[1] 
  gamma<-par[2]
  theta<-par[3]
  a_CHFU<-par[4]
  a_BEMA<-par[5]
  a_LEMA<-par[6]
  a_MEEL<-par[7]
  a_MESU<-par[8]
  a_PUPA<-par[9]
  a_NON_FOCAL<-par[10]
  omega_CHFU<-par[11]
  omega_BEMA<-par[12]
  omega_LEMA<-par[13]
  omega_MEEL<-par[14]
  omega_MESU<-par[15]
  omega_PUPA<-par[16]
  omega_NON_FOCAL<-par[17]
  psi_CHFU<-par[18]
  psi_BEMA<-par[19]
  psi_LEMA<-par[20]
  psi_MEEL<-par[21]
  psi_MESU<-par[22]
  psi_PUPA<-par[23]
  psi_NON_FOCAL<-par[24]
  sigma<-par[25]
  
  pred<- lambda*(1 + theta* d_pupa$salinity + gamma* d_pupa$pol_sum)/ 
    
    (1+ (a_CHFU + omega_CHFU* d_pupa$pol_sum + psi_CHFU* d_pupa$salinity)* d_pupa$CHFU + (a_BEMA + omega_BEMA* d_pupa$pol_sum + psi_BEMA* d_pupa$salinity) * d_pupa$BEMA 
     + (a_LEMA + omega_LEMA* d_pupa$pol_sum + psi_LEMA* d_pupa$salinity)* d_pupa$LEMA  + (a_MEEL + omega_MEEL* d_pupa$pol_sum + psi_MEEL* d_pupa$salinity)* d_pupa$MEEL 
     + (a_MESU + omega_CHFU* d_pupa$pol_sum + psi_CHFU* d_pupa$salinity)* d_pupa$MESU + (a_PUPA + omega_PUPA* d_pupa$pol_sum + psi_PUPA* d_pupa$salinity)* d_pupa$PUPA 
     + (a_NON_FOCAL + omega_NON_FOCAL* d_pupa$pol_sum + psi_NON_FOCAL* d_pupa$salinity)* d_pupa$NON_FOCAL)
  
  
  
  # likelihood as before:
  llik<-dnorm(log_seeds,log(pred), sd=sigma, log=TRUE)
  
  # return sum of negative log likelihoods
  return(sum(-1*llik)) #sum of negative log likelihoods
}

log_seeds<-log(d_pupa$seed_number)

#model fitting using optim and earlier likelihood functions

#############################
## model 1, no competition ----
#############################

###recall parameters are lambda and sigma- initialize these with estimates from the data:
par1<-c(mean(log_seeds), sd(log_seeds))

##repeat optimization until we get convergence (or we try 25 times)
for(k in 1:25){
  result_pupa1<-optim(par1,compmodel1, method="L-BFGS-B", lower=c(1, 0.0000000001), control=list(maxit=1000,parscale=c(10,0.1), trace= T, REPORT= 100))
  ##update start parameters to final estimate to use in next run in case of nonconvergence
  par1<-result_pupa1$par
  if(result_pupa1$convergence==0){
    print(paste("PUPA", "model 1 converged on rep", k, sep=" "))
    break
  }}



#############################
## model 2, one alpha ----
#############################

## pars here are lambda, alpha and sigma- use lambda and sigma from model 1 as starting esimtates

par2<-c(result_pupa1$par[1],0.01,result_pupa1$par[2])

##as before:
for(k in 1:25){
  ##now using a specific method that permits constrained optimization so that alpha has to be nonzero- this is an issue in some of the fits, especially in model 3. lower parameter has same order as par2
  result_pupa2<-optim(par2,compmodel2, method="L-BFGS-B", lower=c(1,0,0.0000000001), control=list(maxit=1000,parscale=c(1000,0.01, 0.1), trace= T, REPORT= 100))
  par2<-result_pupa2$par
  if(result_pupa2$convergence==0){
    print(paste("PUPA", "model 2 converged on rep", k, sep=" "))
    break
  }}


#############################
## model 3, one alpha, salt and polinators ----
#############################

par3<-c(result_pupa2$par[1], 0.01, 0.001, result_pupa2$par[2], 0.01,0.001,result_pupa2$par[3])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_pupa3<-optim(par3,compmodel3, method="L-BFGS-B", lower=c(1, -5,-5, 0, -5, -5, 0.0000000001), control=list(maxit=1000,parscale=c(1000, 0.01, 0.001, 0.01, 0.01, 0.001, 0.1), trace= T, REPORT= 100))
  par3<-result_pupa3$par
  if(result_pupa3$convergence==0){
    print(paste("PUPA", "model 3 converged on rep", k, sep=" "))
    break
  }}



#############################
## model 4, several alphas, no salt and polinators ----
#############################

par4<-c(result_pupa3$par[1], rep(result_pupa3$par[4], times=7), result_pupa3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_pupa4<-optim(par4,compmodel4, method="L-BFGS-B", lower=c(1, rep(0.00001, times=7), 0.0000000001), control=list(maxit=1000,parscale=c(1000, rep(0.01, times=7), 0.1), trace= T, REPORT= 100))
  par4<-result_pupa4$par
  if(result_pupa4$convergence==0){
    print(paste("PUPA", "model 4 converged on rep", k, sep=" "))
    break
  }}



#############################
## model 5, several alphas, salt and polinators ----
#############################

par5<-c(result_pupa3$par[1], result_pupa3$par[2], result_pupa3$par[3], result_pupa4$par[2],
        result_pupa4$par[3], result_pupa4$par[4], result_pupa4$par[5], result_pupa4$par[6], result_pupa4$par[7],
        result_pupa4$par[8], result_pupa3$par[5], result_pupa3$par[6], result_pupa3$par[7])

##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_pupa5<-optim(par5,compmodel5, method="L-BFGS-B", hessian=TRUE, lower=c(1, -5,-5, 0.00001,0.05,0.00001,0.00001,0.05,0.00001,0.00001, -5, -5, 0.0000000001), control=list(maxit=1000,parscale=c(1000, 0.001, 0.01, 0.01, 0.00001, 0.1, 0.01, 0.00001, 0.1, 0.01, 0.001, 0.1 , 0.1), trace= T, REPORT= 100))
  par5<-result_pupa5$par
  if(result_pupa5$convergence==0){
    print(paste("PUPA", "model 5 converged on rep", k, sep=" "))
    break
  }}

inverse <- solve(result_pupa5$hessian)
errors <- sqrt(diag(inverse))  
upper <- result_pupa5$par+1.96*errors
lower <- result_pupa5$par-1.96*errors


#############################
## model 6
#############################

par6<-c(result_pupa3$par[1], result_pupa3$par[2], result_pupa3$par[3], result_pupa4$par[2],
        result_pupa4$par[3], result_pupa4$par[4], result_pupa4$par[5], result_pupa4$par[6], result_pupa4$par[7],
        result_pupa4$par[8], rep(0.0001, times=7), rep(0.0001, times=7), result_pupa3$par[7])


##as before
for(k in 1:25){
  ##now using a specific method that permits constained optimization so that alpha has to be nonzero- 
  result_pupa6<-optim(par6,compmodel6, method="L-BFGS-B", lower=c(1,-5,-5, rep(0.0001, times=7), rep(-5, times=7), rep(-5, times=7), 0.0000000001), control=list(maxit=1000, parscale=c(1000, 0.001, 0.01, 0.01, 0.00001, 0.1, 0.01, 0.00001, 0.1, 0.01, rep(0.0001, times=7), rep(0.0001, times=7), 0.1), trace= T, REPORT= 100))
  par6<-result_pupa6$par  
  if(result_pupa6$convergence==0){
    print(paste("PUPA", "model 6 converged on rep", k, sep=" "))
    break
    
  }}



saveRDS(result_pupa1, file = "results/pupa_results/result_pupa1.rds")
saveRDS(result_pupa2, file = "results/pupa_results/result_pupa2.rds")
saveRDS(result_pupa3, file = "results/pupa_results/result_pupa3.rds")
saveRDS(result_pupa4, file = "results/pupa_results/result_pupa4.rds")
saveRDS(result_pupa5, file = "results/pupa_results/result_pupa5.rds")
saveRDS(result_pupa6, file = "results/pupa_results/result_pupa6.rds")




###Save results
alpha_matrix[6, 1:7]<-par5[4:10]
matrix[6, 1:7]<-par3[1:7]

alpha_matrix<-alpha_matrix[1:6,]

write.csv(alpha_matrix, "alpha_matrix.csv")
write.csv(matrix, "common_effects.csv")

read.csv("alpha_matrix.csv")

read.csv("common_effects.csv")


