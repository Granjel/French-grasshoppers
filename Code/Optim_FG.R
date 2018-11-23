############################################
## Optim analyses --- French grasshoppers ##
############################################
## Rodrigo R. Granjel ## Nov. 2018 #########
############################################



###################################
### OPTIM WITHOUT GRASSHOPPERS: ###
###################################

### predictive function depending on the number of plant species

f_competition <- function(theta, Y, X){
  lambda <- theta[1] #max cover without competition
  alpha <- t(t(theta[2:40])) #as many alphas as plant species + 3 (legumes, grasses, other)
  log_Y_fit <- log(lambda) - log(1 + X %*% alpha) #predictive model
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}


### load dataset:

d <- read.table("Data_FG/FG.txt", header = TRUE, sep = "\t")


### subsample for different dates

d2 <- d[d$time == "1",] #June 2012 ### changing this command changes everything !!!
summary(d2$Focal) #zero ANTODO, GERDIS, TRIFLA and VERPER


### one dataset for each focal species

d_ACHMIL <- d2[d2$Focal == "ACHMIL",]
d_ANTODO <- d2[d2$Focal == "ANTODO",]
d_ARRELA <- d2[d2$Focal == "ARRELA",]
d_BROERE <- d2[d2$Focal == "BROERE",]
d_CENJAC <- d2[d2$Focal == "CENJAC",]
d_CONARV <- d2[d2$Focal == "CONARV",]
d_CREPIS <- d2[d2$Focal == "CREPIS",]
d_DACGLO <- d2[d2$Focal == "DACGLO",]
d_DAUCAR <- d2[d2$Focal == "DAUCAR",]
d_ELYREP <- d2[d2$Focal == "ELYREP",]
d_ERYNGE <- d2[d2$Focal == "ERYNGE",]
d_FESARU <- d2[d2$Focal == "FESARU",]
d_FESRUB <- d2[d2$Focal == "FESRUB",]
d_GALVER <- d2[d2$Focal == "GALVER",]
d_GERDIS <- d2[d2$Focal == "GERDIS",]
d_GERROT <- d2[d2$Focal == "GERROT",]
d_LEUVUL <- d2[d2$Focal == "LEUVUL",]
d_LOLPER <- d2[d2$Focal == "LOLPER",]
d_LOTCOR <- d2[d2$Focal == "LOTCOR",]
d_MEDARA <- d2[d2$Focal == "MEDARA",]
d_ONOREP <- d2[d2$Focal == "ONOREP",]
d_PICECH <- d2[d2$Focal == "PICECH",]
d_PICHIE <- d2[d2$Focal == "PICHIE",]
d_PLALAN <- d2[d2$Focal == "PLALAN",]
d_POAANG <- d2[d2$Focal == "POAANG",]
d_POAPRA <- d2[d2$Focal == "POAPRA",]
d_POATRI <- d2[d2$Focal == "POATRI",]
d_RANACR <- d2[d2$Focal == "RANACR",]
d_RUMACE <- d2[d2$Focal == "RUMACE",]
d_SALPRA <- d2[d2$Focal == "SALPRA",]
d_SONCHU <- d2[d2$Focal == "SONCHU",]
d_TAROFF <- d2[d2$Focal == "TAROFF",]
d_TRIFLA <- d2[d2$Focal == "TRIFLA",]
d_TRIPRA <- d2[d2$Focal == "TRIPRA",]
d_VERBOF <- d2[d2$Focal == "VERBOF",]
d_VERPER <- d2[d2$Focal == "VERPER",]


### X and Y matrices for each plant species:

X_ACHMIL <- as.matrix(d_ACHMIL[, seq(15, 91, by = 2)]) #all the competitors including itself
Y_ACHMIL <- d_ACHMIL$Cover #focal species' cover vector

X_ANTODO <- as.matrix(d_ANTODO[, seq(15, 91, by = 2)])
Y_ANTODO <- d_ANTODO$Cover

X_ARRELA <- as.matrix(d_ARRELA[, seq(15, 91, by = 2)])
Y_ARRELA <- d_ARRELA$Cover

X_BROERE <- as.matrix(d_BROERE[, seq(15, 91, by = 2)])
Y_BROERE <- d_BROERE$Cover

X_CENJAC <- as.matrix(d_CENJAC[, seq(15, 91, by = 2)])
Y_CENJAC <- d_CENJAC$Cover

X_CONARV <- as.matrix(d_CONARV[, seq(15, 91, by = 2)])
Y_CONARV <- d_CONARV$Cover

X_CREPIS <- as.matrix(d_CREPIS[, seq(15, 91, by = 2)])
Y_CREPIS <- d_CREPIS$Cover

X_DACGLO <- as.matrix(d_DACGLO[, seq(15, 91, by = 2)])
Y_DACGLO <- d_DACGLO$Cover

X_DAUCAR <- as.matrix(d_DAUCAR[, seq(15, 91, by = 2)])
Y_DAUCAR <- d_DAUCAR$Cover

X_ELYREP <- as.matrix(d_ELYREP[, seq(15, 91, by = 2)])
Y_ELYREP <- d_ELYREP$Cover 

X_ERYNGE <- as.matrix(d_ERYNGE[, seq(15, 91, by = 2)])
Y_ERYNGE <- d_ERYNGE$Cover

X_FESARU <- as.matrix(d_FESARU[, seq(15, 91, by = 2)])  
Y_FESARU <- d_FESARU$Cover

X_FESRUB <- as.matrix(d_FESRUB[, seq(15, 91, by = 2)])
Y_FESRUB <- d_FESRUB$Cover

X_GALVER <- as.matrix(d_GALVER[, seq(15, 91, by = 2)])
Y_GALVER <- d_GALVER$Cover 

X_GERDIS <- as.matrix(d_GERDIS[, seq(15, 91, by = 2)])
Y_GERDIS <- d_GERDIS$Cover

X_GERROT <- as.matrix(d_GERROT[, seq(15, 91, by = 2)])
Y_GERROT <- d_GERROT$Cover

X_LEUVUL <- as.matrix(d_LEUVUL[, seq(15, 91, by = 2)])
Y_LEUVUL <- d_LEUVUL$Cover

X_LOLPER <- as.matrix(d_LOLPER[, seq(15, 91, by = 2)])
Y_LOLPER <- d_LOLPER$Cover

X_LOTCOR <- as.matrix(d_LOTCOR[, seq(15, 91, by = 2)])
Y_LOTCOR <- d_LOTCOR$Cover

X_MEDARA <- as.matrix(d_MEDARA[, seq(15, 91, by = 2)])
Y_MEDARA <- d_MEDARA$Cover

X_ONOREP <- as.matrix(d_ONOREP[, seq(15, 91, by = 2)])
Y_ONOREP <- d_ONOREP$Cover

X_PICECH <- as.matrix(d_PICECH[, seq(15, 91, by = 2)])
Y_PICECH <- d_PICECH$Cover

X_PICHIE <- as.matrix(d_PICHIE[, seq(15, 91, by = 2)])
Y_PICHIE <- d_PICHIE$Cover

X_PLALAN <- as.matrix(d_PLALAN[, seq(15, 91, by = 2)])
Y_PLALAN <- d_PLALAN$Cover

X_POAANG <- as.matrix(d_POAANG[, seq(15, 91, by = 2)])
Y_POAANG <- d_POAANG$Cover

X_POAPRA <- as.matrix(d_POAPRA[, seq(15, 91, by = 2)])
Y_POAPRA <- d_POAPRA$Cover

X_POATRI <- as.matrix(d_POATRI[, seq(15, 91, by = 2)])
Y_POATRI <- d_POATRI$Cover

X_RANACR <- as.matrix(d_RANACR[, seq(15, 91, by = 2)])
Y_RANACR <- d_RANACR$Cover

X_RUMACE <- as.matrix(d_RUMACE[, seq(15, 91, by = 2)])
Y_RUMACE <- d_RUMACE$Cover

X_SALPRA <- as.matrix(d_SALPRA[, seq(15, 91, by = 2)])
Y_SALPRA <- d_SALPRA$Cover

X_SONCHU <- as.matrix(d_SONCHU[, seq(15, 91, by = 2)])
Y_SONCHU <- d_SONCHU$Cover

X_TAROFF <- as.matrix(d_TAROFF[, seq(15, 91, by = 2)])
Y_TAROFF <- d_TAROFF$Cover

X_TRIFLA <- as.matrix(d_TRIFLA[, seq(15, 91, by = 2)])
Y_TRIFLA <- d_TRIFLA$Cover

X_TRIPRA <- as.matrix(d_TRIPRA[, seq(15, 91, by = 2)])
Y_TRIPRA <- d_TRIPRA$Cover

X_VERBOF <- as.matrix(d_VERBOF[, seq(15, 91, by = 2)])
Y_VERBOF <- d_VERBOF$Cover

X_VERPER <- as.matrix(d_VERPER[, seq(15, 91, by = 2)])
Y_VERPER <- d_VERPER$Cover


### functions for optim, species-specific

f_ACHMIL_optim <- function(theta){f_competition(theta, Y_ACHMIL, X_ACHMIL)}
f_ANTODO_optim <- function(theta){f_competition(theta, Y_ANTODO, X_ANTODO)}
f_ARRELA_optim <- function(theta){f_competition(theta, Y_ARRELA, X_ARRELA)}
f_BROERE_optim <- function(theta){f_competition(theta, Y_BROERE, X_BROERE)}
f_CENJAC_optim <- function(theta){f_competition(theta, Y_CENJAC, X_CENJAC)}
f_CONARV_optim <- function(theta){f_competition(theta, Y_CONARV, X_CONARV)}
f_CREPIS_optim <- function(theta){f_competition(theta, Y_CREPIS, X_CREPIS)}
f_DACGLO_optim <- function(theta){f_competition(theta, Y_DACGLO, X_DACGLO)}
f_DAUCAR_optim <- function(theta){f_competition(theta, Y_DAUCAR, X_DAUCAR)}
f_ELYREP_optim <- function(theta){f_competition(theta, Y_ELYREP, X_ELYREP)}
f_ERYNGE_optim <- function(theta){f_competition(theta, Y_ERYNGE, X_ERYNGE)}
f_FESARU_optim <- function(theta){f_competition(theta, Y_FESARU, X_FESARU)}
f_FESRUB_optim <- function(theta){f_competition(theta, Y_FESRUB, X_FESRUB)}
f_GALVER_optim <- function(theta){f_competition(theta, Y_GALVER, X_GALVER)}
f_GERDIS_optim <- function(theta){f_competition(theta, Y_GERDIS, X_GERDIS)}
f_GERROT_optim <- function(theta){f_competition(theta, Y_GERROT, X_GERROT)}
f_LEUVUL_optim <- function(theta){f_competition(theta, Y_LEUVUL, X_LEUVUL)}
f_LOLPER_optim <- function(theta){f_competition(theta, Y_LOLPER, X_LOLPER)}
f_LOTCOR_optim <- function(theta){f_competition(theta, Y_LOTCOR, X_LOTCOR)}
f_MEDARA_optim <- function(theta){f_competition(theta, Y_MEDARA, X_MEDARA)}
f_ONOREP_optim <- function(theta){f_competition(theta, Y_ONOREP, X_ONOREP)}
f_PICECH_optim <- function(theta){f_competition(theta, Y_PICECH, X_PICECH)}
f_PICHIE_optim <- function(theta){f_competition(theta, Y_PICHIE, X_PICHIE)}
f_PLALAN_optim <- function(theta){f_competition(theta, Y_PLALAN, X_PLALAN)}
f_POAANG_optim <- function(theta){f_competition(theta, Y_POAANG, X_POAANG)}
f_POAPRA_optim <- function(theta){f_competition(theta, Y_POAPRA, X_POAPRA)}
f_POATRI_optim <- function(theta){f_competition(theta, Y_POATRI, X_POATRI)}
f_RANACR_optim <- function(theta){f_competition(theta, Y_RANACR, X_RANACR)}
f_RUMACE_optim <- function(theta){f_competition(theta, Y_RUMACE, X_RUMACE)}
f_SALPRA_optim <- function(theta){f_competition(theta, Y_SALPRA, X_SALPRA)}
f_SONCHU_optim <- function(theta){f_competition(theta, Y_SONCHU, X_SONCHU)}
f_TAROFF_optim <- function(theta){f_competition(theta, Y_TAROFF, X_TAROFF)}
f_TRIFLA_optim <- function(theta){f_competition(theta, Y_TRIFLA, X_TRIFLA)}
f_TRIPRA_optim <- function(theta){f_competition(theta, Y_TRIPRA, X_TRIPRA)}
f_VERBOF_optim <- function(theta){f_competition(theta, Y_VERBOF, X_VERBOF)}
f_VERPER_optim <- function(theta){f_competition(theta, Y_VERPER, X_VERPER)}


### optim wrap for the diff. species ---:

ini <- rep(1, 40) #initial values for all species
low <- rep(0, 40) #lower values for all species

out_ACHMIL <- optim(ini, f_ACHMIL_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_ANTODO <- optim(ini, f_ANTODO_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_ARRELA <- optim(ini, f_ARRELA_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_BROERE <- optim(ini, f_BROERE_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_CENJAC <- optim(ini, f_CENJAC_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_CONARV <- optim(ini, f_CONARV_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_CREPIS <- optim(ini, f_CREPIS_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_DACGLO <- optim(ini, f_DACGLO_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_DAUCAR <- optim(ini, f_DAUCAR_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_ELYREP <- optim(ini, f_ELYREP_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_ERYNGE <- optim(ini, f_ERYNGE_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_FESARU <- optim(ini, f_FESARU_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_FESRUB <- optim(ini, f_FESRUB_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_GALVER <- optim(ini, f_GALVER_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_GERDIS <- optim(ini, f_GERDIS_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_GERROT <- optim(ini, f_GERROT_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_LEUVUL <- optim(ini, f_LEUVUL_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_LOLPER <- optim(ini, f_LOLPER_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_LOTCOR <- optim(ini, f_LOTCOR_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_MEDARA <- optim(ini, f_MEDARA_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_ONOREP <- optim(ini, f_ONOREP_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_PICECH <- optim(ini, f_PICECH_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_PICHIE <- optim(ini, f_PICHIE_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_PLALAN <- optim(ini, f_PLALAN_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_POAANG <- optim(ini, f_POAANG_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_POAPRA <- optim(ini, f_POAPRA_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_POATRI <- optim(ini, f_POATRI_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_RANACR <- optim(ini, f_RANACR_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_RUMACE <- optim(ini, f_RUMACE_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_SALPRA <- optim(ini, f_SALPRA_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_SONCHU <- optim(ini, f_SONCHU_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_TAROFF <- optim(ini, f_TAROFF_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_TRIFLA <- optim(ini, f_TRIFLA_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_TRIPRA <- optim(ini, f_TRIPRA_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_VERBOF <- optim(ini, f_VERBOF_optim, lower = low, method = 'L-BFGS-B', hessian = T)
out_VERPER <- optim(ini, f_VERPER_optim, lower = low, method = 'L-BFGS-B', hessian = T)


### extracting and saving the interesting parameters (lambda and alpha):

lambda <- c(out_ACHMIL$par[1], out_ANTODO$par[1], out_ARRELA$par[1], out_BROERE$par[1], out_CENJAC$par[1], out_CONARV$par[1],
            out_CREPIS$par[1], out_DACGLO$par[1], out_DAUCAR$par[1], out_ELYREP$par[1], out_ERYNGE$par[1], out_FESARU$par[1],
            out_FESRUB$par[1], out_GALVER$par[1], out_GERDIS$par[1], out_GERROT$par[1], out_LEUVUL$par[1], out_LOLPER$par[1],
            out_LOTCOR$par[1], out_MEDARA$par[1], out_ONOREP$par[1], out_PICECH$par[1], out_PICHIE$par[1], out_PLALAN$par[1],
            out_POAANG$par[1], out_POAPRA$par[1], out_POATRI$par[1], out_RANACR$par[1], out_RUMACE$par[1], out_SALPRA$par[1],
            out_SONCHU$par[1], out_TAROFF$par[1], out_TRIFLA$par[1], out_TRIPRA$par[1], out_VERBOF$par[1], out_VERPER$par[1])
alpha <- rbind(out_ACHMIL$par[2:37], out_ANTODO$par[2:37], out_ARRELA$par[2:37], out_BROERE$par[2:37], out_CENJAC$par[2:37], out_CONARV$par[2:37],
               out_CREPIS$par[2:37], out_DACGLO$par[2:37], out_DAUCAR$par[2:37], out_ELYREP$par[2:37], out_ERYNGE$par[2:37], out_FESARU$par[2:37],
               out_FESRUB$par[2:37], out_GALVER$par[2:37], out_GERDIS$par[2:37], out_GERROT$par[2:37], out_LEUVUL$par[2:37], out_LOLPER$par[2:37],
               out_LOTCOR$par[2:37], out_MEDARA$par[2:37], out_ONOREP$par[2:37], out_PICECH$par[2:37], out_PICHIE$par[2:37], out_PLALAN$par[2:37],
               out_POAANG$par[2:37], out_POAPRA$par[2:37], out_POATRI$par[2:37], out_RANACR$par[2:37], out_RUMACE$par[2:37], out_SALPRA$par[2:37],
               out_SONCHU$par[2:37], out_TAROFF$par[2:37], out_TRIFLA$par[2:37], out_TRIPRA$par[2:37], out_VERBOF$par[2:37], out_VERPER$par[2:37])


### naming the rows and columns of the objects lambda and alpha

names(lambda) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                   "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                   "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")

rownames(alpha) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                     "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                     "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")

colnames(alpha) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                     "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                     "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")



### removing ANTODO, GERDIS, TRIFLA and VERPER (WARNING: depending on the date selected!)

lambda <- c(lambda[1], lambda[3:14], lambda[16:32], lambda[34:35])
alpha <- alpha[c(-2, -15, -33, -36),]
alpha <- alpha[, c(-2, -15, -33, -36)]


### save results:

write.table(lambda, file = "Results/lambda_1.txt", sep = "\t", row.names = TRUE)
write.table(alpha, file = "Results/alpha_1.txt", sep = "\t", row.names = FALSE)


### clean environment:

#rm(list=ls())



################################
### OPTIM WITH GRASSHOPPERS: ###
################################

### predictive functions depending on the number of plant and grasshopper species:

f_ACHMIL_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] #max cover without competition
  alpha <- t(t(theta[2:40])) #as many alphas as plant species
  gamma <- matrix(0, nrow = 6, ncol = 1) #grasshopper effects on plants (as many rows as grasshopper species)
  gamma[1:6] <- theta[41:46] #as many thetas as grasshopper species
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha + X_Cb_gh %*% gamma)  #competition function #X_plant = competitive effects between plants
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_ANTODO_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_ARRELA_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_BROERE_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_CENJAC_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_CONARV_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_CREPIS_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_DACGLO_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_DAUCAR_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_ELYREP_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_ERYNGE_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_FESARU_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_FESRUB_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_GALVER_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_GERDIS_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_GERROT_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_LEUVUL_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_LOLPER_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_LOTCOR_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_MEDARA_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_ONOREP_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_PICECH_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_PICHIE_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_PLALAN_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_POAANG_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_POAPRA_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_POATRI_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_RANACR_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_RUMACE_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_SALPRA_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_SONCHU_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_TAROFF_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_TRIFLA_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_TRIPRA_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_VERBOF_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_VERPER_gh <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:40])) 
  gamma <- matrix(0, nrow = 6, ncol = 1) 
  gamma[1:6] <- theta[41:46] 
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) - log(1 + X_gh %*% gamma) 
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
d_ANTODO_gh <- d2[d2$Focal == "ANTODO",]
d_ARRELA_gh <- d2[d2$Focal == "ARRELA",]
d_BROERE_gh <- d2[d2$Focal == "BROERE",]
d_CENJAC_gh <- d2[d2$Focal == "CENJAC",]
d_CONARV_gh <- d2[d2$Focal == "CONARV",]
d_CREPIS_gh <- d2[d2$Focal == "CREPIS",]
d_DACGLO_gh <- d2[d2$Focal == "DACGLO",]
d_DAUCAR_gh <- d2[d2$Focal == "DAUCAR",]
d_ELYREP_gh <- d2[d2$Focal == "ELYREP",]
d_ERYNGE_gh <- d2[d2$Focal == "ERYNGE",]
d_FESARU_gh <- d2[d2$Focal == "FESARU",]
d_FESRUB_gh <- d2[d2$Focal == "FESRUB",]
d_GALVER_gh <- d2[d2$Focal == "GALVER",]
d_GERDIS_gh <- d2[d2$Focal == "GERDIS",]
d_GERROT_gh <- d2[d2$Focal == "GERROT",]
d_LEUVUL_gh <- d2[d2$Focal == "LEUVUL",]
d_LOLPER_gh <- d2[d2$Focal == "LOLPER",]
d_LOTCOR_gh <- d2[d2$Focal == "LOTCOR",]
d_MEDARA_gh <- d2[d2$Focal == "MEDARA",]
d_ONOREP_gh <- d2[d2$Focal == "ONOREP",]
d_PICECH_gh <- d2[d2$Focal == "PICECH",]
d_PICHIE_gh <- d2[d2$Focal == "PICHIE",]
d_PLALAN_gh <- d2[d2$Focal == "PLALAN",]
d_POAANG_gh <- d2[d2$Focal == "POAANG",]
d_POAPRA_gh <- d2[d2$Focal == "POAPRA",]
d_POATRI_gh <- d2[d2$Focal == "POATRI",]
d_RANACR_gh <- d2[d2$Focal == "RANACR",]
d_RUMACE_gh <- d2[d2$Focal == "RUMACE",]
d_SALPRA_gh <- d2[d2$Focal == "SALPRA",]
d_SONCHU_gh <- d2[d2$Focal == "SONCHU",]
d_TAROFF_gh <- d2[d2$Focal == "TAROFF",]
d_TRIFLA_gh <- d2[d2$Focal == "TRIFLA",]
d_TRIPRA_gh <- d2[d2$Focal == "TRIPRA",]
d_VERBOF_gh <- d2[d2$Focal == "VERBOF",]
d_VERPER_gh <- d2[d2$Focal == "VERPER",]


### X_plant, X_gh and Y matrices for each plant species:

X_ACHMIL_plant <- as.matrix(d_ACHMIL_gh[, seq(15, 91, by = 2)]) #competition
X_ACHMIL_gh <- as.matrix(d_ACHMIL_gh[6:11]) #depends on the number of links for every species (suppl. mat. Gross)
Y_ACHMIL_gh <- d_ACHMIL_gh$Cover #cover

X_ANTODO_plant <- as.matrix(d_ANTODO_gh[, seq(15, 91, by = 2)])
X_ANTODO_gh <- as.matrix(d_ANTODO_gh[6:11])
Y_ANTODO_gh <- d_ANTODO_gh$Cover

X_ARRELA_plant <- as.matrix(d_ARRELA_gh[, seq(15, 91, by = 2)])
X_ARRELA_gh <- as.matrix(d_ARRELA_gh[6:11])
Y_ARRELA_gh <- d_ARRELA_gh$Cover

X_BROERE_plant <- as.matrix(d_BROERE_gh[, seq(15, 91, by = 2)])
X_BROERE_gh <- as.matrix(d_BROERE_gh[6:11])
Y_BROERE_gh <- d_BROERE_gh$Cover

X_CENJAC_plant <- as.matrix(d_CENJAC_gh[, seq(15, 91, by = 2)])
X_CENJAC_gh <- as.matrix(d_CENJAC_gh[6:11])
Y_CENJAC_gh <- d_CENJAC_gh$Cover

X_CONARV_plant <- as.matrix(d_CONARV_gh[, seq(15, 91, by = 2)])
X_CONARV_gh <- as.matrix(d_CONARV_gh[6:11])
Y_CONARV_gh <- d_CONARV_gh$Cover

X_CREPIS_plant <- as.matrix(d_CREPIS_gh[, seq(15, 91, by = 2)])
X_CREPIS_gh <- as.matrix(d_CREPIS_gh[6:11])
Y_CREPIS_gh <- d_CREPIS_gh$Cover

X_DACGLO_plant <- as.matrix(d_DACGLO_gh[, seq(15, 91, by = 2)])
X_DACGLO_gh <- as.matrix(d_DACGLO_gh[6:11])
Y_DACGLO_gh <- d_DACGLO_gh$Cover

X_DAUCAR_plant <- as.matrix(d_DAUCAR_gh[, seq(15, 91, by = 2)])
X_DAUCAR_gh <- as.matrix(d_DAUCAR_gh[6:11])
Y_DAUCAR_gh <- d_DAUCAR_gh$Cover

X_ELYREP_plant <- as.matrix(d_ELYREP_gh[, seq(15, 91, by = 2)])
X_ELYREP_gh <- as.matrix(d_ELYREP_gh[6:11])
Y_ELYREP_gh <- d_ELYREP_gh$Cover

X_ERYNGE_plant <- as.matrix(d_ERYNGE_gh[, seq(15, 91, by = 2)])
X_ERYNGE_gh <- as.matrix(d_ERYNGE_gh[6:11])
Y_ERYNGE_gh <- d_ERYNGE_gh$Cover

X_FESARU_plant <- as.matrix(d_FESARU_gh[, seq(15, 91, by = 2)])
X_FESARU_gh <- as.matrix(d_FESARU_gh[6:11])
Y_FESARU_gh <- d_FESARU_gh$Cover

X_FESRUB_plant <- as.matrix(d_FESRUB_gh[, seq(15, 91, by = 2)])
X_FESRUB_gh <- as.matrix(d_FESRUB_gh[6:11])
Y_FESRUB_gh <- d_FESRUB_gh$Cover

X_GALVER_plant <- as.matrix(d_GALVER_gh[, seq(15, 91, by = 2)])
X_GALVER_gh <- as.matrix(d_GALVER_gh[6:11])
Y_GALVER_gh <- d_GALVER_gh$Cover

X_GERDIS_plant <- as.matrix(d_GERDIS_gh[, seq(15, 91, by = 2)])
X_GERDIS_gh <- as.matrix(d_GERDIS_gh[6:11])
Y_GERDIS_gh <- d_GERDIS_gh$Cover

X_GERROT_plant <- as.matrix(d_GERROT_gh[, seq(15, 91, by = 2)])
X_GERROT_gh <- as.matrix(d_GERROT_gh[6:11])
Y_GERROT_gh <- d_GERROT_gh$Cover

X_LEUVUL_plant <- as.matrix(d_LEUVUL_gh[, seq(15, 91, by = 2)])
X_LEUVUL_gh <- as.matrix(d_LEUVUL_gh[6:11])
Y_LEUVUL_gh <- d_LEUVUL_gh$Cover

X_LOLPER_plant <- as.matrix(d_LOLPER_gh[, seq(15, 91, by = 2)])
X_LOLPER_gh <- as.matrix(d_LOLPER_gh[6:11])
Y_LOLPER_gh <- d_LOLPER_gh$Cover

X_LOTCOR_plant <- as.matrix(d_LOTCOR_gh[, seq(15, 91, by = 2)])
X_LOTCOR_gh <- as.matrix(d_LOTCOR_gh[6:11])
Y_LOTCOR_gh <- d_LOTCOR_gh$Cover

X_MEDARA_plant <- as.matrix(d_MEDARA_gh[, seq(15, 91, by = 2)])
X_MEDARA_gh <- as.matrix(d_MEDARA_gh[6:11])
Y_MEDARA_gh <- d_MEDARA_gh$Cover

X_ONOREP_plant <- as.matrix(d_ONOREP_gh[, seq(15, 91, by = 2)])
X_ONOREP_gh <- as.matrix(d_ONOREP_gh[6:11])
Y_ONOREP_gh <- d_ONOREP_gh$Cover

X_PICECH_plant <- as.matrix(d_PICECH_gh[, seq(15, 91, by = 2)])
X_PICECH_gh <- as.matrix(d_PICECH_gh[6:11])
Y_PICECH_gh <- d_PICECH_gh$Cover

X_PICHIE_plant <- as.matrix(d_PICHIE_gh[, seq(15, 91, by = 2)])
X_PICHIE_gh <- as.matrix(d_PICHIE_gh[6:11])
Y_PICHIE_gh <- d_PICHIE_gh$Cover

X_PLALAN_plant <- as.matrix(d_PLALAN_gh[, seq(15, 91, by = 2)])
X_PLALAN_gh <- as.matrix(d_PLALAN_gh[6:11])
Y_PLALAN_gh <- d_PLALAN_gh$Cover

X_POAANG_plant <- as.matrix(d_POAANG_gh[, seq(15, 91, by = 2)])
X_POAANG_gh <- as.matrix(d_POAANG_gh[6:11])
Y_POAANG_gh <- d_POAANG_gh$Cover

X_POAPRA_plant <- as.matrix(d_POAPRA_gh[, seq(15, 91, by = 2)])
X_POAPRA_gh <- as.matrix(d_POAPRA_gh[6:11])
Y_POAPRA_gh <- d_POAPRA_gh$Cover

X_POATRI_plant <- as.matrix(d_POATRI_gh[, seq(15, 91, by = 2)])
X_POATRI_gh <- as.matrix(d_POATRI_gh[6:11])
Y_POATRI_gh <- d_POATRI_gh$Cover

X_RANACR_plant <- as.matrix(d_RANACR_gh[, seq(15, 91, by = 2)])
X_RANACR_gh <- as.matrix(d_RANACR_gh[6:11])
Y_RANACR_gh <- d_RANACR_gh$Cover

X_RUMACE_plant <- as.matrix(d_RUMACE_gh[, seq(15, 91, by = 2)])
X_RUMACE_gh <- as.matrix(d_RUMACE_gh[6:11])
Y_RUMACE_gh <- d_RUMACE_gh$Cover

X_SALPRA_plant <- as.matrix(d_SALPRA_gh[, seq(15, 91, by = 2)])
X_SALPRA_gh <- as.matrix(d_SALPRA_gh[6:11])
Y_SALPRA_gh <- d_SALPRA_gh$Cover

X_SONCHU_plant <- as.matrix(d_SONCHU_gh[, seq(15, 91, by = 2)])
X_SONCHU_gh <- as.matrix(d_SONCHU_gh[6:11])
Y_SONCHU_gh <- d_SONCHU_gh$Cover

X_TAROFF_plant <- as.matrix(d_TAROFF_gh[, seq(15, 91, by = 2)])
X_TAROFF_gh <- as.matrix(d_TAROFF_gh[6:11])
Y_TAROFF_gh <- d_TAROFF_gh$Cover

X_TRIFLA_plant <- as.matrix(d_TRIFLA_gh[, seq(15, 91, by = 2)])
X_TRIFLA_gh <- as.matrix(d_TRIFLA_gh[6:11])
Y_TRIFLA_gh <- d_TRIFLA_gh$Cover

X_TRIPRA_plant <- as.matrix(d_TRIPRA_gh[, seq(15, 91, by = 2)])
X_TRIPRA_gh <- as.matrix(d_TRIPRA_gh[6:11])
Y_TRIPRA_gh <- d_TRIPRA_gh$Cover

X_VERBOF_plant <- as.matrix(d_VERBOF_gh[, seq(15, 91, by = 2)])
X_VERBOF_gh <- as.matrix(d_VERBOF_gh[6:11])
Y_VERBOF_gh <- d_VERBOF_gh$Cover

X_VERPER_plant <- as.matrix(d_VERPER_gh[, seq(15, 91, by = 2)])
X_VERPER_gh <- as.matrix(d_VERPER_gh[6:11])
Y_VERPER_gh <- d_VERPER_gh$Cover


### functions for optim, species-specific

f_ACHMIL_gh_optim <- function(theta){f_ACHMIL_gh(theta, Y_ACHMIL_gh, X_ACHMIL_plant, X_ACHMIL_gh)}
f_ANTODO_gh_optim <- function(theta){f_ANTODO_gh(theta, Y_ANTODO_gh, X_ANTODO_plant, X_ANTODO_gh)}
f_ARRELA_gh_optim <- function(theta){f_ARRELA_gh(theta, Y_ARRELA_gh, X_ARRELA_plant, X_ARRELA_gh)}
f_BROERE_gh_optim <- function(theta){f_BROERE_gh(theta, Y_BROERE_gh, X_BROERE_plant, X_BROERE_gh)}
f_CENJAC_gh_optim <- function(theta){f_CENJAC_gh(theta, Y_CENJAC_gh, X_CENJAC_plant, X_CENJAC_gh)}
f_CONARV_gh_optim <- function(theta){f_CONARV_gh(theta, Y_CONARV_gh, X_CONARV_plant, X_CONARV_gh)}
f_CREPIS_gh_optim <- function(theta){f_CREPIS_gh(theta, Y_CREPIS_gh, X_CREPIS_plant, X_CREPIS_gh)}
f_DACGLO_gh_optim <- function(theta){f_DACGLO_gh(theta, Y_DACGLO_gh, X_DACGLO_plant, X_DACGLO_gh)}
f_DAUCAR_gh_optim <- function(theta){f_DAUCAR_gh(theta, Y_DAUCAR_gh, X_DAUCAR_plant, X_DAUCAR_gh)}
f_ELYREP_gh_optim <- function(theta){f_ELYREP_gh(theta, Y_ELYREP_gh, X_ELYREP_plant, X_ELYREP_gh)}
f_ERYNGE_gh_optim <- function(theta){f_ERYNGE_gh(theta, Y_ERYNGE_gh, X_ERYNGE_plant, X_ERYNGE_gh)}
f_FESARU_gh_optim <- function(theta){f_FESARU_gh(theta, Y_FESARU_gh, X_FESARU_plant, X_FESARU_gh)}
f_FESRUB_gh_optim <- function(theta){f_FESRUB_gh(theta, Y_FESRUB_gh, X_FESRUB_plant, X_FESRUB_gh)}
f_GALVER_gh_optim <- function(theta){f_GALVER_gh(theta, Y_GALVER_gh, X_GALVER_plant, X_GALVER_gh)}
f_GERDIS_gh_optim <- function(theta){f_GERDIS_gh(theta, Y_GERDIS_gh, X_GERDIS_plant, X_GERDIS_gh)}
f_GERROT_gh_optim <- function(theta){f_GERROT_gh(theta, Y_GERROT_gh, X_GERROT_plant, X_GERROT_gh)}
f_LEUVUL_gh_optim <- function(theta){f_LEUVUL_gh(theta, Y_LEUVUL_gh, X_LEUVUL_plant, X_LEUVUL_gh)}
f_LOLPER_gh_optim <- function(theta){f_LOLPER_gh(theta, Y_LOLPER_gh, X_LOLPER_plant, X_LOLPER_gh)}
f_LOTCOR_gh_optim <- function(theta){f_LOTCOR_gh(theta, Y_LOTCOR_gh, X_LOTCOR_plant, X_LOTCOR_gh)}
f_MEDARA_gh_optim <- function(theta){f_MEDARA_gh(theta, Y_MEDARA_gh, X_MEDARA_plant, X_MEDARA_gh)}
f_ONOREP_gh_optim <- function(theta){f_ONOREP_gh(theta, Y_ONOREP_gh, X_ONOREP_plant, X_ONOREP_gh)}
f_PICECH_gh_optim <- function(theta){f_PICECH_gh(theta, Y_PICECH_gh, X_PICECH_plant, X_PICECH_gh)}
f_PICHIE_gh_optim <- function(theta){f_PICHIE_gh(theta, Y_PICHIE_gh, X_PICHIE_plant, X_PICHIE_gh)}
f_PLALAN_gh_optim <- function(theta){f_PLALAN_gh(theta, Y_PLALAN_gh, X_PLALAN_plant, X_PLALAN_gh)}
f_POAANG_gh_optim <- function(theta){f_POAANG_gh(theta, Y_POAANG_gh, X_POAANG_plant, X_POAANG_gh)}
f_POAPRA_gh_optim <- function(theta){f_POAPRA_gh(theta, Y_POAPRA_gh, X_POAPRA_plant, X_POAPRA_gh)}
f_POATRI_gh_optim <- function(theta){f_POATRI_gh(theta, Y_POATRI_gh, X_POATRI_plant, X_POATRI_gh)}
f_RANACR_gh_optim <- function(theta){f_RANACR_gh(theta, Y_RANACR_gh, X_RANACR_plant, X_RANACR_gh)}
f_RUMACE_gh_optim <- function(theta){f_RUMACE_gh(theta, Y_RUMACE_gh, X_RUMACE_plant, X_RUMACE_gh)}
f_SALPRA_gh_optim <- function(theta){f_SALPRA_gh(theta, Y_SALPRA_gh, X_SALPRA_plant, X_SALPRA_gh)}
f_SONCHU_gh_optim <- function(theta){f_SONCHU_gh(theta, Y_SONCHU_gh, X_SONCHU_plant, X_SONCHU_gh)}
f_TAROFF_gh_optim <- function(theta){f_TAROFF_gh(theta, Y_TAROFF_gh, X_TAROFF_plant, X_TAROFF_gh)}
f_TRIFLA_gh_optim <- function(theta){f_TRIFLA_gh(theta, Y_TRIFLA_gh, X_TRIFLA_plant, X_TRIFLA_gh)}
f_TRIPRA_gh_optim <- function(theta){f_TRIPRA_gh(theta, Y_TRIPRA_gh, X_TRIPRA_plant, X_TRIPRA_gh)}
f_VERBOF_gh_optim <- function(theta){f_VERBOF_gh(theta, Y_VERBOF_gh, X_VERBOF_plant, X_VERBOF_gh)}
f_VERPER_gh_optim <- function(theta){f_VERPER_gh(theta, Y_VERPER_gh, X_VERPER_plant, X_VERPER_gh)}


### optim wrap for the diff. species (one 'ini' and 'lower' for each spp.) ---:

ini_ACHMIL_gh <- c(1, rep(0.000001, 45))
lower_ACHMIL_gh <- rep(0.000001, 46)
out_ACHMIL_gh <- optim(ini_ACHMIL_gh, f_ACHMIL_gh_optim, lower = lower_ACHMIL_gh, method = 'L-BFGS-B', hessian = T)

ini_ANTODO_gh <- c(1, rep(0.000001, 45))
lower_ANTODO_gh <- rep(0.000001, 46)
out_ANTODO_gh <- optim(ini_ANTODO_gh, f_ANTODO_gh_optim, lower = lower_ANTODO_gh, method = 'L-BFGS-B', hessian = T)

ini_ARRELA_gh <- c(1, rep(0.000001, 45))
lower_ARRELA_gh <- rep(0.000001, 46)
out_ARRELA_gh <- optim(ini_ARRELA_gh, f_ARRELA_gh_optim, lower = lower_ARRELA_gh, method = 'L-BFGS-B', hessian = T)

ini_BROERE_gh <- c(1, rep(0.000001, 45))
lower_BROERE_gh <- rep(0.000001, 46)
out_BROERE_gh <- optim(ini_BROERE_gh, f_BROERE_gh_optim, lower = lower_BROERE_gh, method = 'L-BFGS-B', hessian = T)

ini_CENJAC_gh <- c(1, rep(0.000001, 45))
lower_CENJAC_gh <- rep(0.000001, 46)
out_CENJAC_gh <- optim(ini_CENJAC_gh, f_CENJAC_gh_optim, lower = lower_CENJAC_gh, method = 'L-BFGS-B', hessian = T)

ini_CONARV_gh <- c(1, rep(0.000001, 45))
lower_CONARV_gh <- rep(0.000001, 46)
out_CONARV_gh <- optim(ini_CONARV_gh, f_CONARV_gh_optim, lower = lower_CONARV_gh, method = 'L-BFGS-B', hessian = T)

ini_CREPIS_gh <- c(1, rep(0.000001, 45))
lower_CREPIS_gh <- rep(0.000001, 46)
out_CREPIS_gh <- optim(ini_CREPIS_gh, f_CREPIS_gh_optim, lower = lower_CREPIS_gh, method = 'L-BFGS-B', hessian = T)

ini_DACGLO_gh <- c(1, rep(0.000001, 45))
lower_DACGLO_gh <- rep(0.000001, 46)
out_DACGLO_gh <- optim(ini_DACGLO_gh, f_DACGLO_gh_optim, lower = lower_DACGLO_gh, method = 'L-BFGS-B', hessian = T)

ini_DAUCAR_gh <- c(1, rep(0.000001, 45))
lower_DAUCAR_gh <- rep(0.000001, 46)
out_DAUCAR_gh <- optim(ini_DAUCAR_gh, f_DAUCAR_gh_optim, lower = lower_DAUCAR_gh, method = 'L-BFGS-B', hessian = T)

ini_ELYREP_gh <- c(1, rep(0.000001, 45))
lower_ELYREP_gh <- rep(0.000001, 46)
out_ELYREP_gh <- optim(ini_ELYREP_gh, f_ELYREP_gh_optim, lower = lower_ELYREP_gh, method = 'L-BFGS-B', hessian = T)

ini_ERYNGE_gh <- c(1, rep(0.000001, 45))
lower_ERYNGE_gh <- rep(0.000001, 46)
out_ERYNGE_gh <- optim(ini_ERYNGE_gh, f_ERYNGE_gh_optim, lower = lower_ERYNGE_gh, method = 'L-BFGS-B', hessian = T)

ini_FESARU_gh <- c(1, rep(0.000001, 45))
lower_FESARU_gh <- rep(0.000001, 46)
out_FESARU_gh <- optim(ini_FESARU_gh, f_FESARU_gh_optim, lower = lower_FESARU_gh, method = 'L-BFGS-B', hessian = T)

ini_FESRUB_gh <- c(1, rep(0.000001, 45))
lower_FESRUB_gh <- rep(0.000001, 46)
out_FESRUB_gh <- optim(ini_FESRUB_gh, f_FESRUB_gh_optim, lower = lower_FESRUB_gh, method = 'L-BFGS-B', hessian = T)

ini_GALVER_gh <- c(1, rep(0.000001, 45))
lower_GALVER_gh <- rep(0.000001, 46)
out_GALVER_gh <- optim(ini_GALVER_gh, f_GALVER_gh_optim, lower = lower_GALVER_gh, method = 'L-BFGS-B', hessian = T)

ini_GERDIS_gh <- c(1, rep(0.000001, 45))
lower_GERDIS_gh <- rep(0.000001, 46)
out_GERDIS_gh <- optim(ini_GERDIS_gh, f_GERDIS_gh_optim, lower = lower_GERDIS_gh, method = 'L-BFGS-B', hessian = T)

ini_GERROT_gh <- c(1, rep(0.000001, 45))
lower_GERROT_gh <- rep(0.000001, 46)
out_GERROT_gh <- optim(ini_GERROT_gh, f_GERROT_gh_optim, lower = lower_GERROT_gh, method = 'L-BFGS-B', hessian = T)

ini_LEUVUL_gh <- c(1, rep(0.000001, 45))
lower_LEUVUL_gh <- rep(0.000001, 46)
out_LEUVUL_gh <- optim(ini_LEUVUL_gh, f_LEUVUL_gh_optim, lower = lower_LEUVUL_gh, method = 'L-BFGS-B', hessian = T)

ini_LOLPER_gh <- c(1, rep(0.000001, 45))
lower_LOLPER_gh <- rep(0.000001, 46)
out_LOLPER_gh <- optim(ini_LOLPER_gh, f_LOLPER_gh_optim, lower = lower_LOLPER_gh, method = 'L-BFGS-B', hessian = T)

ini_LOTCOR_gh <- c(1, rep(0.000001, 45))
lower_LOTCOR_gh <- rep(0.000001, 46)
out_LOTCOR_gh <- optim(ini_LOTCOR_gh, f_LOTCOR_gh_optim, lower = lower_LOTCOR_gh, method = 'L-BFGS-B', hessian = T)

ini_MEDARA_gh <- c(1, rep(0.000001, 45))
lower_MEDARA_gh <- rep(0.000001, 46)
out_MEDARA_gh <- optim(ini_MEDARA_gh, f_MEDARA_gh_optim, lower = lower_MEDARA_gh, method = 'L-BFGS-B', hessian = T)

ini_ONOREP_gh <- c(1, rep(0.000001, 45))
lower_ONOREP_gh <- rep(0.000001, 46)
out_ONOREP_gh <- optim(ini_ONOREP_gh, f_ONOREP_gh_optim, lower = lower_ONOREP_gh, method = 'L-BFGS-B', hessian = T)

ini_PICECH_gh <- c(1, rep(0.000001, 45))
lower_PICECH_gh <- rep(0.000001, 46)
out_PICECH_gh <- optim(ini_PICECH_gh, f_PICECH_gh_optim, lower = lower_PICECH_gh, method = 'L-BFGS-B', hessian = T)

ini_PICHIE_gh <- c(1, rep(0.000001, 45))
lower_PICHIE_gh <- rep(0.000001, 46)
out_PICHIE_gh <- optim(ini_PICHIE_gh, f_PICHIE_gh_optim, lower = lower_PICHIE_gh, method = 'L-BFGS-B', hessian = T)

ini_PLALAN_gh <- c(1, rep(0.000001, 45))
lower_PLALAN_gh <- rep(0.000001, 46)
out_PLALAN_gh <- optim(ini_PLALAN_gh, f_PLALAN_gh_optim, lower = lower_PLALAN_gh, method = 'L-BFGS-B', hessian = T)

ini_POAANG_gh <- c(1, rep(0.000001, 45))
lower_POAANG_gh <- rep(0.000001, 46)
out_POAANG_gh <- optim(ini_POAANG_gh, f_POAANG_gh_optim, lower = lower_POAANG_gh, method = 'L-BFGS-B', hessian = T)

ini_POAPRA_gh <- c(1, rep(0.000001, 45))
lower_POAPRA_gh <- rep(0.000001, 46)
out_POAPRA_gh <- optim(ini_POAPRA_gh, f_POAPRA_gh_optim, lower = lower_POAPRA_gh, method = 'L-BFGS-B', hessian = T)

ini_POATRI_gh <- c(1, rep(0.000001, 45))
lower_POATRI_gh <- rep(0.000001, 46)
out_POATRI_gh <- optim(ini_POATRI_gh, f_POATRI_gh_optim, lower = lower_POATRI_gh, method = 'L-BFGS-B', hessian = T)

ini_RANACR_gh <- c(1, rep(0.000001, 45))
lower_RANACR_gh <- rep(0.000001, 46)
out_RANACR_gh <- optim(ini_RANACR_gh, f_RANACR_gh_optim, lower = lower_RANACR_gh, method = 'L-BFGS-B', hessian = T)

ini_RUMACE_gh <- c(1, rep(0.000001, 45))
lower_RUMACE_gh <- rep(0.000001, 46)
out_RUMACE_gh <- optim(ini_RUMACE_gh, f_RUMACE_gh_optim, lower = lower_RUMACE_gh, method = 'L-BFGS-B', hessian = T)

ini_SALPRA_gh <- c(1, rep(0.000001, 45))
lower_SALPRA_gh <- rep(0.000001, 46)
out_SALPRA_gh <- optim(ini_SALPRA_gh, f_SALPRA_gh_optim, lower = lower_SALPRA_gh, method = 'L-BFGS-B', hessian = T)

ini_SONCHU_gh <- c(1, rep(0.000001, 45))
lower_SONCHU_gh <- rep(0.000001, 46)
out_SONCHU_gh <- optim(ini_SONCHU_gh, f_SONCHU_gh_optim, lower = lower_SONCHU_gh, method = 'L-BFGS-B', hessian = T)

ini_TAROFF_gh <- c(1, rep(0.000001, 45))
lower_TAROFF_gh <- rep(0.000001, 46)
out_TAROFF_gh <- optim(ini_TAROFF_gh, f_TAROFF_gh_optim, lower = lower_TAROFF_gh, method = 'L-BFGS-B', hessian = T)

ini_TRIFLA_gh <- c(1, rep(0.000001, 45))
lower_TRIFLA_gh <- rep(0.000001, 46)
out_TRIFLA_gh <- optim(ini_TRIFLA_gh, f_TRIFLA_gh_optim, lower = lower_TRIFLA_gh, method = 'L-BFGS-B', hessian = T)

ini_TRIPRA_gh <- c(1, rep(0.000001, 45))
lower_TRIPRA_gh <- rep(0.000001, 46)
out_TRIPRA_gh <- optim(ini_TRIPRA_gh, f_TRIPRA_gh_optim, lower = lower_TRIPRA_gh, method = 'L-BFGS-B', hessian = T)

ini_VERBOF_gh <- c(1, rep(0.000001, 45))
lower_VERBOF_gh <- rep(0.000001, 46)
out_VERBOF_gh <- optim(ini_VERBOF_gh, f_VERBOF_gh_optim, lower = lower_VERBOF_gh, method = 'L-BFGS-B', hessian = T)

ini_VERPER_gh <- c(1, rep(0.000001, 45))
lower_VERPER_gh <- rep(0.000001, 46)
out_VERPER_gh <- optim(ini_VERPER_gh, f_VERPER_gh_optim, lower = lower_VERPER_gh, method = 'L-BFGS-B', hessian = T)


### extracting and saving the interesting parameters (lambda, alpha and gamma):

lambda_gh <- c(out_ACHMIL_gh$par[1], out_ANTODO_gh$par[1], out_ARRELA_gh$par[1], out_BROERE_gh$par[1], out_CENJAC_gh$par[1], out_CONARV_gh$par[1],
               out_CREPIS_gh$par[1], out_DACGLO_gh$par[1], out_DAUCAR_gh$par[1], out_ELYREP_gh$par[1], out_ERYNGE_gh$par[1], out_FESARU_gh$par[1],
               out_FESRUB_gh$par[1], out_GALVER_gh$par[1], out_GERDIS_gh$par[1], out_GERROT_gh$par[1], out_LEUVUL_gh$par[1], out_LOLPER_gh$par[1],
               out_LOTCOR_gh$par[1], out_MEDARA_gh$par[1], out_ONOREP_gh$par[1], out_PICECH_gh$par[1], out_PICHIE_gh$par[1], out_PLALAN_gh$par[1],
               out_POAANG_gh$par[1], out_POAPRA_gh$par[1], out_POATRI_gh$par[1], out_RANACR_gh$par[1], out_RUMACE_gh$par[1], out_SALPRA_gh$par[1],
               out_SONCHU_gh$par[1], out_TAROFF_gh$par[1], out_TRIFLA_gh$par[1], out_TRIPRA_gh$par[1], out_VERBOF_gh$par[1], out_VERPER_gh$par[1])

alpha_gh <- rbind(out_ACHMIL_gh$par[2:37], out_ANTODO_gh$par[2:37], out_ARRELA_gh$par[2:37], out_BROERE_gh$par[2:37], out_CENJAC_gh$par[2:37], out_CONARV_gh$par[2:37],
                  out_CREPIS_gh$par[2:37], out_DACGLO_gh$par[2:37], out_DAUCAR_gh$par[2:37], out_ELYREP_gh$par[2:37], out_ERYNGE_gh$par[2:37], out_FESARU_gh$par[2:37],
                  out_FESRUB_gh$par[2:37], out_GALVER_gh$par[2:37], out_GERDIS_gh$par[2:37], out_GERROT_gh$par[2:37], out_LEUVUL_gh$par[2:37], out_LOLPER_gh$par[2:37],
                  out_LOTCOR_gh$par[2:37], out_MEDARA_gh$par[2:37], out_ONOREP_gh$par[2:37], out_PICECH_gh$par[2:37], out_PICHIE_gh$par[2:37], out_PLALAN_gh$par[2:37],
                  out_POAANG_gh$par[2:37], out_POAPRA_gh$par[2:37], out_POATRI_gh$par[2:37], out_RANACR_gh$par[2:37], out_RUMACE_gh$par[2:37], out_SALPRA_gh$par[2:37],
                  out_SONCHU_gh$par[2:37], out_TAROFF_gh$par[2:37], out_TRIFLA_gh$par[2:37], out_TRIPRA_gh$par[2:37], out_VERBOF_gh$par[2:37], out_VERPER_gh$par[2:37])

gamma_gh <- rbind(out_ACHMIL_gh$par[41:46], out_ANTODO_gh$par[41:46], out_ARRELA_gh$par[41:46], out_BROERE_gh$par[41:46], out_CENJAC_gh$par[41:46], out_CONARV_gh$par[41:46],
                  out_CREPIS_gh$par[41:46], out_DACGLO_gh$par[41:46], out_DAUCAR_gh$par[41:46], out_ELYREP_gh$par[41:46], out_ERYNGE_gh$par[41:46], out_FESARU_gh$par[41:46],
                  out_FESRUB_gh$par[41:46], out_GALVER_gh$par[41:46], out_GERDIS_gh$par[41:46], out_GERROT_gh$par[41:46], out_LEUVUL_gh$par[41:46], out_LOLPER_gh$par[41:46],
                  out_LOTCOR_gh$par[41:46], out_MEDARA_gh$par[41:46], out_ONOREP_gh$par[41:46], out_PICECH_gh$par[41:46], out_PICHIE_gh$par[41:46], out_PLALAN_gh$par[41:46],
                  out_POAANG_gh$par[41:46], out_POAPRA_gh$par[41:46], out_POATRI_gh$par[41:46], out_RANACR_gh$par[41:46], out_RUMACE_gh$par[41:46], out_SALPRA_gh$par[41:46],
                  out_SONCHU_gh$par[41:46], out_TAROFF_gh$par[41:46], out_TRIFLA_gh$par[41:46], out_TRIPRA_gh$par[41:46], out_VERBOF_gh$par[41:46], out_VERPER_gh$par[41:46])


### naming the rows and columns of the objects lambda, alpha and gamma

names(lambda_gh) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                      "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                      "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")

rownames(alpha_gh) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                        "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                        "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")

colnames(alpha_gh) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                        "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                        "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")

rownames(gamma_gh) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                        "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                        "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")

colnames(gamma_gh) <- c("Cb", "Cd", "Ci", "Ee", "Pg", "Pp")


### removing ANTODO, GERDIS, TRIFLA and VERPER (WARNING: depending on the date selected!)

#lambda_gh <- c(lambda_gh[1], lambda_gh[3:14], lambda_gh[16:32], lambda_gh[34:35])
#alpha_gh <- alpha_gh[c(-2, -15, -33, -36),]
#alpha_gh <- alpha_gh[, c(-2, -15, -33, -36)]
#gamma_gh <- gamma_gh[c(-2, -15, -33, -36),]


### save results:

write.table(lambda_gh, file = "Results/lambda_gh_2.txt", sep = "\t", row.names = TRUE)
write.table(alpha_gh, file = "Results/alpha_gh_2.txt", sep = "\t", row.names = FALSE)
write.table(gamma_gh, file = "Results/alpha_gh_2.txt", sep = "\t", row.names = FALSE)


### clean environment:

#rm(list=ls())

