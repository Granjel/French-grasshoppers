############################################
## Optim analyses --- French grasshoppers ##
############################################
## Rodrigo R. Granjel ## Nov. 2018 #########
############################################

f_link_R <- function(theta, Y, X_plant, X_gh){
  lambda <- theta[1] #max cover without competition
  alpha <- t(t(theta[2:37])) #as many alphas as plant species
  gamma <- matrix(0, nrow = 36,ncol = 1) #grasshopper effects on plants (as many rows as grasshopper species)
  gamma[1:6] <- theta[38:43] #as many thetas as grasshopper species
  log_Y_fit <- log(lambda) - log(1 + X_plant %*% alpha) + log(1 + X_gh %*% gamma) #competition function #X_plant = competitive effects between plants
  #X_pol <- X_gh (nº grasshoppers)
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}
#1 function per plant species, varying the number of parameters of gamma depending on the plant-grasshopper interaction



##################################################################
# OPTIM WITHOUT GRASSHOPPERS:


# 1. Loading the dataset:

d <- read.table("Data_Fg/FG.txt", header = TRUE, sep = "\t")


# 2. Subsample for different dates

d2 <- d[d$time == "1",] #June 2012 ### changing this command changes everything !!!
summary(d2$Focal) #zero ANTODO, GERDIS, TRIFLA and VERPER


# 3. Predictive function depending on the number of plant species

f_competition <- function(theta, Y, X){
  lambda <- theta[1] #max cover without competition
  alpha <- t(t(theta[2:40])) #as many alphas as plant species + 3 (legumes, grasses, other)
  log_Y_fit <- log(lambda) - log(1 + X %*% alpha) #predictive model
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}


# one dataset for each focal species

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


# X and Y matrices for each plant species:

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


# functions for optim, species-specific
f_ACHMIL <- function(theta){f_competition(theta, Y_ACHMIL, X_ACHMIL)}
f_ANTODO <- function(theta){f_competition(theta, Y_ANTODO, X_ANTODO)}
f_ARRELA <- function(theta){f_competition(theta, Y_ARRELA, X_ARRELA)}
f_BROERE <- function(theta){f_competition(theta, Y_BROERE, X_BROERE)}
f_CENJAC <- function(theta){f_competition(theta, Y_CENJAC, X_CENJAC)}
f_CONARV <- function(theta){f_competition(theta, Y_CONARV, X_CONARV)}
f_CREPIS <- function(theta){f_competition(theta, Y_CREPIS, X_CREPIS)}
f_DACGLO <- function(theta){f_competition(theta, Y_DACGLO, X_DACGLO)}
f_DAUCAR <- function(theta){f_competition(theta, Y_DAUCAR, X_DAUCAR)}
f_ELYREP <- function(theta){f_competition(theta, Y_ELYREP, X_ELYREP)}
f_ERYNGE <- function(theta){f_competition(theta, Y_ERYNGE, X_ERYNGE)}
f_FESARU <- function(theta){f_competition(theta, Y_FESARU, X_FESARU)}
f_FESRUB <- function(theta){f_competition(theta, Y_FESRUB, X_FESRUB)}
f_GALVER <- function(theta){f_competition(theta, Y_GALVER, X_GALVER)}
f_GERDIS <- function(theta){f_competition(theta, Y_GERDIS, X_GERDIS)}
f_GERROT <- function(theta){f_competition(theta, Y_GERROT, X_GERROT)}
f_LEUVUL <- function(theta){f_competition(theta, Y_LEUVUL, X_LEUVUL)}
f_LOLPER <- function(theta){f_competition(theta, Y_LOLPER, X_LOLPER)}
f_LOTCOR <- function(theta){f_competition(theta, Y_LOTCOR, X_LOTCOR)}
f_MEDARA <- function(theta){f_competition(theta, Y_MEDARA, X_MEDARA)}
f_ONOREP <- function(theta){f_competition(theta, Y_ONOREP, X_ONOREP)}
f_PICECH <- function(theta){f_competition(theta, Y_PICECH, X_PICECH)}
f_PICHIE <- function(theta){f_competition(theta, Y_PICHIE, X_PICHIE)}
f_PLALAN <- function(theta){f_competition(theta, Y_PLALAN, X_PLALAN)}
f_POAANG <- function(theta){f_competition(theta, Y_POAANG, X_POAANG)}
f_POAPRA <- function(theta){f_competition(theta, Y_POAPRA, X_POAPRA)}
f_POATRI <- function(theta){f_competition(theta, Y_POATRI, X_POATRI)}
f_RANACR <- function(theta){f_competition(theta, Y_RANACR, X_RANACR)}
f_RUMACE <- function(theta){f_competition(theta, Y_RUMACE, X_RUMACE)}
f_SALPRA <- function(theta){f_competition(theta, Y_SALPRA, X_SALPRA)}
f_SONCHU <- function(theta){f_competition(theta, Y_SONCHU, X_SONCHU)}
f_TAROFF <- function(theta){f_competition(theta, Y_TAROFF, X_TAROFF)}
f_TRIFLA <- function(theta){f_competition(theta, Y_TRIFLA, X_TRIFLA)}
f_TRIPRA <- function(theta){f_competition(theta, Y_TRIPRA, X_TRIPRA)}
f_VERBOF <- function(theta){f_competition(theta, Y_VERBOF, X_VERBOF)}
f_VERPER <- function(theta){f_competition(theta, Y_VERPER, X_VERPER)}


# optim wrap for the diff. species ---:

ini <- rep(1, 40)
low <- rep(0, 40)

out_ACHMIL <- optim(ini, f_ACHMIL, lower = low, method = 'L-BFGS-B', hessian = T)
out_ANTODO <- optim(ini, f_ANTODO, lower = low, method = 'L-BFGS-B', hessian = T)
out_ARRELA <- optim(ini, f_ARRELA, lower = low, method = 'L-BFGS-B', hessian = T)
out_BROERE <- optim(ini, f_BROERE, lower = low, method = 'L-BFGS-B', hessian = T)
out_CENJAC <- optim(ini, f_CENJAC, lower = low, method = 'L-BFGS-B', hessian = T)
out_CONARV <- optim(ini, f_CONARV, lower = low, method = 'L-BFGS-B', hessian = T)
out_CREPIS <- optim(ini, f_CREPIS, lower = low, method = 'L-BFGS-B', hessian = T)
out_DACGLO <- optim(ini, f_DACGLO, lower = low, method = 'L-BFGS-B', hessian = T)
out_DAUCAR <- optim(ini, f_DAUCAR, lower = low, method = 'L-BFGS-B', hessian = T)
out_ELYREP <- optim(ini, f_ELYREP, lower = low, method = 'L-BFGS-B', hessian = T)
out_ERYNGE <- optim(ini, f_ERYNGE, lower = low, method = 'L-BFGS-B', hessian = T)
out_FESARU <- optim(ini, f_FESARU, lower = low, method = 'L-BFGS-B', hessian = T)
out_FESRUB <- optim(ini, f_FESRUB, lower = low, method = 'L-BFGS-B', hessian = T)
out_GALVER <- optim(ini, f_GALVER, lower = low, method = 'L-BFGS-B', hessian = T)
out_GERDIS <- optim(ini, f_GERDIS, lower = low, method = 'L-BFGS-B', hessian = T)
out_GERROT <- optim(ini, f_GERROT, lower = low, method = 'L-BFGS-B', hessian = T)
out_LEUVUL <- optim(ini, f_LEUVUL, lower = low, method = 'L-BFGS-B', hessian = T)
out_LOLPER <- optim(ini, f_LOLPER, lower = low, method = 'L-BFGS-B', hessian = T)
out_LOTCOR <- optim(ini, f_LOTCOR, lower = low, method = 'L-BFGS-B', hessian = T)
out_MEDARA <- optim(ini, f_MEDARA, lower = low, method = 'L-BFGS-B', hessian = T)
out_ONOREP <- optim(ini, f_ONOREP, lower = low, method = 'L-BFGS-B', hessian = T)
out_PICECH <- optim(ini, f_PICECH, lower = low, method = 'L-BFGS-B', hessian = T)
out_PICHIE <- optim(ini, f_PICHIE, lower = low, method = 'L-BFGS-B', hessian = T)
out_PLALAN <- optim(ini, f_PLALAN, lower = low, method = 'L-BFGS-B', hessian = T)
out_POAANG <- optim(ini, f_POAANG, lower = low, method = 'L-BFGS-B', hessian = T)
out_POAPRA <- optim(ini, f_POAPRA, lower = low, method = 'L-BFGS-B', hessian = T)
out_POATRI <- optim(ini, f_POATRI, lower = low, method = 'L-BFGS-B', hessian = T)
out_RANACR <- optim(ini, f_RANACR, lower = low, method = 'L-BFGS-B', hessian = T)
out_RUMACE <- optim(ini, f_RUMACE, lower = low, method = 'L-BFGS-B', hessian = T)
out_SALPRA <- optim(ini, f_SALPRA, lower = low, method = 'L-BFGS-B', hessian = T)
out_SONCHU <- optim(ini, f_SONCHU, lower = low, method = 'L-BFGS-B', hessian = T)
out_TAROFF <- optim(ini, f_TAROFF, lower = low, method = 'L-BFGS-B', hessian = T)
out_TRIFLA <- optim(ini, f_TRIFLA, lower = low, method = 'L-BFGS-B', hessian = T)
out_TRIPRA <- optim(ini, f_TRIPRA, lower = low, method = 'L-BFGS-B', hessian = T)
out_VERBOF <- optim(ini, f_VERBOF, lower = low, method = 'L-BFGS-B', hessian = T)
out_VERPER <- optim(ini, f_VERPER, lower = low, method = 'L-BFGS-B', hessian = T)



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

rownames(alpha) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                     "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                     "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")
colnames(alpha) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                     "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                     "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")
names(lambda) <- c("ACHMIL", "ANTODO", "ARRELA", "BROERE", "CENJAC", "CONARV", "CREPIS", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FESARU",
                   "FESRUB", "GALVER", "GERDIS", "GERROT", "LEUVUL", "LOLPER", "LOTCOR", "MEDARA", "ONOREP", "PICECH", "PICHIE", "PLALAN",
                   "POAANG", "POAPRA", "POATRI", "RANACR", "RUMACE", "SALPRA", "SONCHU", "TAROFF", "TRIFLA", "TRIPRA", "VERBOF", "VERPER")

#removing ANTODO, GERDIS, TRIFLA and VERPER
alpha <- alpha[c(-2, -15, -33, -36),]
alpha <- alpha[, c(-2, -15, -33, -36)]
lambda <- c(lambda[1], lambda[3:14], lambda[16:32], lambda[34:35])


######################################################################

d2 <- d[(d$treatment == 'no_link' | d$treatment == 'NO_Pol'),] #remove

levels(d2$focal_plant)
str(d2)


d_T <- d2[d2$focal_plant == 'T',]
d_R <- d2[d2$focal_plant == 'R',]
d_H <- d2[d2$focal_plant == 'H',]

X_T_plant <- as.matrix(d_T[4:6]) #competition
X_T_pol <- as.matrix(d_T[7:9]) #7:9 == depends on the number of links for every species (suppl. mat. Gross)
Y_T <- d_T$Y #cover
#repeat for each spp

X_R_plant <- as.matrix(d_R[4:6])
X_R_pol <- as.matrix(d_R[7:9])
Y_R <- d_R$Y

X_H_plant <- as.matrix(d_H[4:6])
X_H_pol <- as.matrix(d_H[7:9])
Y_H <- d_H$Y


f_T <- function(theta){f_no_link_T(theta,Y_T,X_T_plant,X_T_pol)}
f_R <- function(theta){f_no_link_R(theta,Y_R,X_R_plant,X_R_pol)}
f_H <- function(theta){f_no_link_H(theta,Y_H,X_H_plant,X_H_pol)}

out_T <- optim(c(1,1,1,1,1),f_T,lower = c(0,0,0,0,-0.08),method = 'L-BFGS-B',hessian = T)
out_T
out_R <- optim(c(1,1,1,1,1,1),f_R,lower = c(0,0,0,0,0,0),method = 'L-BFGS-B',hessian = T)
out_R
out_H <- optim(c(1,1,1,1,1,1),f_H,lower = c(0,0,0,0,-0.01,-0.03),method = 'L-BFGS-B',hessian = T)
out_H


lambda_no_link <- c(out_T$par[1],out_R$par[1],out_H$par[1])
alpha_no_link <- rbind(out_T$par[2:4],out_R$par[2:4],out_H$par[2:4])
gamma_no_link <- matrix(0,nrow = 3,ncol = 3)
gamma_no_link[c(4,2,8,3,6)] <- c(out_T$par[5],out_R$par[5:6],out_H$par[5:6])

rownames(alpha_no_link) <- c('T','R','H')
colnames(alpha_no_link) <- c('T','R','H')
rownames(gamma_no_link) <- c('T','R','H')
colnames(gamma_no_link) <- c('O','B','F')
names(lambda_no_link) <- c('T','R','H')



#####################################################################


d2 <- d[(d$treatment == 'link' | d$treatment == 'NO_Pol'),]
d2 <- d2[d2$Y > 0,]

levels(d2$focal_plant)
str(d2)


d_T <- d2[d2$focal_plant == 'T',]
d_R <- d2[d2$focal_plant == 'R',]
d_H <- d2[d2$focal_plant == 'H',]

X_T_plant <- as.matrix(d_T[4:6])
X_T_pol <- as.matrix(d_T[7:9])
Y_T <- d_T$Y

X_R_plant <- as.matrix(d_R[4:6])
X_R_pol <- as.matrix(d_R[7:9])
Y_R <- d_R$Y

X_H_plant <- as.matrix(d_H[4:6])
X_H_pol <- as.matrix(d_H[7:9])
Y_H <- d_H$Y


f_T <- function(theta){f_link_T(theta,Y_T,X_T_plant,X_T_pol)}
f_R <- function(theta){f_link_R(theta,Y_R,X_R_plant,X_R_pol)}
f_H <- function(theta){f_link_H(theta,Y_H,X_H_plant,X_H_pol)}

out_T <- optim(c(1,1,1,1,1),f_T,lower = c(0,0,0,0,-0.05),method = 'L-BFGS-B',hessian = T)
out_T
out_R <- optim(c(1,1,1,1,1,1,1),f_R,lower = c(0,0,0,0,0,-0.01,-0.001),method = 'L-BFGS-B',hessian = T)
out_R
out_H <- optim(c(1,1,1,1,1,1),f_H,lower = c(0,0,0,0,-0.01,-0.03),method = 'L-BFGS-B',hessian = T)
out_H

lambda_link <- c(out_T$par[1],out_R$par[1],out_H$par[1])
alpha_link <- rbind(out_T$par[2:4],out_R$par[2:4],out_H$par[2:4])
gamma_link <- matrix(0,nrow = 3,ncol = 3)
gamma_link[c(4,2,5,8,3,6)] <- c(out_T$par[5],out_R$par[5:7],out_H$par[5:6])

rownames(alpha_link) <- c('T','R','H')
colnames(alpha_link) <- c('T','R','H')
rownames(gamma_link) <- c('T','R','H')
colnames(gamma_link) <- c('O','B','F')
names(lambda_link) <- c('T','R','H')

lambda_link
alpha_link
gamma_link

