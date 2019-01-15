rm(list=ls())

f_competition <- function(theta, Y ,X){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  log_Y_fit <- log(lambda) - log( 1 + X %*% alpha)
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}


f_link_T <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[2] <- theta[5]
  log_Y_fit <- log(lambda) - log( 1 + X_plant %*% alpha) + log( 1 + X_pol %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_link_R <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[1:3] <- theta[5:7]
  log_Y_fit <- log(lambda) - log( 1 + X_plant %*% alpha) + log( 1 + X_pol %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_link_H <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[1:2] <- theta[5:6]
  log_Y_fit <- log(lambda) - log( 1 + X_plant %*% alpha) + log( 1 + X_pol %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}


f_no_link_T <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[2] <- theta[5]
  log_Y_fit <- log(lambda) - log( 1 + X_plant %*% alpha) + log( 1 + X_pol %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_no_link_R <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[c(1,3)] <- theta[5:6]
  log_Y_fit <- log(lambda) - log( 1 + X_plant %*% alpha) + log( 1 + X_pol %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_no_link_H <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[1:2] <- theta[5:6]
  log_Y_fit <- log(lambda) - log( 1 + X_plant %*% alpha) + log( 1 + X_pol %*% gamma) 
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}





########################################################

d <- read.csv(file='/Users/oscargodoy/Desktop/Rcode for coexistence/EXPLORA/data_for_rudolf.csv',header = T)
str(d)
d$Y <- (d$seeds  - (1-d$g) * d$s) / d$g
#d$Y <- d$seeds
levels(d$treatment)

########################################################

d2 <- d[d$treatment == 'NO_Pol',]

levels(d2$focal_plant)
str(d2)


d_T <- d2[d2$focal_plant == 'T',]
d_R <- d2[d2$focal_plant == 'R',]
d_H <- d2[d2$focal_plant == 'H',]

X_T <- as.matrix(d_T[4:6])
Y_T <- d_T$Y

X_R <- as.matrix(d_R[4:6])
Y_R <- d_R$Y

X_H <- as.matrix(d_H[4:6])
Y_H <- d_H$Y


f_T <- function(theta){f_competition(theta,Y_T,X_T)}
f_R <- function(theta){f_competition(theta,Y_R,X_R)}
f_H <- function(theta){f_competition(theta,Y_H,X_H)}

out_T <- optim(c(1,1,1,1),f_T,lower = c(0,0,0,0),method = 'L-BFGS-B',hessian = T)
out_T
out_R <- optim(c(1,1,1,1),f_R,lower = c(0,0,0,0),method = 'L-BFGS-B',hessian = T)
out_R
out_H <- optim(c(1,1,1,1),f_H,lower = c(0,0,0,0),method = 'L-BFGS-B',hessian = T)
out_H


lambda <- c(out_T$par[1],out_R$par[1],out_H$par[1])
alpha <- rbind(out_T$par[2:4],out_R$par[2:4],out_H$par[2:4])

rownames(alpha) <- c('T','R','H')
colnames(alpha) <- c('T','R','H')
names(lambda) <- c('T','R','H')

######################################################################

d2 <- d[(d$treatment == 'no_link' | d$treatment == 'NO_Pol'),]

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

