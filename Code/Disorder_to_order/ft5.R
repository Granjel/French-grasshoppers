rm(list=ls())

f_competition <- function(theta, Y ,X){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  log_Y_fit <- lambda - X %*% alpha
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}


f_link_T <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[2] <- theta[5]
  log_Y_fit <- lambda - X_plant %*% alpha + X_pol %*% gamma
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_link_R <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[1:3] <- theta[5:7]
  log_Y_fit <- lambda - X_plant %*% alpha + X_pol %*% gamma
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_link_H <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[1:2] <- theta[5:6]
  log_Y_fit <- lambda - X_plant %*% alpha + X_pol %*% gamma
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}


f_no_link_T <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[2] <- theta[5]
  log_Y_fit <- lambda - X_plant %*% alpha + X_pol %*% gamma
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_no_link_R <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[c(1,3)] <- theta[5:6]
  log_Y_fit <- lambda - X_plant %*% alpha + X_pol %*% gamma
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}

f_no_link_H <- function(theta, Y ,X_plant,X_pol){
  lambda <- theta[1] 
  alpha <- t(t(theta[2:4]))
  gamma <- matrix(0,nrow = 3,ncol = 1)
  gamma[1:2] <- theta[5:6]
  log_Y_fit <- lambda - X_plant %*% alpha + X_pol %*% gamma
  SS <- sum((log(Y)-log_Y_fit)^2)
  return(SS)
}





########################################################

d <- read.csv(file='data_for_rudolf.csv',header = T)
str(d)
#d$Y <- (d$seeds  - (1-d$g) * d$s) / d$g
d$Y <- d$seeds
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

out_T <- lm(log(Y_T) ~ X_T)
summary(out_T)
qqnorm(residuals(out_T)); qqline(residuals(out_T))
out_R <- lm(log(Y_R) ~ X_R)
summary(out_R)
qqnorm(residuals(out_R)); qqline(residuals(out_R))
out_H <- lm(log(Y_H) ~ X_H)
summary(out_H)
qqnorm(residuals(out_H)); qqline(residuals(out_H))

r <- c(out_T$coefficients[1],out_R$coefficients[1],out_H$coefficients[1])
alpha <- -rbind(out_T$coefficients[2:4],out_R$coefficients[2:4],out_H$coefficients[2:4])

rownames(alpha) <- c('T','R','H')
colnames(alpha) <- c('T','R','H')
names(r) <- c('T','R','H')

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


out_T <- lm(log(Y_T) ~ X_T_plant + X_T_pol[,2])
summary(out_T)
qqnorm(residuals(out_T)); qqline(residuals(out_T))
out_R <- lm(log(Y_R) ~ X_R_plant + X_R_pol[,c(1,3)])
summary(out_R)
qqnorm(residuals(out_R)); qqline(residuals(out_R))
out_H <- lm(log(Y_H) ~ X_H_plant + X_H_pol[,c(1,2)])
summary(out_H)
qqnorm(residuals(out_H)); qqline(residuals(out_H))

r_no_link <- c(out_T$coefficients[1],out_R$coefficients[1],out_H$coefficients[1])
alpha_no_link <- -rbind(out_T$coefficients[2:4],out_R$coefficients[2:4],out_H$coefficients[2:4])
gamma_no_link <- matrix(0,nrow = 3,ncol = 3)
gamma_no_link[c(4,2,8,3,6)] <- c(out_T$coefficients[5],out_R$coefficients[5:6],out_H$coefficients[5:6])

rownames(alpha_no_link) <- c('T','R','H')
colnames(alpha_no_link) <- c('T','R','H')
rownames(gamma_no_link) <- c('T','R','H')
colnames(gamma_no_link) <- c('O','B','F')
names(r_no_link) <- c('T','R','H')



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


out_T <- lm(log(Y_T) ~ X_T_plant + X_T_pol[,2])
summary(out_T)
qqnorm(residuals(out_T)); qqline(residuals(out_T))
out_R <- lm(log(Y_R) ~ X_R_plant + X_R_pol[,c(1,2,3)])
summary(out_R)
qqnorm(residuals(out_R)); qqline(residuals(out_R))
out_H <- lm(log(Y_H) ~ X_H_plant + X_H_pol[,c(1,2)])
summary(out_H)
qqnorm(residuals(out_H)); qqline(residuals(out_H))


r_link <- c(out_T$coefficients[1],out_R$coefficients[1],out_H$coefficients[1])
alpha_link <- -rbind(out_T$coefficients[2:4],out_R$coefficients[2:4],out_H$coefficients[2:4])
gamma_link <- matrix(0,nrow = 3,ncol = 3)
gamma_link[c(4,2,5,8,3,6)] <- c(out_T$coefficients[5],out_R$coefficients[5:7],out_H$coefficients[5:6])

rownames(alpha_link) <- c('T','R','H')
colnames(alpha_link) <- c('T','R','H')
rownames(gamma_link) <- c('T','R','H')
colnames(gamma_link) <- c('O','B','F')
names(r_link) <- c('T','R','H')

#######################

r
r_link
r_no_link

alpha
alpha_link
alpha_no_link

gamma_link
gamma_no_link
