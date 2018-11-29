install.packages("glinternet")
require(glinternet)
library(glinternet)


### load dataset:

d <- read.table("Data_Fg/FG.txt", header = TRUE, sep = "\t")


### subsample for different dates

d2 <- d[d$time == "2",] #June 2012 ### changing this command changes everything !!!
summary(d2$Focal)


### one dataset for each focal species

d_CENJAC_gh <- d2[d2$Focal == "CENJAC",]


### X_plant, X_gh and Y matrices for each plant species:

X_CENJAC_plant <- as.matrix(d_CENJAC_gh[, seq(15, 91, by = 2)]) #competition
#X_CENJAC_plant <- X_CENJAC_plant[c(-2, -15, -20, -27, -33, -36),]
X_CENJAC_gh <- as.matrix(d_CENJAC_gh[6:11]) #depends on the number of links for every species (suppl. mat. Gross)
Y_CENJAC_gh <- d_CENJAC_gh$Cover #cover
X_CENJAC <- cbind(X_CENJAC_plant, X_CENJAC_gh)


MODEL <- glinternet(X_CENJAC, Y_CENJAC_gh, numLevels = rep(1, 45), lambda = NULL, nLambda = 500, lambdaMinRatio = 0.01,
                    interactionCandidates = c(40:45), screenLimit = NULL, numToFind = NULL, family = "gaussian", tol = 1e-05,
                    maxIter = 5000, verbose = TRUE, numCores = 1)


bet <- NULL
for(i in 1:500){
  bet <- c(bet, length(MODEL$betahat[[i]]))
}
plot(MODEL$lambda, bet, ylim = c(0,500))

which(MODEL$lambda == min(MODEL$lambda))


coef(MODEL, lambdaIndex = 500)



length(MODEL$betahat[[500]])

coef(MODEL, lambdaIndex = 400)

MODEL$lambda[2]


str(MODEL)

lenght(MODEL$betahat[250])

