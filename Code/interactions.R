#script to calculate the coexistence outputs for the interactions between VIM and network metrics

library(bipartite)

#nested(sigma)
#mean(as.numeric(lapply(deletions25, nested)))
#sd(as.numeric(lapply(deletions25, nested)))
#mean(as.numeric(lapply(deletions50, nested)))
#sd(as.numeric(lapply(deletions50, nested)))
#mean(as.numeric(lapply(deletions75, nested)))
#sd(as.numeric(lapply(deletions75, nested)))

a <- nested(sigma, method = "weighted NODF", rescale = FALSE, normalised = FALSE)

n <- NULL

for (i in 1:reps){
  n <- c(n, nested(deletions25[[i]], method = "weighted NODF", rescale = FALSE, normalised = FALSE),
         nested(deletions50[[i]], method = "weighted NODF", rescale = FALSE, normalised = FALSE),
         nested(deletions75[[i]], method = "weighted NODF", rescale = FALSE, normalised = FALSE))
}

deletions_percentage <- c(1, rep(25, 20), rep(50, 20), rep(75, 20))
nestedness <- c(a, n)

df <- data.frame(deletions_percentage, nestedness)

plot(df)

mean(df$nestedness[df$deletions_percentage == 25])
mean(df$nestedness[df$deletions_percentage == 50])
mean(df$nestedness[df$deletions_percentage == 75])




networklevel(t(sigma))

wine(sigma, nreps = 10)$'win'
plot(wine(sigma, nreps = 1000))

wine(deletions25[[4]], nreps = 10)$'win'
wine(deletions25[[4]], nreps = 10)$wine




V.ratio(sigma) #Calculates the variance-ratio as suggested by Schluter (1984)
computeModules(sigma)
nested
degreedistr(sigma)















































