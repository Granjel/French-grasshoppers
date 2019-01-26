#load bipartite for network metrics: conectance, nestedness and moudlarity
library(bipartite)
source("Code/sigma_to_gamma.R")


#1. field conditions

#load bipartite matrix
sigma <- read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1) #negative interaction: herbivory
d <- sigma

#load dataset of outcomes
field <- read.table("Results/structural_coex_field_conditions.txt", header = TRUE, sep = "\t")
f <- field

#define different features of the bipartite network
# triplets: 3 plant species
tri <- apply(t(combn(rownames(d), 3)), 1, paste, collapse = "_")
loops <- as.numeric(nrow(f))/1540
tri <- rep(tri, loops)
# scenario considered: f conditions = 1
sce <- rep(1, nrow(f)) #rep = number of rows
# variability of the interaction strength = 1 (untouched)
vis <- rep(1, nrow(f))
# conectance (applies function conect())
con <- rep(conect(d), nrow(f))
# nestedness (applies function nest() within bipartite)
nes <- rep(as.numeric(nested(d, method = "weighted NODF", rescale = FALSE, normalised = FALSE)), nrow(f))
# modularity (applies function computeModules() within bipartite)
mod <- rep(as.numeric(computeModules(d, method = "Beckett")@likelihood), nrow(f))

#define a dataframe for the bipartite network properties and else
df <- data.frame("triplet" = tri, "scenario" = sce, "VIS" = vis, "conectance" = con, "nestedness" = nes, "modularity" = mod)

#put together factors and data
DATA <- cbind(df, field)

#rm and end of this scenario
rm(sigma, field, df)



#1. field conditions

#load bipartite matrix: there is no need, it is a matrix of zeroes
#d <- NULL (actually we need to keep it from scenario 1 for the plant names)

#load dataset of outcomes
no_herb <- read.table("Results/structural_coex_no_herbivory.txt", header = TRUE, sep = "\t")
f <- no_herb

#define different features of the bipartite network
# triplets: 3 plant species
tri <- apply(t(combn(rownames(d), 3)), 1, paste, collapse = "_")
loops <- as.numeric(nrow(f))/1540
tri <- rep(tri, loops)
# scenario considered: divory = 2
sce <- rep(2, nrow(f)) #rep = number of rows
# variability of the interaction strength = 0 (none)
vis <- rep(0, nrow(f))
# conectance (applies function conect())
con <- rep(0, nrow(f))
# nestedness (applies function nest() within bipartite)
nes <- rep(0, nrow(f))
# modularity (applies function computeModules() within bipartite)
mod <- rep(0, nrow(f))

#define a dataframe for the bipartite network properties and else
df <- data.frame("triplet" = tri, "scenario" = sce, "VIS" = vis, "conectance" = con, "nestedness" = nes, "modularity" = mod)

#put together factors and data
sc2 <- cbind(df, f)

#add to the previous data frame
DATA <- rbind(DATA, sc2)

#rm and end of this scenario
rm(no_herb, df, sc2)



#3. VIS

#load bipartite matrix
sigma <- read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1) #negative interaction: herbivory
d <- sigma

#the following code is not really reproducible, but I can't do it properly right now. CHANGE.
strength_seq <- c((1/3), 3, 9, 27) #the four strength values selected. Not reproducible
sigma_vis <- list() #save as a list and fill it
for (i in 1:length(strength_seq)){
  sigma_vis[[i]] <- sigma * strength_seq[i]
}

#load dataset of outcomes
vis_res <- read.table("Results/structural_coex_VIS_2.txt", header = TRUE, sep = "\t")
f <- vis_res

#define different features of the bipartite network
# triplets: 3 plant species
tri <- apply(t(combn(rownames(d), 3)), 1, paste, collapse = "_")
loops <- 4 #not reproducible
tri <- rep(tri, loops)

# scenario considered: VIS = 3
sce <- rep(3, nrow(f)) #rep = number of rows
# variability of the interaction strength = 1 (untouched)
vis <- c(rep(strength_seq[1], nrow(f)/4), rep(strength_seq[2], nrow(f)/4),
         rep(strength_seq[3], nrow(f)/4), rep(strength_seq[4], nrow(f)/4))

d <- sigma_vis #define 'd' again...

# nestedness (applies function nest() within bipartite)
con <- NULL
nes <- NULL
mod <- NULL
for (i in 1:length(d)){
  # conectance (applies function conect())
  con <- c(con, rep(conect(d[[i]]), 1540))
  # nestedness (applies function nest() within bipartite)
  nes <- c(nes, rep(as.numeric(nested(d[[i]], method = "weighted NODF", rescale = FALSE, normalised = FALSE)), 1540))
  # modularity (applies function computeModules() within bipartite)
  mod <- c(mod, rep(as.numeric(computeModules(d[[i]], method = "Beckett")@likelihood), 1540))
}

#define a dataframe for the bipartite network properties and else
df <- data.frame("triplet" = tri, "scenario" = sce, "VIS" = vis,
                 "conectance" = con, "nestedness" = nes, "modularity" = mod)

#put together factors and data
sc3 <- cbind(df, f)

#add to the previous data frame
DATA <- rbind(DATA, sc3)

#rm and end of this scenario
rm(sigma, vis_res, sigma_vis, df, sc3, i, strength_seq)



#4. Metrics

#load bipartite matrix
matrices <- read.table("Results/metrics_4rep_matrices.txt", header = TRUE, sep = "\t")
matrices <- recover_sigma(matrices, 22, 6)
d <- matrices

#load dataset of outcomes
vis_res <- read.table("Results/structural_coex_VIS_2.txt", header = TRUE, sep = "\t")
f <- vis_res

#define different features of the bipartite network
# triplets: 3 plant species
tri <- apply(t(combn(rownames(d), 3)), 1, paste, collapse = "_")
loops <- 4 #not reproducible
tri <- rep(tri, loops)

# scenario considered: VIS = 3
sce <- rep(3, nrow(f)) #rep = number of rows
# variability of the interaction strength = 1 (untouched)
vis <- c(rep(strength_seq[1], nrow(f)/4), rep(strength_seq[2], nrow(f)/4),
         rep(strength_seq[3], nrow(f)/4), rep(strength_seq[4], nrow(f)/4))

d <- sigma_vis #define 'd' again...

# nestedness (applies function nest() within bipartite)
con <- NULL
nes <- NULL
mod <- NULL
for (i in 1:length(d)){
  # conectance (applies function conect())
  con <- c(con, rep(conect(d[[i]]), 1540))
  # nestedness (applies function nest() within bipartite)
  nes <- c(nes, rep(as.numeric(nested(d[[i]], method = "weighted NODF", rescale = FALSE, normalised = FALSE)), 1540))
  # modularity (applies function computeModules() within bipartite)
  mod <- c(mod, rep(as.numeric(computeModules(d[[i]], method = "Beckett")@likelihood), 1540))
}

#define a dataframe for the bipartite network properties and else
df <- data.frame("triplet" = tri, "scenario" = sce, "VIS" = vis,
                 "conectance" = con, "nestedness" = nes, "modularity" = mod)

#put together factors and data
sc3 <- cbind(df, f)

#add to the previous data frame
DATA <- rbind(DATA, sc3)

#rm and end of this scenario
rm(sigma, vis_res, sigma_vis, df, sc3, i, strength_seq)


#metrics_4rep_matrices.txt!!!!!!!!!!




metrics1 <- read.table("Results/structural_coex_metrics.txt", header = TRUE, sep = "\t")
metrics2 <- read.table("Results/structural_coex_metrics_sobremesa.txt", header = TRUE, sep = "\t")
metrics3 <- read.table("Results/structural_coex_metrics_Teresa.txt", header = TRUE, sep = "\t")
metrics <- rbind(metrics1, metrics2, metrics3)
rm(metrics1, metrics2, metrics3)







