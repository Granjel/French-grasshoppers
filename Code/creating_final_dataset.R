#######################################################################
### Legumes, grasses and other plants - reshaping the dataset, again ##
### Rodrigo R. Granjel ------- 28 June 2018 ###########################
#######################################################################

## load database
French_grasshoppers <- read.table(file = "C:/Users/Granjel RR/Desktop/Nico Gross/French_grasshoppers_converted_step1.txt", header = TRUE, sep = "\t")

## new dataframe to save the target rows
FG <- data.frame()
OTHER <- data.frame() #this is to save the non-focal species, "just in case"

## list of non-focal species, classified by category
legumes <- c("MEDLUP", "TRICAM", "TRIREP", "VICSAT")
grasses <- c("BROMOL", "BROSTE", "POA.SP", "POAPOI")
other <- c("AGREUP", "CARCAR", "CRULAE", "FALVUL", "FRAEXE", "GALAPA", "GALMOL", "HIMHIR", "MALSYL", "MYORAM", "ORCHID", "PLAMAJ", "POTREP", "PRIVUL", "PRUVUL", "RANREP", "RUBFRU", "SENJAC", "TRAPRA")

## there we go -- 1st, removing the rows with non-focal plants
for (i in 1:nrow(French_grasshoppers)){
  if (French_grasshoppers$Focal[i] %in% legumes || French_grasshoppers$Focal[i] %in% grasses ||
      French_grasshoppers$Focal[i] %in% other){
    OTHER <- rbind(OTHER, French_grasshoppers[i,])
  } else {
    FG <- rbind(FG, French_grasshoppers[i,])
  }
  print((i/nrow(French_grasshoppers))*100)
}

## locating the position of the columns with legumes, grasses and other plants
pos_leg <- NULL
pos_gra <- NULL
pos_oth <- NULL

nms <- colnames(French_grasshoppers) #column names

for (i in 1:length(nms)){
  if (nms[i] %in% legumes){
    pos_leg <- c(pos_leg, i)
  }
  if (nms[i] %in% grasses){
    pos_gra <- c(pos_gra, i)
  }
  if (nms[i] %in% other){
    pos_oth <- c(pos_oth, i)
  }
}

## summing the cover and individuals of the non-focal species
leg <- rep(0, nrow(FG))
leg_i <- rep(0, nrow(FG))
gra <- rep(0, nrow(FG))
gra_i <- rep(0, nrow(FG))
oth <- rep(0, nrow(FG))
oth_i <- rep(0, nrow(FG))

for (i in 1:nrow(FG)){
  #legumes
  for (j in pos_leg){
    leg[i] <- leg[i] + FG[i,j]
    leg_i[i] <- leg_i[i] + FG[i,j-1]
  }
  #grasses
  for (w in pos_gra){
    gra[i] <- gra[w] + FG[i,w]
    gra_i[i] <- gra_i[w] + FG[i,w-1]
  }
  #other plants
  for (z in pos_oth){
    oth[i] <- oth[z] + FG[i,z]
    oth_i[i] <- oth_i[z] + FG[i,z-1]
  }
  print((i/nrow(FG))*100)
}

## removing the columns with undesired species from the database
pos <- c(pos_leg, (pos_leg - 1), pos_gra, (pos_gra - 1), pos_oth, (pos_oth - 1))
FG <- FG[-pos]

## saving the new categories into the database and writing the database in a txt file
FG <- cbind(FG, "legumes_i" = leg_i, "legumes" = leg, "grasses_i" = gra_i, "grasses" = gra, "other_i" = oth_i, "other" = oth)
write.table(FG, file = "C:/Users/Granjel RR/Desktop/Nico Gross/FG.txt", sep = "\t")

## clean up
rm(gra, gra_i, grasses, i, j, leg, leg_i, legumes, nms, oth, oth_i, other, pos, pos_gra, pos_leg, pos_oth, w, z)