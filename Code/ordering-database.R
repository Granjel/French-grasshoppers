# Initializing

species <- c("ACHMIL", "AGREUP", "ANTODO", "ARRELA", "BROERE", "BROMOL", "BROSTE", "CARCAR", "CENJAC",
             "CONARV", "CREPIS", "CRULAE", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FALVUL", "FESARU",
             "FESRUB", "FRAEXE", "GALAPA", "GALMOL", "GALVER", "GERDIS", "GERROT", "HIMHIR", "LEUVUL",
             "LOLPER", "LOTCOR", "MALSYL", "MEDARA", "MEDLUP", "MYORAM", "ONOREP", "ORCHID", "PICECH",
             "PICHIE", "PLALAN", "PLAMAJ", "POAANG", "POAPOI", "POAPRA", "POA.SP", "POATRI", "POTREP",
             "PRIVUL", "PRUVUL", "RANACR", "RANREP", "RUBFRU", "RUMACE", "SALPRA", "SENJAC", "SONCHU",
             "TAROFF", "TRAPRA", "TRICAM", "TRIFLA", "TRIPRA", "TRIREP", "VERBOF", "VERPER", "VICSAT")


# JUNE 2012

jun12.raw <- read.table(file = "D:/jun12.txt", header = TRUE, sep = "\t")
jun12.oka <- read.table(file = "D:/skeletonjun12.txt", header = TRUE, sep = "\t")
copyjun12 <- jun12.oka

## filling the skeleton of the new dataframe with data
for (i in 1:nrow(jun12.oka)){
  for (j in 1:nrow(jun12.raw)){
    if (jun12.oka$block[i] == jun12.raw$block[j] && jun12.oka$treatment[i] == jun12.raw$treatment[j] &&
        jun12.oka$datapoint[i] == jun12.raw$datapoint[j]){
      for (z in 1:length(species)){
        if (jun12.raw$species[j] == species[z]){
          copyjun12[i, z+11] <- jun12.raw$cover[j]
        }
      }
    }
  }
}

## writes a txt file that must be carefully opened with Excel (separator = space!)
jun12.oka <- copyjun12
write.table(jun12.oka, file = "jun12.txt", na = "NA", row.names = FALSE)
### remember to change the location of the dataframes so we don't post them on GitHub!

## cleaning behind
rm(i, j, z, copyjun12)


# SEPTEMBER 2012

sep12.raw <- read.table(file = "D:/sep12.txt", header = TRUE, sep = "\t")
sep12.oka <- read.table(file = "D:/skeletonsep12.txt", header = TRUE, sep = "\t")
copysep12 <- sep12.oka

## filling the skeleton of the new dataframe with data
for (i in 1:nrow(sep12.oka)){
  for (j in 1:nrow(sep12.raw)){
    if (sep12.oka$block[i] == sep12.raw$block[j] && sep12.oka$treatment[i] == sep12.raw$treatment[j] &&
        sep12.oka$datapoint[i] == sep12.raw$datapoint[j]){
      for (z in 1:length(species)){
        if (sep12.raw$species[j] == species[z]){
          copysep12[i, z+11] <- sep12.raw$cover[j]
        }
      }
    }
  }
}

## writes a txt file that must be carefully opened with Excel (separator = space!)
sep12.oka <- copysep12
write.table(sep12.oka, file = "sep12.txt", na = "NA", row.names = FALSE)
### remember to change the location of the dataframes so we don't post them on GitHub!

## cleaning behind
rm(i, j, z, copysep12)


# MAY 2013

may13.raw <- read.table(file = "D:/may13.txt", header = TRUE, sep = "\t")
may13.oka <- read.table(file = "D:/skeletonmay13.txt", header = TRUE, sep = "\t")
copymay13 <- may13.oka

## filling the skeleton of the new dataframe with data
for (i in 1:nrow(may13.oka)){
  for (j in 1:nrow(may13.raw)){
    if (may13.oka$block[i] == may13.raw$block[j] && may13.oka$treatment[i] == may13.raw$treatment[j] &&
        may13.oka$datapoint[i] == may13.raw$datapoint[j]){
      for (z in 1:length(species)){
        if (may13.raw$species[j] == species[z]){
          copymay13[i, z+11] <- may13.raw$cover[j]
        }
      }
    }
  }
}

## writes a txt file that must be carefully opened with Excel (separator = space!)
may13.oka <- copymay13
write.table(may13.oka, file = "may13.txt", na = "NA", row.names = FALSE)
### remember to change the location of the dataframes so we don't post them on GitHub!

## cleaning behind
rm(i, j, z, copymay13)


# JUNE 2013

jun13.raw <- read.table(file = "D:/jun13.txt", header = TRUE, sep = "\t")
jun13.oka <- read.table(file = "D:/skeletonjun13.txt", header = TRUE, sep = "\t")
copyjun13 <- jun13.oka

## filling the skeleton of the new dataframe with data
for (i in 1:nrow(jun13.oka)){
  if (jun13.oka$datapoint[i] == 1){
    for (j in 1:nrow(jun13.raw)){
      if (jun13.oka$block[i] == jun13.raw$block[j] && jun13.oka$treatment[i] == jun13.raw$treatment[j]){
        for (z in 1:length(species)){
          if (species[z] == jun13.raw$species[j]){
            copyjun13[i, z+11] <- jun13.raw[j, 4]
            copyjun13[i+1, z+11] <- jun13.raw[j, 5]
            copyjun13[i+2, z+11] <- jun13.raw[j, 6]
            copyjun13[i+3, z+11] <- jun13.raw[j, 7]
            copyjun13[i+4, z+11] <- jun13.raw[j, 8]
            copyjun13[i+5, z+11] <- jun13.raw[j, 9]
            copyjun13[i+6, z+11] <- jun13.raw[j, 10]
            copyjun13[i+7, z+11] <- jun13.raw[j, 11]
            copyjun13[i+8, z+11] <- jun13.raw[j, 12]
          }
        }
      }
    }
  }
}

## writes a txt file that must be carefully opened with Excel (separator = space!)
jun13.oka <- copyjun13
write.table(jun13.oka, file = "jun13.txt", na = "NA", row.names = FALSE)
### remember to change the location of the dataframes so we don't post them on GitHub!

## cleaning behind
rm(i, j, z, copyjun13)


# OCTOBER 2013 - COVER

oct_cover13.raw <- read.table(file = "D:/oct_cover13.txt", header = TRUE, sep = "\t")
oct_cover13.oka <- read.table(file = "D:/skeletonoct_cover13.txt", header = TRUE, sep = "\t")
copyoct_cover13 <- oct_cover13.oka

## filling the skeleton of the new dataframe with data
for (i in 1:nrow(oct_cover13.oka)){
  if (oct_cover13.oka$datapoint[i] == 1){
    for (j in 1:nrow(oct_cover13.raw)){
      if (oct_cover13.oka$block[i] == oct_cover13.raw$block[j] &&
          oct_cover13.oka$treatment[i] == oct_cover13.raw$treatment[j]){
        for (z in 1:length(species)){
          if (species[z] == oct_cover13.raw$species[j]){
            copyoct_cover13[i, z+11] <- oct_cover13.raw[j, 4]
            copyoct_cover13[i+1, z+11] <- oct_cover13.raw[j, 5]
            copyoct_cover13[i+2, z+11] <- oct_cover13.raw[j, 6]
            copyoct_cover13[i+3, z+11] <- oct_cover13.raw[j, 7]
            copyoct_cover13[i+4, z+11] <- oct_cover13.raw[j, 8]
            copyoct_cover13[i+5, z+11] <- oct_cover13.raw[j, 9]
            copyoct_cover13[i+6, z+11] <- oct_cover13.raw[j, 10]
            copyoct_cover13[i+7, z+11] <- oct_cover13.raw[j, 11]
            copyoct_cover13[i+8, z+11] <- oct_cover13.raw[j, 12]
          }
        }
      }
    }
  }
}

## writes a txt file that must be carefully opened with Excel (separator = space!)
oct_cover13.oka <- copyoct_cover13
write.table(oct_cover13.oka, file = "oct_cover13.txt", na = "NA", row.names = FALSE)
### remember to change the location of the dataframes so we don't post them on GitHub!

##cleaning behind
rm(i, j, z, copyoct_cover13)


# OCTOBER 2013 - DAMAGE

oct_damage13.raw <- read.table(file = "D:/oct_damage13.txt", header = TRUE, sep = "\t")
oct_damage13.oka <- read.table(file = "D:/skeletonoct_damage13.txt", header = TRUE, sep = "\t")
copyoct_damage13 <- oct_damage13.oka

## filling the skeleton of the new dataframe with data
for (i in 1:nrow(oct_damage13.oka)){
  if (oct_damage13.oka$datapoint[i] == 1){
    for (j in 1:nrow(oct_damage13.raw)){
      if (oct_damage13.oka$block[i] == oct_damage13.raw$block[j] &&
          oct_damage13.oka$treatment[i] == oct_damage13.raw$treatment[j]){
        for (z in 1:length(species)){
          if (species[z] == oct_damage13.raw$species[j]){
            copyoct_damage13[i, z+11] <- oct_damage13.raw[j, 4]
            copyoct_damage13[i+1, z+11] <- oct_damage13.raw[j, 5]
            copyoct_damage13[i+2, z+11] <- oct_damage13.raw[j, 6]
            copyoct_damage13[i+3, z+11] <- oct_damage13.raw[j, 7]
            copyoct_damage13[i+4, z+11] <- oct_damage13.raw[j, 8]
            copyoct_damage13[i+5, z+11] <- oct_damage13.raw[j, 9]
            copyoct_damage13[i+6, z+11] <- oct_damage13.raw[j, 10]
            copyoct_damage13[i+7, z+11] <- oct_damage13.raw[j, 11]
            copyoct_damage13[i+8, z+11] <- oct_damage13.raw[j, 12]
          }
        }
      }
    }
  }
}

## writes a txt file that must be carefully opened with Excel (separator = space!)
oct_damage13.oka <- copyoct_damage13
write.table(oct_damage13.oka, file = "oct_damage13.txt", na = "NA", row.names = FALSE)
### remember to change the location of the dataframes so we don't post them on GitHub!

##cleaning behind
rm(i, j, z, copyoct_damage13)


# MAY 2014

may14.raw <- read.table(file = "D:/may14.txt", header = TRUE, sep = "\t")
may14.oka <- read.table(file = "D:/skeletonmay14.txt", header = TRUE, sep = "\t")
copymay14 <- may14.oka

## filling the skeleton of the new dataframe with data
for (i in 1:nrow(may14.oka)){
  if (may14.oka$datapoint[i] == 1){
    for (j in 1:nrow(may14.raw)){
      if (may14.oka$block[i] == may14.raw$block[j] && may14.oka$treatment[i] == may14.raw$treatment[j]){
        for (z in 1:length(species)){
          if (species[z] == may14.raw$species[j]){
            copymay14[i, z+11] <- may14.raw[j, 4]
            copymay14[i+1, z+11] <- may14.raw[j, 5]
            copymay14[i+2, z+11] <- may14.raw[j, 6]
            copymay14[i+3, z+11] <- may14.raw[j, 7]
            copymay14[i+4, z+11] <- may14.raw[j, 8]
            copymay14[i+5, z+11] <- may14.raw[j, 9]
            copymay14[i+6, z+11] <- may14.raw[j, 10]
            copymay14[i+7, z+11] <- may14.raw[j, 11]
            copymay14[i+8, z+11] <- may14.raw[j, 12]
          }
        }
      }
    }
  }
}

## writes a txt file that must be carefully opened with Excel (separator = space!)
may14.oka <- copymay14
write.table(may14.oka, file = "may14.txt", na = "NA", row.names = FALSE)
### remember to change the location of the dataframes so we don't post them on GitHub!

## cleaning behind
rm(i, j, z, copymay14)


# COMBINING ALL DATASETS AND EXPORTING

full_dataset <- rbind(jun12.oka, sep12.oka, may13.oka, jun13.oka, oct_cover13.oka, oct_damage13.oka, may14.oka)
write.table(full_dataset, file = "French_grasshoppers.txt", na = "NA", row.names = FALSE)




# CHANGING THE WHOLE DATABASE AGAIN
## We need to create a new database in order to analyse the data properly. We need a column for the response in cover (in time t + 1)
## and another column for the explanative variable (cover in time t), plus the amount of grasshoppers in time t.

old <- read.table(file = "C:/Users/Granjel/Desktop/Nico_Gross/Data_French_grasshoppers.txt", header = TRUE, sep = "\t")
old_data <- subset(old, type == "cover")

# "d_time", "block", "datapoint", "Cb", "Cd", "Ci", "Ee", "Pg", "Pp",
#                       "ACHMIL", "AGREUP", "ANTODO", "ARRELA", "BROERE", "BROMOL", "BROSTE", "CARCAR", "CENJAC",
#                       "CONARV", "CREPIS", "CRULAE", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FALVUL", "FESARU",
#                       "FESRUB", "FRAEXE", "GALAPA", "GALMOL", "GALVER", "GERDIS", "GERROT", "HIMHIR", "LEUVUL",
#                       "LOLPER", "LOTCOR", "MALSYL", "MEDARA", "MEDLUP", "MYORAM", "ONOREP", "ORCHID", "PICECH",
#                       "PICHIE", "PLALAN", "PLAMAJ", "POAANG", "POAPOI", "POAPRA", "POA.SP", "POATRI", "POTREP",
#                       "PRIVUL", "PRUVUL", "RANACR", "RANREP", "RUBFRU", "RUMACE", "SALPRA", "SENJAC", "SONCHU",
#                       "TAROFF", "TRAPRA", "TRICAM", "TRIFLA", "TRIPRA", "TRIREP", "VERBOF", "VERPER", "VICSAT",
#                       "d_ACHMIL", "d_AGREUP", "d_ANTODO", "d_ARRELA", "d_BROERE", "d_BROMOL", "d_BROSTE", "d_CARCAR", "d_CENJAC",
#                       "d_CONARV", "d_CREPIS", "d_CRULAE", "d_DACGLO", "d_DAUCAR", "d_ELYREP", "d_ERYNGE", "d_FALVUL", "d_FESARU",
#                       "d_FESRUB", "d_FRAEXE", "d_GALAPA", "d_GALMOL", "d_GALVER", "d_GERDIS", "d_GERROT", "d_HIMHIR", "d_LEUVUL",
#                       "d_LOLPER", "d_LOTCOR", "d_MALSYL", "d_MEDARA", "d_MEDLUP", "d_MYORAM", "d_ONOREP", "d_ORCHID", "d_PICECH",
#                       "d_PICHIE", "d_PLALAN", "d_PLAMAJ", "d_POAANG", "d_POAPOI", "d_POAPRA", "d_POA.SP", "d_POATRI", "d_POTREP",
#                       "d_PRIVUL", "d_PRUVUL", "d_RANACR", "d_RANREP", "d_RUBFRU", "d_RUMACE", "d_SALPRA", "d_SENJAC", "d_SONCHU",
#                       "d_TAROFF", "d_TRAPRA", "d_TRICAM", "d_TRIFLA", "d_TRIPRA", "d_TRIREP", "d_VERBOF", "d_VERPER", "d_VICSAT")

before <- subset(old_data, time < 6, select = c(time, block, treatment, datapoint, ACHMIL, AGREUP, ANTODO, ARRELA, BROERE,
                                                BROMOL, BROSTE, CARCAR, CENJAC, CONARV, CREPIS, CRULAE, DACGLO, DAUCAR, ELYREP,
                                                ERYNGE, FALVUL, FESARU, FESRUB, FRAEXE, GALAPA, GALMOL, GALVER, GERDIS, GERROT,
                                                HIMHIR, LEUVUL, LOLPER, LOTCOR, MALSYL, MEDARA, MEDLUP, MYORAM, ONOREP, ORCHID,
                                                PICECH, PICHIE, PLALAN, PLAMAJ, POAANG, POAPOI, POAPRA, POA.SP, POATRI, POTREP,
                                                PRIVUL, PRUVUL, RANACR, RANREP, RUBFRU, RUMACE, SALPRA, SENJAC, SONCHU, TAROFF,
                                                TRAPRA, TRICAM, TRIFLA, TRIPRA, TRIREP, VERBOF, VERPER, VICSAT))

d_ACHMIL <- subset(old_data, time > 1, select = ACHMIL)
d_AGREUP <- subset(old_data, time > 1, select = AGREUP)
d_ANTODO <- subset(old_data, time > 1, select = ANTODO)
d_ARRELA <- subset(old_data, time > 1, select = ARRELA)
d_BROERE <- subset(old_data, time > 1, select = BROERE)
d_BROMOL <- subset(old_data, time > 1, select = BROMOL)
d_BROSTE <- subset(old_data, time > 1, select = BROSTE)
d_CARCAR <- subset(old_data, time > 1, select = CARCAR)
d_CENJAC <- subset(old_data, time > 1, select = CENJAC)
d_CONARV <- subset(old_data, time > 1, select = CONARV)
d_CREPIS <- subset(old_data, time > 1, select = CREPIS)
d_CRULAE <- subset(old_data, time > 1, select = CRULAE)
d_DACGLO <- subset(old_data, time > 1, select = DACGLO)
d_DAUCAR <- subset(old_data, time > 1, select = DAUCAR)
d_ELYREP <- subset(old_data, time > 1, select = ELYREP)
d_ERYNGE <- subset(old_data, time > 1, select = ERYNGE)
d_FALVUL <- subset(old_data, time > 1, select = FALVUL)
d_FESARU <- subset(old_data, time > 1, select = FESARU)
d_FESRUB <- subset(old_data, time > 1, select = FESRUB)
d_FRAEXE <- subset(old_data, time > 1, select = FRAEXE)
d_GALAPA <- subset(old_data, time > 1, select = GALAPA)
d_GALMOL <- subset(old_data, time > 1, select = GALMOL)
d_GALVER <- subset(old_data, time > 1, select = GALVER)
d_GERDIS <- subset(old_data, time > 1, select = GERDIS)
d_GERROT <- subset(old_data, time > 1, select = GERROT)
d_HIMHIR <- subset(old_data, time > 1, select = HIMHIR)
d_LEUVUL <- subset(old_data, time > 1, select = LEUVUL)
d_LOLPER <- subset(old_data, time > 1, select = LOLPER)
d_LOTCOR <- subset(old_data, time > 1, select = LOTCOR)
d_MALSYL <- subset(old_data, time > 1, select = MALSYL)
d_MEDARA <- subset(old_data, time > 1, select = MEDARA)
d_MEDLUP <- subset(old_data, time > 1, select = MEDLUP)
d_MYORAM <- subset(old_data, time > 1, select = MYORAM)
d_ONOREP <- subset(old_data, time > 1, select = ONOREP)
d_ORCHID <- subset(old_data, time > 1, select = ORCHID)
d_PICECH <- subset(old_data, time > 1, select = PICECH)
d_PICHIE <- subset(old_data, time > 1, select = PICHIE)
d_PLALAN <- subset(old_data, time > 1, select = PLALAN)
d_PLAMAJ <- subset(old_data, time > 1, select = PLAMAJ)
d_POAANG <- subset(old_data, time > 1, select = POAANG)
d_POAPOI <- subset(old_data, time > 1, select = POAPOI)
d_POAPRA <- subset(old_data, time > 1, select = POAPRA)
d_POA.SP <- subset(old_data, time > 1, select = POA.SP)
d_POATRI <- subset(old_data, time > 1, select = POATRI)
d_POTREP <- subset(old_data, time > 1, select = POTREP)
d_PRIVUL <- subset(old_data, time > 1, select = PRIVUL)
d_PRUVUL <- subset(old_data, time > 1, select = PRUVUL)
d_RANACR <- subset(old_data, time > 1, select = RANACR)
d_RANREP <- subset(old_data, time > 1, select = RANREP)
d_RUBFRU <- subset(old_data, time > 1, select = RUBFRU)
d_RUMACE <- subset(old_data, time > 1, select = RUMACE)
d_SALPRA <- subset(old_data, time > 1, select = SALPRA)
d_SENJAC <- subset(old_data, time > 1, select = SENJAC)
d_SONCHU <- subset(old_data, time > 1, select = SONCHU)
d_TAROFF <- subset(old_data, time > 1, select = TAROFF)
d_TRAPRA <- subset(old_data, time > 1, select = TRAPRA)
d_TRICAM <- subset(old_data, time > 1, select = TRICAM)
d_TRIFLA <- subset(old_data, time > 1, select = TRIFLA)
d_TRIPRA <- subset(old_data, time > 1, select = TRIPRA)
d_TRIREP <- subset(old_data, time > 1, select = TRIREP)
d_VERBOF <- subset(old_data, time > 1, select = VERBOF)
d_VERPER <- subset(old_data, time > 1, select = VERPER)
d_VICSAT <- subset(old_data, time > 1, select = VICSAT)


new_data <- cbind(before, d_ACHMIL , d_AGREUP , d_ANTODO , d_ARRELA , d_BROERE , d_BROMOL , d_BROSTE , d_CARCAR , d_CENJAC , d_CONARV , d_CREPIS , d_CRULAE , d_DACGLO , d_DAUCAR , d_ELYREP , d_ERYNGE , d_FALVUL , d_FESARU , d_FESRUB , d_FRAEXE , d_GALAPA , d_GALMOL , d_GALVER , d_GERDIS , d_GERROT , d_HIMHIR , d_LEUVUL , d_LOLPER , d_LOTCOR , d_MALSYL , d_MEDARA , d_MEDLUP , d_MYORAM , d_ONOREP , d_ORCHID , d_PICECH , d_PICHIE , d_PLALAN , d_PLAMAJ , d_POAANG , d_POAPOI , d_POAPRA , d_POA.SP , d_POATRI , d_POTREP , d_PRIVUL , d_PRUVUL , d_RANACR , d_RANREP , d_RUBFRU , d_RUMACE , d_SALPRA , d_SENJAC , d_SONCHU , d_TAROFF , d_TRAPRA , d_TRICAM , d_TRIFLA , d_TRIPRA , d_TRIREP , d_VERBOF , d_VERPER , d_VICSAT)

names(new_data) <- c("d_time", "block", "treatment", "datapoint", "ACHMIL", "AGREUP", "ANTODO", "ARRELA", "BROERE", "BROMOL", "BROSTE", "CARCAR", "CENJAC", "CONARV", "CREPIS", "CRULAE", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FALVUL", "FESARU", "FESRUB", "FRAEXE", "GALAPA", "GALMOL", "GALVER", "GERDIS", "GERROT", "HIMHIR", "LEUVUL", "LOLPER", "LOTCOR", "MALSYL", "MEDARA", "MEDLUP", "MYORAM", "ONOREP", "ORCHID", "PICECH", "PICHIE", "PLALAN", "PLAMAJ", "POAANG", "POAPOI", "POAPRA", "POA.SP", "POATRI", "POTREP", "PRIVUL", "PRUVUL", "RANACR", "RANREP", "RUBFRU", "RUMACE", "SALPRA", "SENJAC", "SONCHU", "TAROFF", "TRAPRA", "TRICAM", "TRIFLA", "TRIPRA", "TRIREP", "VERBOF", "VERPER", "VICSAT", "d_ACHMIL", "d_AGREUP", "d_ANTODO", "d_ARRELA", "d_BROERE", "d_BROMOL", "d_BROSTE", "d_CARCAR", "d_CENJAC", "d_CONARV", "d_CREPIS", "d_CRULAE", "d_DACGLO", "d_DAUCAR", "d_ELYREP", "d_ERYNGE", "d_FALVUL", "d_FESARU", "d_FESRUB", "d_FRAEXE", "d_GALAPA", "d_GALMOL", "d_GALVER", "d_GERDIS", "d_GERROT", "d_HIMHIR", "d_LEUVUL", "d_LOLPER", "d_LOTCOR", "d_MALSYL", "d_MEDARA", "d_MEDLUP", "d_MYORAM", "d_ONOREP", "d_ORCHID", "d_PICECH", "d_PICHIE", "d_PLALAN", "d_PLAMAJ", "d_POAANG", "d_POAPOI", "d_POAPRA", "d_POA.SP", "d_POATRI", "d_POTREP", "d_PRIVUL", "d_PRUVUL", "d_RANACR", "d_RANREP", "d_RUBFRU", "d_RUMACE", "d_SALPRA", "d_SENJAC", "d_SONCHU", "d_TAROFF", "d_TRAPRA", "d_TRICAM", "d_TRIFLA", "d_TRIPRA", "d_TRIREP", "d_VERBOF", "d_VERPER", "d_VICSAT")

write.table(new_data, file = "New_French_grashoppers.txt", na = "NA", row.names = FALSE, sep = "\t")

summary(new_data)





