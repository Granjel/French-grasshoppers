# Initializing

species <- c("ACHMIL", "AGREUP", "ANTODO", "ARRELA", "BROERE", "BROMOL", "BROSTE", "CARCAR", "CENJAC",
             "CONARV", "CREPIS", "CRULAE", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FALVUL", "FESARU",
             "FESRUB", "FRAEXE", "GALAPA", "GALMOL", "GALVER", "GERDIS", "GERROT", "HIMHIR", "LEUVUL",
             "LOLPER", "LOTCOR", "MALSYL", "MEDARA", "MEDLUP", "MYORAM", "ONOREP", "ORCHID", "PICECH",
             "PICHIE", "PLALAN", "PLAMAJ", "POAANG", "POAPOI", "POAPRA", "POA.SP", "POATRI", "POTREP",
             "PRIVUL", "PRUVUL", "RANACR", "RANREP", "RUBFRU", "RUMACE", "SALPRA", "SENJAC", "SONCHU",
             "TAROFF", "TRAPRA", "TRICAM", "TRIFLA", "TRIPRA", "TRIREP", "TRISET", "VERBOF", "VERPER",
             "VICSAT")


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
