# Reshaping the French grasshoppers dataset - adding the effect of neighbours
French_grasshoppers <- read.table(file = "C:/Users/Granjel RR/Desktop/Nico Gross/Data_French_grasshoppers.txt", header = TRUE, sep = "\t")
summary(French_grasshoppers)

species <- c("ACHMIL", "AGREUP", "ANTODO", "ARRELA", "BROERE", "BROMOL", "BROSTE", "CARCAR", "CENJAC",
             "CONARV", "CREPIS", "CRULAE", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FALVUL", "FESARU",
             "FESRUB", "FRAEXE", "GALAPA", "GALMOL", "GALVER", "GERDIS", "GERROT", "HIMHIR", "LEUVUL",
             "LOLPER", "LOTCOR", "MALSYL", "MEDARA", "MEDLUP", "MYORAM", "ONOREP", "ORCHID", "PICECH",
             "PICHIE", "PLALAN", "PLAMAJ", "POAANG", "POAPOI", "POAPRA", "POA.SP", "POATRI", "POTREP",
             "PRIVUL", "PRUVUL", "RANACR", "RANREP", "RUBFRU", "RUMACE", "SALPRA", "SENJAC", "SONCHU",
             "TAROFF", "TRAPRA", "TRICAM", "TRIFLA", "TRIPRA", "TRIREP", "VERBOF", "VERPER", "VICSAT")


time <- NULL
date <- NULL
block <- NULL
treatment <- NULL
datapoint <- NULL
Cb <- NULL
Cd <- NULL
Ci <- NULL
Ee <- NULL
Pg <- NULL
Pp <- NULL
spp <- NULL
cover <- NULL

#creating the "heading" of the dataframe (factors)
#this is still not correct because I don't have the missing rows corresponding to the dead/disappeared indivs.
#but, at least, it is working
for (i in 1:nrow(French_grasshoppers)){
  for (j in 13:ncol(French_grasshoppers)){
    if (!is.na(French_grasshoppers[i,j])){
      time <- c(time, French_grasshoppers$time[i])
      date <- c(date, French_grasshoppers$date[i])
      block <- c(block, French_grasshoppers$block[i])
      treatment <- c(treatment, French_grasshoppers$treatment[i])
      datapoint <- c(datapoint, French_grasshoppers$datapoint[i])
      Cb <- c(Cb, French_grasshoppers$Cb[i])
      Cd <- c(Cd, French_grasshoppers$Cd[i])
      Ci <- c(Ci, French_grasshoppers$Ci[i])
      Ee <- c(Ee, French_grasshoppers$Ee[i])
      Pg <- c(Pg, French_grasshoppers$Pg[i])
      Pp <- c(Pp, French_grasshoppers$Pp[i])
      spp <- c(spp, species[j-12])
      cover <- c(cover, French_grasshoppers[i,j])
    }
  }
}
new <- data.frame(cbind(time, date, block, treatment, datapoint, Cb, Cd, Ci, Ee, Pg, Pp, spp, cover))


##
shell <- data.frame(matrix(0, 23944, 126))

#this shit is not working... I don't actually know what I'm doing
for (i in 1:nrow(French_grasshoppers)){
  #datapoint1
  if (French_grasshoppers$datapoint[i] == 1){
    for (j in 13:ncol(French_grasshoppers)){
      if (!is.na(French_grasshoppers[i,j])){
        for (w in 1:nrow(French_grasshoppers)){
          if (French_grasshoppers$time[i] == French_grasshoppers$time[w] &&
              French_grasshoppers$date[i] == French_grasshoppers$date[w] &&
              French_grasshoppers$block[i] == French_grasshoppers$block[w] &&
              French_grasshoppers$treatment[i] == French_grasshoppers$treatment[w]){
            if (French_grasshoppers$datapoint[w] == 1 ||
                French_grasshoppers$datapoint[w] == 2 ||
                French_grasshoppers$datapoint[w] == 6){
              for (z in 13:ncol(French_grasshoppers)){
                if (!is.na(French_grasshoppers[w,z])){
                  shell[i,z-12] <- shell[i,z-12] + 1
                  shell[i,z-11] <- shell[i,z-11] + (French_grasshoppers[w,z])
                }
              }
            }
          }
        }
      }
    }
  }
  print(i) #cheater, to see the process and how long it has left
}


#HELP!