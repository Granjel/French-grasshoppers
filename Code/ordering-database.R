# Initializing

species <- c("ACHMIL", "AGREUP", "ANTODO", "ARRELA", "BROERE", "BROMOL", "BROSTE", "CARCAR", "CENJAC",
             "CONARV", "CREPIS", "CRULAE", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FALVUL", "FESARU",
             "FESRUB", "FRAEXE", "GALAPA", "GALMOL", "GALVER", "GERDIS", "GERROT", "HIMHIR", "LEUVUL",
             "LOLPER", "LOTCOR", "MALSYL", "MEDARA", "MYORAM", "ONOREP", "ORCHID", "PICECH", "PICHIE",
             "PLALAN", "PLAMAJ", "POAANG", "POAPOI", "POAPRA", "POA.SP", "POATRI", "POTREP", "PRIVUL",
             "PRUVUL", "RANACR", "RANREP", "RUBFRU", "RUMACE", "SALPRA", "SENJAC", "SONCHU", "TAROFF",
             "TRAPRA", "TRICAM", "TRIFLA", "TRIPRA", "TRIREP", "VERBOF", "VERPER", "VICSAT")

spp <- read.table(file = "C:/Users/Granjel/Desktop/Nico_Gross/Raw_data_Gross/plantilla_species.txt", header = TRUE, sep = "\t")
template <- read.table(file = "C:/Users/Granjel/Desktop/Nico_Gross/Raw_data_Gross/plantilla.txt", header = TRUE, sep = "\t")
jun12.raw <- read.table(file = "C:/Users/Granjel/Desktop/Nico_Gross/Raw_data_Gross/jun12.txt", header = TRUE, sep = "\t")
sep12.raw <- read.table(file = "C:/Users/Granjel/Desktop/Nico_Gross/Raw_data_Gross/sep12.txt", header = TRUE, sep = "\t")
may13.raw <- read.table(file = "C:/Users/Granjel/Desktop/Nico_Gross/Raw_data_Gross/may13.txt", header = TRUE, sep = "\t")
jul13.raw <- read.table(file = "C:/Users/Granjel/Desktop/Nico_Gross/Raw_data_Gross/jul13.txt", header = TRUE, sep = "\t")
oct13.raw <- read.table(file = "C:/Users/Granjel/Desktop/Nico_Gross/Raw_data_Gross/oct13_cover.txt", header = TRUE, sep = "\t")
may14.raw <- read.table(file = "C:/Users/Granjel/Desktop/Nico_Gross/Raw_data_Gross/may14.txt", header = TRUE, sep = "\t")


# RESHAPING DATASET

for (i in 1:nrow(template)){
  #JUNE 2012
  if (template$time[i] == 1){
    for (ja in 1:nrow(jun12.raw)){
      if (template$block[i] == jun12.raw$block[ja] &&
          template$treatment[i] == jun12.raw$treatment[ja] &&
          template$datapoint[i] == jun12.raw$datapoint[ja]){
        for (za in 1:length(species)){
          if (species[za] == jun12.raw$species[ja]){
            spp[i, za] <- jun12.raw$cover[ja]
          }
        }
      }
    }
  } else {
    #SEPTEMBER 2012
    if (template$time[i] == 2){
      for (jb in 1:nrow(sep12.raw)){
        if (template$block[i] == sep12.raw$block[jb] &&
            template$treatment[i] == sep12.raw$treatment[jb] &&
            template$datapoint[i] == sep12.raw$datapoint[jb]){
          for (zb in 1:length(species)){
            if (species[zb] == sep12.raw$species[jb]){
              spp[i, zb] <- sep12.raw$cover[jb]
            }
          }
        }
      }
    } else {
      #MAY 2013
      if (template$time[i] == 3){
        for (jc in 1:nrow(may13.raw)){
          if (template$block[i] == may13.raw$block[jc] &&
              template$treatment[i] == may13.raw$treatment[jc] &&
              template$datapoint[i] == may13.raw$datapoint[jc]){
            for (zc in 1:length(species)){
              if (species[zc] == may13.raw$species[jc]){
                spp[i, zc] <- may13.raw$cover[jc]
              }
            }
          }
        }
      } else {
        # JULY 2013
        if (template$time[i] == 4){
          for (jd in 1:nrow(jul13.raw)){
            if (template$block[i] == jul13.raw$block[jd] &&
                template$treatment[i] == jul13.raw$treatment[jd] &&
                template$datapoint[i] == 1){
              for (zd in 1:length(species)){
                if (species[zd] == jul13.raw$species[jd]){
                  spp[i, zd] <- jul13.raw[jd, 4]
                  spp[i+1, zd] <- jul13.raw[jd, 5]
                  spp[i+2, zd] <- jul13.raw[jd, 6]
                  spp[i+3, zd] <- jul13.raw[jd, 7]
                  spp[i+4, zd] <- jul13.raw[jd, 8]
                  spp[i+5, zd] <- jul13.raw[jd, 9]
                  spp[i+6, zd] <- jul13.raw[jd, 10]
                  spp[i+7, zd] <- jul13.raw[jd, 11]
                  spp[i+8, zd] <- jul13.raw[jd, 12]
                }
              }
            }
          }
        } else {
          # OCTOBER 2013
          if (template$time[i] == 5){
            for (je in 1:nrow(oct13.raw)){
              if (template$block[i] == oct13.raw$block[je] &&
                  template$treatment[i] == oct13.raw$treatment[je] &&
                  template$datapoint[i] == 1){
                for (ze in 1:length(species)){
                  if (species[ze] == oct13.raw$species[je]){
                    spp[i, ze] <- oct13.raw[je, 4]
                    spp[i+1, ze] <- oct13.raw[je, 5]
                    spp[i+2, ze] <- oct13.raw[je, 6]
                    spp[i+3, ze] <- oct13.raw[je, 7]
                    spp[i+4, ze] <- oct13.raw[je, 8]
                    spp[i+5, ze] <- oct13.raw[je, 9]
                    spp[i+6, ze] <- oct13.raw[je, 10]
                    spp[i+7, ze] <- oct13.raw[je, 11]
                    spp[i+8, ze] <- oct13.raw[je, 12]
                  }
                }
              }
            }
          } else {
            # MAY 2014
            if (template$time[i] == 6){
              for (jf in 1:nrow(may14.raw)){
                if (template$block[i] == may14.raw$block[jf] &&
                    template$treatment[i] == may14.raw$treatment[jf] &&
                    template$datapoint[i] == 1){
                  for (zf in 1:length(species)){
                    if (species[zf] == may14.raw$species[jf]){
                      spp[i, zf] <- may14.raw[jf, 4]
                      spp[i+1, zf] <- may14.raw[jf, 5]
                      spp[i+2, zf] <- may14.raw[jf, 6]
                      spp[i+3, zf] <- may14.raw[jf, 7]
                      spp[i+4, zf] <- may14.raw[jf, 8]
                      spp[i+5, zf] <- may14.raw[jf, 9]
                      spp[i+6, zf] <- may14.raw[jf, 10]
                      spp[i+7, zf] <- may14.raw[jf, 11]
                      spp[i+8, zf] <- may14.raw[jf, 12]
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
} #end


#Combining tables and exporting -- FULL DATASET

French_grasshoppers <- cbind(template, spp)
write.table(French_grasshoppers, file = "French_grasshoppers.txt", na = "NA", row.names = FALSE, sep = "\t")


### DON'T RUN IF NOT NEEDED TO CHECK:
####### beggining - CODE TO CHECK IF I MADE THIS IS RIGHT!

jun12 <- subset(French_grasshoppers, time == 1)
sep12 <- subset(French_grasshoppers, time == 2)
may13 <- subset(French_grasshoppers, time == 3)
jul13 <- subset(French_grasshoppers, time == 4)
oct13 <- subset(French_grasshoppers, time == 5)
may14 <- subset(French_grasshoppers, time == 6)

####### from time 1 to time 3
check <- 0
for (i in 1:nrow(may13.raw)){
  if (may13.raw$species[i] == "VICSAT"){
    check <- check + may13.raw$cover[i]
  }
}
check

sum(may13$VICSAT, na.rm = T)

####### from time 4 to time 6
check <- 0
for (i in 1:nrow(may14.raw)){
  if (may14.raw$species[i] == "VICSAT"){
    check <- check + (sum(may14.raw$X1[i], may14.raw$X2[i], may14.raw$X3[i], may14.raw$X4[i], 
      may14.raw$X5[i], may14.raw$X6[i], may14.raw$X7[i], may14.raw$X8[i], may14.raw$X9[i], na.rm = T))
  }
}
check

sum(may14$VICSAT, na.rm = T)
####### end - CODE TO CHECK IF I MADE THIS IS RIGHT!



# CHANGING THE WHOLE DATABASE AGAIN:
## We need to create a new database in order to analyse the data properly. We need a column for the response in cover (in time t + 1)
## and another column for the explanative variable (cover in time t), plus the amount of grasshoppers in time t.
old_data <- French_grasshoppers

### selecting the time t
before <- subset(old_data, time < 6, select = c(time, block, treatment, datapoint, ACHMIL, AGREUP, ANTODO, ARRELA, BROERE,
                                                BROMOL, BROSTE, CARCAR, CENJAC, CONARV, CREPIS, CRULAE, DACGLO, DAUCAR, ELYREP,
                                                ERYNGE, FALVUL, FESARU, FESRUB, FRAEXE, GALAPA, GALMOL, GALVER, GERDIS, GERROT,
                                                HIMHIR, LEUVUL, LOLPER, LOTCOR, MALSYL, MEDARA, MYORAM, ONOREP, ORCHID,
                                                PICECH, PICHIE, PLALAN, PLAMAJ, POAANG, POAPOI, POAPRA, POA.SP, POATRI, POTREP,
                                                PRIVUL, PRUVUL, RANACR, RANREP, RUBFRU, RUMACE, SALPRA, SENJAC, SONCHU, TAROFF,
                                                TRAPRA, TRICAM, TRIFLA, TRIPRA, TRIREP, VERBOF, VERPER, VICSAT))

### selecting the time t+1
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

### creating the new dataframe...
new_data <- cbind(before, d_ACHMIL , d_AGREUP , d_ANTODO , d_ARRELA , d_BROERE , d_BROMOL , d_BROSTE , d_CARCAR , d_CENJAC , d_CONARV , d_CREPIS , d_CRULAE , d_DACGLO , d_DAUCAR , d_ELYREP , d_ERYNGE , d_FALVUL , d_FESARU , d_FESRUB , d_FRAEXE , d_GALAPA , d_GALMOL , d_GALVER , d_GERDIS , d_GERROT , d_HIMHIR , d_LEUVUL , d_LOLPER , d_LOTCOR , d_MALSYL , d_MEDARA , d_MYORAM , d_ONOREP , d_ORCHID , d_PICECH , d_PICHIE , d_PLALAN , d_PLAMAJ , d_POAANG , d_POAPOI , d_POAPRA , d_POA.SP , d_POATRI , d_POTREP , d_PRIVUL , d_PRUVUL , d_RANACR , d_RANREP , d_RUBFRU , d_RUMACE , d_SALPRA , d_SENJAC , d_SONCHU , d_TAROFF , d_TRAPRA , d_TRICAM , d_TRIFLA , d_TRIPRA , d_TRIREP , d_VERBOF , d_VERPER , d_VICSAT)
### changing the names...
names(new_data) <- c("d_time", "block", "treatment", "datapoint", "ACHMIL", "AGREUP", "ANTODO", "ARRELA", "BROERE", "BROMOL", "BROSTE", "CARCAR", "CENJAC", "CONARV", "CREPIS", "CRULAE", "DACGLO", "DAUCAR", "ELYREP", "ERYNGE", "FALVUL", "FESARU", "FESRUB", "FRAEXE", "GALAPA", "GALMOL", "GALVER", "GERDIS", "GERROT", "HIMHIR", "LEUVUL", "LOLPER", "LOTCOR", "MALSYL", "MEDARA", "MYORAM", "ONOREP", "ORCHID", "PICECH", "PICHIE", "PLALAN", "PLAMAJ", "POAANG", "POAPOI", "POAPRA", "POA.SP", "POATRI", "POTREP", "PRIVUL", "PRUVUL", "RANACR", "RANREP", "RUBFRU", "RUMACE", "SALPRA", "SENJAC", "SONCHU", "TAROFF", "TRAPRA", "TRICAM", "TRIFLA", "TRIPRA", "TRIREP", "VERBOF", "VERPER", "VICSAT", "d_ACHMIL", "d_AGREUP", "d_ANTODO", "d_ARRELA", "d_BROERE", "d_BROMOL", "d_BROSTE", "d_CARCAR", "d_CENJAC", "d_CONARV", "d_CREPIS", "d_CRULAE", "d_DACGLO", "d_DAUCAR", "d_ELYREP", "d_ERYNGE", "d_FALVUL", "d_FESARU", "d_FESRUB", "d_FRAEXE", "d_GALAPA", "d_GALMOL", "d_GALVER", "d_GERDIS", "d_GERROT", "d_HIMHIR", "d_LEUVUL", "d_LOLPER", "d_LOTCOR", "d_MALSYL", "d_MEDARA", "d_MYORAM", "d_ONOREP", "d_ORCHID", "d_PICECH", "d_PICHIE", "d_PLALAN", "d_PLAMAJ", "d_POAANG", "d_POAPOI", "d_POAPRA", "d_POA.SP", "d_POATRI", "d_POTREP", "d_PRIVUL", "d_PRUVUL", "d_RANACR", "d_RANREP", "d_RUBFRU", "d_RUMACE", "d_SALPRA", "d_SENJAC", "d_SONCHU", "d_TAROFF", "d_TRAPRA", "d_TRICAM", "d_TRIFLA", "d_TRIPRA", "d_TRIREP", "d_VERBOF", "d_VERPER", "d_VICSAT")
### and saving it!
write.table(new_data, file = "New_French_grashoppers.txt", na = "NA", row.names = FALSE, sep = "\t")

