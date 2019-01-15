library(nlme)

# ACHMIL + AGREUP + ANTODO + ARRELA + BROERE + BROMOL + BROSTE + CARCAR + CENJAC + CONARV + CREPIS + CRULAE + DACGLO + DAUCAR + ELYREP + ERYNGE + FALVUL + FESARU + FESRUB + FRAEXE + GALAPA + GALMOL + GALVER + GERDIS + GERROT + HIMHIR + LEUVUL + LOLPER + LOTCOR + MALSYL + MEDARA + MYORAM + ONOREP + ORCHID + PICECH + PICHIE + PLALAN + PLAMAJ + POAANG + POAPOI + POAPRA + POA.SP + POATRI + POTREP + PRIVUL + PRUVUL + RANACR + RANREP + RUBFRU + RUMACE + SALPRA + SENJAC + SONCHU + TAROFF + TRAPRA + TRICAM + TRIFLA + TRIPRA + TRIREP + VERBOF + VERPER + VICSAT

French_grasshoppers_dt <- read.table(file = "C:/Users/Granjel RR/Desktop/Nico Gross/French_grasshoppers_dt.txt", header = TRUE, sep = "\t")

summary(lme(d_BROERE ~ log(Cb+1) + log(Cd+1) + log(Ci+1) + log(Ee+1) + log(Pg+1) + log(Pp+1),
            data = French_grasshoppers_dt, random = ~ 1 | block, method = "REML", na.action = na.omit))

summary(lme(d_ACHMIL ~ ACHMIL + BROERE,
            data = French_grasshoppers_dt, random = ~ 1 | block, method = "REML", na.action = na.omit))

cover_oct <- subset(French_grasshoppers_dt, d_time == 5)
cor.test(try, herbivory$Total, na.action = na.omit)

try <- NULL
for (i in 1:nrow(cover_oct)){
  x <- sum(cover_oct[i, (11:73)], na.rm = TRUE)
  try <- c(try, x)
}

#####



data_BROERE <- subset(French_grasshoppers_dt, !is.na(d_BROERE))
for (i in 1:nrow(data_BROERE)){
  for (j in 1:ncol(data_BROERE)){
    if (is.na(data_BROERE[i,j])){
      data_BROERE[i,j] <- 0
    }
  }
}

summary(lme(d_BROERE ~ BROERE + ARRELA + DACGLO + DAUCAR + PLALAN + POAANG + RANACR + GALVER +
              TRIPRA + SALPRA + CONARV + MEDARA + GERDIS + FESRUB + TRIFLA + ERYNGE + LEUVUL +
              TAROFF + ONOREP + PICECH + CENJAC + ELYREP + ACHMIL + CREPIS + POAPRA + POATRI +
              RUMACE + SONCHU + PICHIE + FESARU + ANTODO + VERBOF + VERPER + GERROT + LOTCOR,
            data = data_BROERE, random = ~ 1 | block, method = "REML"))
## need to add the temporal correlation

# 1st, 





