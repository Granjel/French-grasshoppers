# BROERE
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


# ARRELA
data_ARRELA <- subset(French_grasshoppers_dt, !is.na(d_ARRELA))
for (i in 1:nrow(data_ARRELA)){
  for (j in 1:ncol(data_ARRELA)){
    if (is.na(data_ARRELA[i,j])){
      data_ARRELA[i,j] <- 0
    }
  }
}

summary(lme(d_ARRELA ~ BROERE + ARRELA + DACGLO + DAUCAR + PLALAN + POAANG + RANACR + GALVER +
              TRIPRA + SALPRA + CONARV + MEDARA + GERDIS + FESRUB + TRIFLA + ERYNGE + LEUVUL +
              TAROFF + ONOREP + PICECH + CENJAC + ELYREP + ACHMIL + CREPIS + POAPRA + POATRI +
              RUMACE + SONCHU + PICHIE + FESARU + ANTODO + VERBOF + VERPER + GERROT + LOTCOR,
            data = data_ARRELA, random = ~ 1 | block, method = "REML"))


# DACGLO
data_DACGLO <- subset(French_grasshoppers_dt, !is.na(d_DACGLO))
for (i in 1:nrow(data_DACGLO)){
  for (j in 1:ncol(data_DACGLO)){
    if (is.na(data_DACGLO[i,j])){
      data_DACGLO[i,j] <- 0
    }
  }
}

summary(lme(d_DACGLO ~ BROERE + ARRELA + DACGLO + DAUCAR + PLALAN + POAANG + RANACR + GALVER +
              TRIPRA + SALPRA + CONARV + MEDARA + GERDIS + FESRUB + TRIFLA + ERYNGE + LEUVUL +
              TAROFF + ONOREP + PICECH + CENJAC + ELYREP + ACHMIL + CREPIS + POAPRA + POATRI +
              RUMACE + SONCHU + PICHIE + FESARU + ANTODO + VERBOF + VERPER + GERROT + LOTCOR,
            data = data_DACGLO, random = ~ 1 | block, method = "REML"))


# DAUCAR
data_DAUCAR <- subset(French_grasshoppers_dt, !is.na(d_DAUCAR))
for (i in 1:nrow(data_DAUCAR)){
  for (j in 1:ncol(data_DAUCAR)){
    if (is.na(data_DAUCAR[i,j])){
      data_DAUCAR[i,j] <- 0
    }
  }
}

summary(lme(d_DAUCAR ~ BROERE + ARRELA + DACGLO + DAUCAR + PLALAN + POAANG + RANACR + GALVER +
              TRIPRA + SALPRA + CONARV + MEDARA + GERDIS + FESRUB + TRIFLA + ERYNGE + LEUVUL +
              TAROFF + ONOREP + PICECH + CENJAC + ELYREP + ACHMIL + CREPIS + POAPRA + POATRI +
              RUMACE + SONCHU + PICHIE + FESARU + ANTODO + VERBOF + VERPER + GERROT + LOTCOR,
            data = data_DAUCAR, random = ~ 1 | block, method = "REML"))


# PLALAN
data_PLALAN <- subset(French_grasshoppers_dt, !is.na(d_PLALAN))
for (i in 1:nrow(data_PLALAN)){
  for (j in 1:ncol(data_PLALAN)){
    if (is.na(data_PLALAN[i,j])){
      data_PLALAN[i,j] <- 0
    }
  }
}

summary(lme(d_PLALAN ~ BROERE + ARRELA + DACGLO + DAUCAR + PLALAN + POAANG + RANACR + GALVER +
              TRIPRA + SALPRA + CONARV + MEDARA + GERDIS + FESRUB + TRIFLA + ERYNGE + LEUVUL +
              TAROFF + ONOREP + PICECH + CENJAC + ELYREP + ACHMIL + CREPIS + POAPRA + POATRI +
              RUMACE + SONCHU + PICHIE + FESARU + ANTODO + VERBOF + VERPER + GERROT + LOTCOR,
            data = data_PLALAN, random = ~ 1 | block, method = "REML"))


# POAANG
data_POAANG <- subset(French_grasshoppers_dt, !is.na(d_POAANG))
for (i in 1:nrow(data_POAANG)){
  for (j in 1:ncol(data_POAANG)){
    if (is.na(data_POAANG[i,j])){
      data_POAANG[i,j] <- 0
    }
  }
}

summary(lme(d_POAANG ~ BROERE + ARRELA + DACGLO + DAUCAR + PLALAN + POAANG + RANACR + GALVER +
              TRIPRA + SALPRA + CONARV + MEDARA + GERDIS + FESRUB + TRIFLA + ERYNGE + LEUVUL +
              TAROFF + ONOREP + PICECH + CENJAC + ELYREP + ACHMIL + CREPIS + POAPRA + POATRI +
              RUMACE + SONCHU + PICHIE + FESARU + ANTODO + VERBOF + VERPER + GERROT + LOTCOR,
            data = data_POAANG, random = ~ 1 | block, method = "REML"))


# RANACR
data_RANACR <- subset(French_grasshoppers_dt, !is.na(d_RANACR))
for (i in 1:nrow(data_RANACR)){
  for (j in 1:ncol(data_RANACR)){
    if (is.na(data_RANACR[i,j])){
      data_RANACR[i,j] <- 0
    }
  }
}

summary(lme(d_RANACR ~ BROERE + ARRELA + DACGLO + DAUCAR + PLALAN + POAANG + RANACR + GALVER +
              TRIPRA + SALPRA + CONARV + MEDARA + GERDIS + FESRUB + TRIFLA + ERYNGE + LEUVUL +
              TAROFF + ONOREP + PICECH + CENJAC + ELYREP + ACHMIL + CREPIS + POAPRA + POATRI +
              RUMACE + SONCHU + PICHIE + FESARU + ANTODO + VERBOF + VERPER + GERROT + LOTCOR,
            data = data_RANACR, random = ~ 1 | block, method = "REML"))


# GALVER
data_GALVER <- subset(French_grasshoppers_dt, !is.na(d_GALVER))
for (i in 1:nrow(data_GALVER)){
  for (j in 1:ncol(data_GALVER)){
    if (is.na(data_GALVER[i,j])){
      data_GALVER[i,j] <- 0
    }
  }
}

summary(lme(d_GALVER ~ BROERE + ARRELA + DACGLO + DAUCAR + PLALAN + POAANG + RANACR + GALVER +
              TRIPRA + SALPRA + CONARV + MEDARA + GERDIS + FESRUB + TRIFLA + ERYNGE + LEUVUL +
              TAROFF + ONOREP + PICECH + CENJAC + ELYREP + ACHMIL + CREPIS + POAPRA + POATRI +
              RUMACE + SONCHU + PICHIE + FESARU + ANTODO + VERBOF + VERPER + GERROT + LOTCOR,
            data = data_GALVER, random = ~ 1 | block, method = "REML"))

# IMPORTANT NOTE:
# all the intra-specific coefficients are positive. Does this sound weird? Isn't the intra-specific competition
# thought to be very dependent on the plant species and, thus, it could be both positive or negative?
# Yep, but we are not having something into account - we are not analysing individuals here, we are analysing
# the cover of a certain species present at a certain spot. I don't see how we can obtain the 'intras' here...

