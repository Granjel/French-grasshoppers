
# Oscar, copia y pega esto en R y cambia el directorio desde el que quieres que se cargue
# el archivo .txt y después corre todo.
full_dataset <- read.table(file = "/Users/oscargodoy/Dropbox/French grasshoppers/Data_French_grasshoppers.txt", header = TRUE, sep = "\t")
library(nlme)
cover <- subset(full_dataset, type=="cover")
ACHMIL <- lme(ARRELA ~ ACHMIL + AGREUP, data = cover,
              random = ~ 1 | block, correlation = corAR1(form = ~ time),
              method='ML', na.action = na.omit)
### No sé por qué no corre -- échale un ojo, porfa; y, además, ¿cuál es el objeto del comando 'control'?
### Porque también me dice que en 'control = lCtr' no se reconoce el objeto 'lCtr'...
### Yo entiendo que lCtr es un objeto de un análisis anterior tuyo, pero habría que elegir otro ahora, ¿no?

 + ANTODO + ARRELA + BROERE + BROMOL + BROSTE + CARCAR + CENJAC +
  CONARV + CREPIS + CRULAE + DACGLO + DAUCAR + ELYREP + ERYNGE + FALVUL + FESARU + FESRUB +
  FRAEXE + GALAPA + GALMOL + GALVER + GERDIS + GERROT + HIMHIR + LEUVUL + LOLPER + LOTCOR +
  MALSYL + MEDARA + MEDLUP + MYORAM + ONOREP + ORCHID + PICECH + PICHIE + PLALAN + PLAMAJ +
  POAANG + POAPOI + POAPRA + POA.SP + POATRI + POTREP + PRIVUL + PRUVUL + RANACR + RANREP +
  RUBFRU + RUMACE + SALPRA + SENJAC + SONCHU + TAROFF + TRAPRA + TRICAM + TRIFLA + TRIPRA +
  TRIREP + VERBOF + VERPER + VICSAT + Cb + Cd + Ci + Ee + Pg + Pp,