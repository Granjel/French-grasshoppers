### load dataset:

d <- read.table("Data_Fg/FG.txt", header = TRUE, sep = "\t")


### subsample for different dates

d2 <- d[d$time == "2",] #June 2012 ### changing this command changes everything !!!
summary(d2$Focal)
#d2 <- d2[-c(16:17, 42:43, 52:53, 66:67, 78:79, 84:85)]


### one dataset for each focal species

d_BROERE_gh <- d2[d2$Focal == "BROERE",]
d_ANTODO_gh <- d2[d2$Focal == "ANTODO",]
d_ARRELA_gh <- d2[d2$Focal == "ARRELA",]
d_BROERE_gh <- d2[d2$Focal == "BROERE",]
d_CENJAC_gh <- d2[d2$Focal == "CENJAC",]
d_CONARV_gh <- d2[d2$Focal == "CONARV",]
d_CREPIS_gh <- d2[d2$Focal == "CREPIS",]
d_DACGLO_gh <- d2[d2$Focal == "DACGLO",]
d_DAUCAR_gh <- d2[d2$Focal == "DAUCAR",]
d_ELYREP_gh <- d2[d2$Focal == "ELYREP",]
d_ERYNGE_gh <- d2[d2$Focal == "ERYNGE",]
d_FESARU_gh <- d2[d2$Focal == "FESARU",]
d_FESRUB_gh <- d2[d2$Focal == "FESRUB",]
d_GALVER_gh <- d2[d2$Focal == "GALVER",]
d_GERDIS_gh <- d2[d2$Focal == "GERDIS",]
d_GERROT_gh <- d2[d2$Focal == "GERROT",]
d_LEUVUL_gh <- d2[d2$Focal == "LEUVUL",]
d_LOLPER_gh <- d2[d2$Focal == "LOLPER",]
d_LOTCOR_gh <- d2[d2$Focal == "LOTCOR",]
d_MEDARA_gh <- d2[d2$Focal == "MEDARA",]
d_ONOREP_gh <- d2[d2$Focal == "ONOREP",]
d_PICECH_gh <- d2[d2$Focal == "PICECH",]
d_PICHIE_gh <- d2[d2$Focal == "PICHIE",]
d_PLALAN_gh <- d2[d2$Focal == "PLALAN",]
d_POAANG_gh <- d2[d2$Focal == "POAANG",]
d_POAPRA_gh <- d2[d2$Focal == "POAPRA",]
d_POATRI_gh <- d2[d2$Focal == "POATRI",]
d_RANACR_gh <- d2[d2$Focal == "RANACR",]
d_RUMACE_gh <- d2[d2$Focal == "RUMACE",]
d_SALPRA_gh <- d2[d2$Focal == "SALPRA",]
d_SONCHU_gh <- d2[d2$Focal == "SONCHU",]
d_TAROFF_gh <- d2[d2$Focal == "TAROFF",]
d_TRIFLA_gh <- d2[d2$Focal == "TRIFLA",]
d_TRIPRA_gh <- d2[d2$Focal == "TRIPRA",]
d_VERBOF_gh <- d2[d2$Focal == "VERBOF",]
d_VERPER_gh <- d2[d2$Focal == "VERPER",]


### X_plant, X_gh and Y matrices for each plant species:

X_BROERE_plant <- as.matrix(d_BROERE_gh[, seq(15, 91, by = 2)]) #competition
X_BROERE_plant <- X_BROERE_plant[c(-2, -15, -20, -27, -33, -36),]
X_BROERE_gh <- as.matrix(d_BROERE_gh[6:11]) #depends on the number of links for every species (suppl. mat. Gross)
Y_BROERE_gh <- d_BROERE_gh$Cover #cover

X_ANTODO_plant <- as.matrix(d_ANTODO_gh[, seq(15, 91, by = 2)])
X_ANTODO_gh <- as.matrix(d_ANTODO_gh[6:11])
Y_ANTODO_gh <- d_ANTODO_gh$Cover

X_ARRELA_plant <- as.matrix(d_ARRELA_gh[, seq(15, 91, by = 2)])
X_ARRELA_gh <- as.matrix(d_ARRELA_gh[6:11])
Y_ARRELA_gh <- d_ARRELA_gh$Cover

X_BROERE_plant <- as.matrix(d_BROERE_gh[, seq(15, 91, by = 2)])
X_BROERE_gh <- as.matrix(d_BROERE_gh[6:11])
Y_BROERE_gh <- d_BROERE_gh$Cover

X_CENJAC_plant <- as.matrix(d_CENJAC_gh[, seq(15, 91, by = 2)])
X_CENJAC_gh <- as.matrix(d_CENJAC_gh[6:11])
Y_CENJAC_gh <- d_CENJAC_gh$Cover

X_CONARV_plant <- as.matrix(d_CONARV_gh[, seq(15, 91, by = 2)])
X_CONARV_gh <- as.matrix(d_CONARV_gh[6:11])
Y_CONARV_gh <- d_CONARV_gh$Cover

X_CREPIS_plant <- as.matrix(d_CREPIS_gh[, seq(15, 91, by = 2)])
X_CREPIS_gh <- as.matrix(d_CREPIS_gh[6:11])
Y_CREPIS_gh <- d_CREPIS_gh$Cover

X_DACGLO_plant <- as.matrix(d_DACGLO_gh[, seq(15, 91, by = 2)])
X_DACGLO_gh <- as.matrix(d_DACGLO_gh[6:11])
Y_DACGLO_gh <- d_DACGLO_gh$Cover

X_DAUCAR_plant <- as.matrix(d_DAUCAR_gh[, seq(15, 91, by = 2)])
X_DAUCAR_gh <- as.matrix(d_DAUCAR_gh[6:11])
Y_DAUCAR_gh <- d_DAUCAR_gh$Cover

X_ELYREP_plant <- as.matrix(d_ELYREP_gh[, seq(15, 91, by = 2)])
X_ELYREP_gh <- as.matrix(d_ELYREP_gh[6:11])
Y_ELYREP_gh <- d_ELYREP_gh$Cover

X_ERYNGE_plant <- as.matrix(d_ERYNGE_gh[, seq(15, 91, by = 2)])
X_ERYNGE_gh <- as.matrix(d_ERYNGE_gh[6:11])
Y_ERYNGE_gh <- d_ERYNGE_gh$Cover

X_FESARU_plant <- as.matrix(d_FESARU_gh[, seq(15, 91, by = 2)])
X_FESARU_gh <- as.matrix(d_FESARU_gh[6:11])
Y_FESARU_gh <- d_FESARU_gh$Cover

X_FESRUB_plant <- as.matrix(d_FESRUB_gh[, seq(15, 91, by = 2)])
X_FESRUB_gh <- as.matrix(d_FESRUB_gh[6:11])
Y_FESRUB_gh <- d_FESRUB_gh$Cover

X_GALVER_plant <- as.matrix(d_GALVER_gh[, seq(15, 91, by = 2)])
X_GALVER_gh <- as.matrix(d_GALVER_gh[6:11])
Y_GALVER_gh <- d_GALVER_gh$Cover

X_GERDIS_plant <- as.matrix(d_GERDIS_gh[, seq(15, 91, by = 2)])
X_GERDIS_gh <- as.matrix(d_GERDIS_gh[6:11])
Y_GERDIS_gh <- d_GERDIS_gh$Cover

X_GERROT_plant <- as.matrix(d_GERROT_gh[, seq(15, 91, by = 2)])
X_GERROT_gh <- as.matrix(d_GERROT_gh[6:11])
Y_GERROT_gh <- d_GERROT_gh$Cover

X_LEUVUL_plant <- as.matrix(d_LEUVUL_gh[, seq(15, 91, by = 2)])
X_LEUVUL_gh <- as.matrix(d_LEUVUL_gh[6:11])
Y_LEUVUL_gh <- d_LEUVUL_gh$Cover

X_LOLPER_plant <- as.matrix(d_LOLPER_gh[, seq(15, 91, by = 2)])
X_LOLPER_gh <- as.matrix(d_LOLPER_gh[6:11])
Y_LOLPER_gh <- d_LOLPER_gh$Cover

X_LOTCOR_plant <- as.matrix(d_LOTCOR_gh[, seq(15, 91, by = 2)])
X_LOTCOR_gh <- as.matrix(d_LOTCOR_gh[6:11])
Y_LOTCOR_gh <- d_LOTCOR_gh$Cover

X_MEDARA_plant <- as.matrix(d_MEDARA_gh[, seq(15, 91, by = 2)])
X_MEDARA_gh <- as.matrix(d_MEDARA_gh[6:11])
Y_MEDARA_gh <- d_MEDARA_gh$Cover

X_ONOREP_plant <- as.matrix(d_ONOREP_gh[, seq(15, 91, by = 2)])
X_ONOREP_gh <- as.matrix(d_ONOREP_gh[6:11])
Y_ONOREP_gh <- d_ONOREP_gh$Cover

X_PICECH_plant <- as.matrix(d_PICECH_gh[, seq(15, 91, by = 2)])
X_PICECH_gh <- as.matrix(d_PICECH_gh[6:11])
Y_PICECH_gh <- d_PICECH_gh$Cover

X_PICHIE_plant <- as.matrix(d_PICHIE_gh[, seq(15, 91, by = 2)])
X_PICHIE_gh <- as.matrix(d_PICHIE_gh[6:11])
Y_PICHIE_gh <- d_PICHIE_gh$Cover

X_PLALAN_plant <- as.matrix(d_PLALAN_gh[, seq(15, 91, by = 2)])
X_PLALAN_gh <- as.matrix(d_PLALAN_gh[6:11])
Y_PLALAN_gh <- d_PLALAN_gh$Cover

X_POAANG_plant <- as.matrix(d_POAANG_gh[, seq(15, 91, by = 2)])
X_POAANG_gh <- as.matrix(d_POAANG_gh[6:11])
Y_POAANG_gh <- d_POAANG_gh$Cover

X_POAPRA_plant <- as.matrix(d_POAPRA_gh[, seq(15, 91, by = 2)])
X_POAPRA_gh <- as.matrix(d_POAPRA_gh[6:11])
Y_POAPRA_gh <- d_POAPRA_gh$Cover

X_POATRI_plant <- as.matrix(d_POATRI_gh[, seq(15, 91, by = 2)])
X_POATRI_gh <- as.matrix(d_POATRI_gh[6:11])
Y_POATRI_gh <- d_POATRI_gh$Cover

X_RANACR_plant <- as.matrix(d_RANACR_gh[, seq(15, 91, by = 2)])
X_RANACR_gh <- as.matrix(d_RANACR_gh[6:11])
Y_RANACR_gh <- d_RANACR_gh$Cover

X_RUMACE_plant <- as.matrix(d_RUMACE_gh[, seq(15, 91, by = 2)])
X_RUMACE_gh <- as.matrix(d_RUMACE_gh[6:11])
Y_RUMACE_gh <- d_RUMACE_gh$Cover

X_SALPRA_plant <- as.matrix(d_SALPRA_gh[, seq(15, 91, by = 2)])
X_SALPRA_gh <- as.matrix(d_SALPRA_gh[6:11])
Y_SALPRA_gh <- d_SALPRA_gh$Cover

X_SONCHU_plant <- as.matrix(d_SONCHU_gh[, seq(15, 91, by = 2)])
X_SONCHU_gh <- as.matrix(d_SONCHU_gh[6:11])
Y_SONCHU_gh <- d_SONCHU_gh$Cover

X_TAROFF_plant <- as.matrix(d_TAROFF_gh[, seq(15, 91, by = 2)])
X_TAROFF_gh <- as.matrix(d_TAROFF_gh[6:11])
Y_TAROFF_gh <- d_TAROFF_gh$Cover

X_TRIFLA_plant <- as.matrix(d_TRIFLA_gh[, seq(15, 91, by = 2)])
X_TRIFLA_gh <- as.matrix(d_TRIFLA_gh[6:11])
Y_TRIFLA_gh <- d_TRIFLA_gh$Cover

X_TRIPRA_plant <- as.matrix(d_TRIPRA_gh[, seq(15, 91, by = 2)])
X_TRIPRA_gh <- as.matrix(d_TRIPRA_gh[6:11])
Y_TRIPRA_gh <- d_TRIPRA_gh$Cover

X_VERBOF_plant <- as.matrix(d_VERBOF_gh[, seq(15, 91, by = 2)])
X_VERBOF_gh <- as.matrix(d_VERBOF_gh[6:11])
Y_VERBOF_gh <- d_VERBOF_gh$Cover

X_VERPER_plant <- as.matrix(d_VERPER_gh[, seq(15, 91, by = 2)])
X_VERPER_gh <- as.matrix(d_VERPER_gh[6:11])
Y_VERPER_gh <- d_VERPER_gh$Cover


out_BROERE_gh <- lm(log(Y_BROERE_gh) ~ X_BROERE_plant * X_BROERE_gh)
summary(out_BROERE_gh)
qqnorm(residuals(out_BROERE_gh)); qqline(residuals(out_BROERE_gh))

out_CONARV_gh <- lm(log(Y_CONARV_gh) ~ X_CONARV_plant * X_CONARV_gh)
summary(out_CONARV_gh)
qqnorm(residuals(out_CONARV_gh)); qqline(residuals(out_CONARV_gh))


out_BROERE_gh <- lm(log(Y_BROERE_gh) ~ X_BROERE_plant + X_BROERE_gh)
summary(out_BROERE_gh)
qqnorm(residuals(out_BROERE_gh)); qqline(residuals(out_BROERE_gh))




