Andrea <- read.table(file = "C:/Users/Granjel RR/Desktop/Nico Gross/French_grasshoppers_converted_step1.txt", header = TRUE, sep = "\t")

point <- NULL
cheat <- NULL

#los NAs en cero, simplemente para trabajar mejor,
#y de paso guardamos los datapoints y las sumas totales de cada fila
for (i in 1:nrow(Andrea)){
  for (j in 1:ncol(Andrea)){
    if (is.na(Andrea[i,j])){
      Andrea[i,j] <- 0
      
    }
  }
  point <- c(point, Andrea$datapoint[i])
  cheat <- c(cheat, sum (Andrea[i, seq(13,139,by=2)])) #sumas por fila
}

tabla <- data.frame(point, cheat) #tabla easy

#ahora vamos a comprobar que las coberturas no sean mayores de lo que se supone deben ser
#por ejemplo, si el datapoint es 1, tiene en cuenta 3 datapoints (1, 2 y 6), y por lo tanto se
#supone que no puede sumar más de 30; así con el resto también

fuck <- NULL #para guardar la posición de los que sobrepasen el máximo
for (w in 1:nrow(tabla)){
  #maximo de 30
  if (tabla$point[i] == 1 || tabla$point[i] == 3 || tabla$point[i] == 7 || tabla$point[i] == 9){
    if (tabla$cheat[i] > 30){
      fuck <- c(fuck, w)
    }
  }
  #maximo de 40
  if (tabla$point[i] == 2 || tabla$point[i] == 4 || tabla$point[i] == 6 || tabla$point[i] == 8){
    if (tabla$cheat[i] > 40){
      fuck <- c(fuck, w)
    }
  }
  #máximo de 50
  if (tabla$point[i] == 5){
    if (tabla$cheat[i] > 50){
      fuck <- c(fuck, w)
    }
  }
}

print(fuck) #me sale que TODOS superan el máximo...


#entonces voy a ver la base de datos original;
#cárgala desde tu ruta y tal, yo ya la tenía cargada
#estoy cambiando los NAs por ceros
for (i in 1:nrow(French_grasshoppers)){
  for (j in 1:ncol(French_grasshoppers)){
    if (is.na(French_grasshoppers[i,j])){
      French_grasshoppers[i,j] <- 0
    }
  }
}
rm(i, j)

#y ahora voy a guardar las posiciones de los que sumen más de 10, que se supone que es el máximo
ok <- NULL
wtf <- NULL
mama <- NULL
for (z in 1:(nrow(French_grasshoppers))){
  ok <- sum(French_grasshoppers[z, 13:75])
  wtf <- c(wtf, ok)
  if (ok > 10){
    mama <- c(mama, z)
  }
}

#resulta que muchas veces suma más de 10... incluso en el summary de abajo se ve que, en alguna ocasión,
#la suma alcanza 23... así que no sé!
print(wtf)
print(mama)
print(cheat)

#voy a hacer un boxplot. Se supone que el resultado de la transformación debería dar 3.67 veces los valores
#de la base de datos original, ya que 1/9 de las veces se suman 5 datapoints, 4/9 se suman 4...
#...y 4/9 se suman 3, lo cual da 3.67; por eso, multiplico los valores de las sumas iniciales por ese valor
try <- wtf*3.67

summary(cheat)
summary(wtf)

#y ahora los represento todos en un boxplot
boxplot(wtf, try, cheat, names = c("Original", "Orig. x3.67", "Andrea"))

#no sé qué pensar... puede que la nueva base de datos transformada esté bien? que los valores por encima
#del máximo se deban a que en la propia base de datos original hay muchos valores por encima? Ni idea,
#imagino que la mejor manera de saberlo es "calculando a mano" algunas filas