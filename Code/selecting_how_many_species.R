########################################################################################
## Transforming the database (x10) and exploring how many species should be selected ##
########################################################################################

#load dataframe
French_grasshoppers <- read.table(file = "C:/Users/Granjel RR/Desktop/Nico Gross/French_grasshoppers_converted_step1_nona.txt", header = TRUE, sep = "\t")

## STEP 1
#multiplying all the cover data x10
for (i in 1:nrow(French_grasshoppers)){
  for (j in seq(13,139,by=2)){
    French_grasshoppers[i,j] <- French_grasshoppers[i,j] * 10
  }
  print((i/23944)*100)
}
#saving dataset
write.table(French_grasshoppers, file = "C:/Users/Granjel RR/Desktop/Nico Gross/French_grasshoppers_converted_step1.txt", sep = "\t")

##STEP 2
# creating a dataframe to save the number of times a species appears in the dataset
time <- c(rep(1, 63), rep(2, 63), rep(3, 63), rep(4, 63), rep(5, 63), rep(6, 63)) #6 dates
spp <- names(French_grasshoppers[,seq(15,139,by=2)]) #species names
spp <- rep(spp, 6) #x6 dates
try <- data.frame(cbind(time, spp)) #binding them all
colnames(try) <- c("time", "species") #changing the col names
rm(spp, time)

#saving how many times per date a species appears in the dataset
pff <- NULL
for (i in 1:nrow(try)){
  for (j in 1:nrow(French_grasshoppers)){
    if ((try$time[i] == French_grasshoppers$time[j]) && (try$species[i] == French_grasshoppers$Focal[j])){
      cuenta <- cuenta + 1
    }
  }
  pff <- c(pff, cuenta)
  cuenta <- 0
  print((i/378)*100)
}

#creating the final dataset with appearance information
sel <- data.frame(cbind(try, pff)) #called sel
colnames(sel) <- c("time", "species", "number")

#subset of dates to calculate the relative abundance of the species within each date
t1 <- subset(sel, sel$time == 1)
t2 <- subset(sel, sel$time == 2)
t3 <- subset(sel, sel$time == 3)
t4 <- subset(sel, sel$time == 4)
t5 <- subset(sel, sel$time == 5)
t6 <- subset(sel, sel$time == 6)
s1 <- sum(t1$number) #sum of individuals t1
s2 <- sum(t2$number) # " t2
s3 <- sum(t3$number) # " t3
s4 <- sum(t4$number) # " t4
s5 <- sum(t5$number) # " t5
s6 <- sum(t6$number) # " t6

#calculating the relative abundance per date of each species
for (i in 1:nrow(sel)){
  if (sel$time[i] == 1){
    sel$number[i] <- (sel$number[i] / s1) * 100
  }
  if (sel$time[i] == 2){
    sel$number[i] <- (sel$number[i] / s2) * 100
  }
  if (sel$time[i] == 3){
    sel$number[i] <- (sel$number[i] / s3) * 100
  }
  if (sel$time[i] == 4){
    sel$number[i] <- (sel$number[i] / s4) * 100
  }
  if (sel$time[i] == 5){
    sel$number[i] <- (sel$number[i] / s5) * 100
  }
  if (sel$time[i] == 6){
    sel$number[i] <- (sel$number[i] / s6) * 100
  }
}

# graph! species relative abundance per date!
p<-ggplot(sel, aes(x=species, y=number, fill=time)) +
  geom_bar(stat = "identity", position="dodge") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  geom_hline(yintercept = 0.5)
print(p)