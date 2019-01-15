
#   .: NEW SCRIPT Grasshoppers :.  
##################################### 26.06.18

French_grasshoppers <- read.table(file = "/Users/An/Documents/DB_Rodri/Data_French_grasshoppers.txt", header = TRUE, sep = "\t")

# Datapoints: combinations
datapoint1_inc <- c(1,2,6)
datapoint2_inc <- c(2,1,3,5)
datapoint3_inc <- c(3,2,4)
datapoint4_inc <- c(4,3,5,9)
datapoint5_inc <- c(5,2,4,6,8)
datapoint6_inc <- c(6,1,5,7)
datapoint7_inc <- c(7,6,8)
datapoint8_inc <- c(8,5,7,9)
datapoint9_inc <- c(9,4,8)

species <- colnames(French_grasshoppers)[13:length(colnames(French_grasshoppers))]
blocks <- unique(French_grasshoppers$block)
treatments <- unique(French_grasshoppers$treatment)


# Empty dataframe to fill in
French_grasshoppers_converted <- as.data.frame(matrix(ncol=139 , nrow=50000))

# Headers
species_i <- paste0(species, "_i")
vector_sps <- c()
for(i in 1:length(species_i)){
  vector_sps <- c(vector_sps,species_i[i])
  vector_sps <- c(vector_sps,species[i])
}
colnames(French_grasshoppers_converted) <- c(colnames(French_grasshoppers)[c(1:5,7:12)],"Focal","Cover",vector_sps)

n = 0

# Loop per time in the dataframe (from 1 to 6)
for(time_i in unique(French_grasshoppers$time)){
  subtable <- French_grasshoppers[which(French_grasshoppers$time==time_i),]
  date_i <- unique(subtable$date)
  print(time_i)
  
  # Loop per blocks (from 1 to 5)
  for(block_i in blocks){
    subtable_block <- subtable[which(subtable$block==block_i),]
    print(block_i)
    
    # Loop per treatment (from 1 to 14)
    for(treatment_i in treatments){
      subtable_block_treat <- subtable_block[which(subtable_block$treatment==treatment_i),]
      print(treatment_i)
      
      # Loop per datapoint (from 1 to 9)
      for(row_i in 1:nrow(subtable_block_treat)){
        actual_datapoint <- French_grasshoppers$datapoint[row_i]
        
        # Dataframe restricted to plants:
        plants_df <- subtable_block_treat[,13:ncol(subtable_block_treat)]
        plants_avail_in_datapoint <- species[which(is.na(plants_df[row_i,])==FALSE)]
        
        #Considering each plant available in the datapoint as focal:
        for(plant in plants_avail_in_datapoint){
          n = n+1
          
          # Filling general info
          French_grasshoppers_converted$time[n] <- time_i
          French_grasshoppers_converted$date[n] <- as.character(date_i)
          French_grasshoppers_converted$block[n] <- block_i
          French_grasshoppers_converted$treatment[n] <- as.character(treatment_i)
          French_grasshoppers_converted$datapoint[n] <- actual_datapoint
          French_grasshoppers_converted$Cb[n] <- unique(subtable_block_treat$Cb)
          French_grasshoppers_converted$Cd[n] <- unique(subtable_block_treat$Cd)
          French_grasshoppers_converted$Ci[n] <- unique(subtable_block_treat$Ci)
          French_grasshoppers_converted$Ee[n] <- unique(subtable_block_treat$Ee)
          French_grasshoppers_converted$Pg[n] <- unique(subtable_block_treat$Pg)
          French_grasshoppers_converted$Pp[n] <- unique(subtable_block_treat$Pp)
          
          # Filling plats info
          focal <- plant
          French_grasshoppers_converted$Focal[n] <- focal
          French_grasshoppers_converted$Cover[n] <- plants_df[row_i,plant]
          
          # Selecting the correponding datapoints
          if(actual_datapoint == 1){
            datapoints_involved <- datapoint1_inc
          }else if(actual_datapoint == 2){
            datapoints_involved <- datapoint2_inc
          }else if(actual_datapoint == 3){
            datapoints_involved <- datapoint3_inc
          }else if(actual_datapoint == 4){
            datapoints_involved <- datapoint4_inc
          }else if(actual_datapoint == 5){
            datapoints_involved <- datapoint5_inc
          }else if(actual_datapoint == 6){
            datapoints_involved <- datapoint6_inc
          }else if(actual_datapoint == 7){
            datapoints_involved <- datapoint7_inc
          }else if(actual_datapoint == 8){
            datapoints_involved <- datapoint8_inc
          }else if(actual_datapoint == 9){
            datapoints_involved <- datapoint9_inc
          }
          
          
          # Loop to fill information of each species
          for(species_i in species){
            
            # Subtable with the datapoints to consider in this case
            subtable_block_treat_datapoints <- subtable_block_treat[which(subtable_block_treat$datapoint %in% datapoints_involved),]
            
            if(species_i == focal){
              # Case A. "The species is Focal", so we don't count the actual_datapoint for this
              subtable_block_treat_datap_without_actual_datapoint <- subtable_block_treat_datapoints[which(subtable_block_treat_datapoints$datapoint != actual_datapoint),]
              num_sp <- length(which(is.na(subtable_block_treat_datap_without_actual_datapoint[,species_i]) == FALSE))
              subtable_block_treat_datap_without_actual_datapoint[is.na(subtable_block_treat_datap_without_actual_datapoint)] <- 0
              percent_sp<- sum(subtable_block_treat_datap_without_actual_datapoint[,species_i])
              
              name <- paste0(species_i,"_i")
              French_grasshoppers_converted[n, name] <- num_sp
              French_grasshoppers_converted[n, species_i] <- percent_sp
              
            }else{
              # Case B. "The species is NOT Focal", so we count the information of all data_points involved
              num_sp <- length(which(is.na(subtable_block_treat_datapoints[,species_i]) == FALSE))
              subtable_block_treat_datapoints[is.na(subtable_block_treat_datapoints)] <- 0
              percent_sp<- sum(subtable_block_treat_datapoints[,species_i])
              
              name <- paste0(species_i,"_i")
              French_grasshoppers_converted[n, name] <- num_sp
              French_grasshoppers_converted[n, species_i] <- percent_sp
            }
            
          } #end loop to fill information of each species
          
          
        } #end loop considering as focal each plant available in the data point
        
        
      } #end loop going along the different datapoints (1 to 9)
      
      
    } #end loop going along the different treatments
    
    
  } #end loop along the different blocks
  
  
} #end loop different times


dim(French_grasshoppers_converted[which(is.na(French_grasshoppers_converted$time)==FALSE),]) #23944   139
French_grasshoppers_converted_step1 <- French_grasshoppers_converted[which(is.na(French_grasshoppers_converted$time)==FALSE),]

write.table(French_grasshoppers_converted_step1, paste0("/Users/An/Documents/DB_Rodri/French_grasshoppers_converted_step1_nona.txt"), sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)



