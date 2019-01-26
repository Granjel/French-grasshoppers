alpha <- as.matrix(read.table("Results/output_glinternet/alpha.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: competition
lambda <- as.matrix(read.table("Results/output_glinternet/lambda.txt", header = TRUE, sep = "\t"))
sigma <- as.matrix(read.table("Results/output_glinternet/sigma.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma1 <- as.matrix(read.table("Results/output_glinternet/gamma1.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma2 <- as.matrix(read.table("Results/output_glinternet/gamma2.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma3 <- as.matrix(read.table("Results/output_glinternet/gamma3.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma4 <- as.matrix(read.table("Results/output_glinternet/gamma4.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma5 <- as.matrix(read.table("Results/output_glinternet/gamma5.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory
gamma6 <- as.matrix(read.table("Results/output_glinternet/gamma6.txt", header = TRUE, sep = "\t") * (-1)) #negative interaction: herbivory

source("Code/functions_structural_coex_outputs.R") #functions to compute the structural coexistence outputs

#FIELD CONDITIONS
#sigma -> lambda
intr <- lambda + rowSums(sigma) #intrinsic growth rate vector

#gamma(s) -> alpha
comp_alpha <- alpha + (gamma1 + gamma2 + gamma3 + gamma4 + gamma5 + gamma6) #competition alpha coefficients

#output
field_conditions <- structural_coex_3spp(alpha = as.data.frame(comp_alpha), intrinsic = as.matrix(intr)) #calculate the coex. outputs
write.table(field_conditions, file = "Results/structural_coex_field_conditions.txt", sep = "\t", row.names = FALSE)

#NO HERBIVORY
no_herbivory <- structural_coex_3spp(alpha = as.data.frame(alpha), intrinsic = as.matrix(lambda)) #calculate the coex. outputs
write.table(no_herbivory, file = "Results/structural_coex_no_herbivory.txt", sep = "\t", row.names = FALSE)