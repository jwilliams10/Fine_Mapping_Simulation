rm(list = ls())
set.seed(1330)
gc()

library(bigsnpr)
library(dplyr)

fam_file <- read.table("/data/DCEG_shared/data/simulated/simulated_multi_ances_genotype_600K/EUR/all_chr.tag.fam", quote="\"", comment.char="")


if(!file.exists("Fine_Mapping_Simulation/Data/Sim1.rds")){
  rds <- bigsnpr::snp_readBed("Fine_Mapping_Simulation/Data/Sim1.bed")  
}

# Loading the data from backing files
common_variants_Sim1 <- snp_attach("Fine_Mapping_Simulation/Data/Sim1.rds")

Sim1_Bim <- read.delim("Fine_Mapping_Simulation/Data/Sim1.bim", header=FALSE)

causal_SNPs_Sim1 <- which(Sim1_Bim[,4] == 56432339)

beta_snps_Sim1 <- rnorm(1,mean = 0,sqrt(0.01/1))

g_snps_Sim1 <- common_variants_Sim1$genotypes[,causal_SNPs_Sim1,drop = FALSE]
scaled_causal_snps_Sim1 <- scale(g_snps_Sim1)


if(!file.exists("Fine_Mapping_Simulation/Data/Sim2.rds")){
  rds <- bigsnpr::snp_readBed("Fine_Mapping_Simulation/Data/Sim2.bed")  
}

# Loading the data from backing files
common_variants_Sim2 <- snp_attach("Fine_Mapping_Simulation/Data/Sim2.rds")

Sim2_Bim <- read.delim("Fine_Mapping_Simulation/Data/Sim2.bim", header=FALSE)

causal_SNPs_Sim2 <- which(Sim2_Bim[,4] %in% c(167372978,167370971,167442994,167442994,167404689,167412048,167411788))

beta_snps_Sim2 <- rnorm(length(causal_SNPs_Sim2),mean = 0,sqrt(0.01/length(causal_SNPs_Sim2)))

g_snps_Sim2 <- common_variants_Sim2$genotypes[,causal_SNPs_Sim2,drop = FALSE]
scaled_causal_snps_Sim2 <- scale(g_snps_Sim2)


Y_Sim1 <- list()
Y_Sim2 <- list()

for(l in 1:50){
  
  v <- var(scaled_causal_snps_Sim1%*%matrix(beta_snps_Sim1,ncol = 1))/0.01
  Y_hat_1 <- data.frame(IDs = fam_file[,2],Y_hat_Common = scaled_causal_snps_Sim1%*%matrix(beta_snps_Sim1/sqrt(as.numeric(v)),ncol = 1)) 
  
  epsilon <- rnorm(nrow(Y_hat_1),mean = 0,sd = sqrt(1 - 0.01))
  
  Y_tmp <- Y_hat_1$Y_hat_Common + epsilon
  
  Y_Sim1[[l]] <- data.frame(IDs = Y_hat_1$IDs,Y=Y_tmp)
  
  
  v <- var(scaled_causal_snps_Sim2%*%matrix(beta_snps_Sim2,ncol = 1))/0.01
  Y_hat_1 <- data.frame(IDs = fam_file[,2],Y_hat_Common = scaled_causal_snps_Sim2%*%matrix(beta_snps_Sim2/sqrt(as.numeric(v)),ncol = 1)) 
  
  epsilon <- rnorm(nrow(Y_hat_1),mean = 0,sd = sqrt(1 - 0.01))
  
  Y_tmp <- Y_hat_1$Y_hat_Common + epsilon
  
  Y_Sim2[[l]] <- data.frame(IDs = Y_hat_1$IDs,Y=Y_tmp)
  
  print(l)
}

save(Y_Sim1,file = "Fine_Mapping_Simulation/Data/Sim1_Y.Rdata")
save(Y_Sim2,file = "Fine_Mapping_Simulation/Data/Sim2_Y.Rdata")

