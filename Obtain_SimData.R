rm(list = ls())

library(readxl)
FunctionalScores_Sim1 <- read_excel("Fine_Mapping_Simulation/Data/FunctionalScores.xlsx",sheet = "lung cancer 1 var")
FunctionalScores_Sim2 <- read_excel("Fine_Mapping_Simulation/Data/FunctionalScores.xlsx",sheet = "lung cancer 6 vars")

FunctionalScores_Sim1$Chr <- unlist(lapply(strsplit(FunctionalScores_Sim1$Variant_ID,":"),function(x){as.numeric(x[1])}))
FunctionalScores_Sim1$Position <- unlist(lapply(strsplit(FunctionalScores_Sim1$Variant_ID,":"),function(x){as.numeric(x[2])}))

FunctionalScores_Sim2$Chr <- unlist(lapply(strsplit(FunctionalScores_Sim2$Variant_ID,":"),function(x){as.numeric(x[1])}))
FunctionalScores_Sim2$Position <- unlist(lapply(strsplit(FunctionalScores_Sim2$Variant_ID,":"),function(x){as.numeric(x[2])}))

all_chr.tag <- read.delim("/data/DCEG_shared/data/simulated/simulated_multi_ances_genotype_600K/EUR/all_chr.tag.bim", header=FALSE)

bim_sim1 <- all_chr.tag[(all_chr.tag[,1] %in% FunctionalScores_Sim1$Chr) & (all_chr.tag[,4] %in% FunctionalScores_Sim1$Position),]
bim_sim2 <- all_chr.tag[(all_chr.tag[,1] %in% FunctionalScores_Sim2$Chr) & (all_chr.tag[,4] %in% FunctionalScores_Sim2$Position),]

write.table(bim_sim1[,2],file="Fine_Mapping_Simulation/Data/extract_list_sim1.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(bim_sim2[,2],file="Fine_Mapping_Simulation/Data/extract_list_sim2.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)

system("/data/williamsjacr/software/plink --bfile /data/DCEG_shared/data/simulated/simulated_multi_ances_genotype_600K/EUR/all_chr.tag --extract Fine_Mapping_Simulation/Data/extract_list_sim1.txt --make-bed --out Fine_Mapping_Simulation/Data/Sim1")
system("/data/williamsjacr/software/plink --bfile /data/DCEG_shared/data/simulated/simulated_multi_ances_genotype_600K/EUR/all_chr.tag --extract Fine_Mapping_Simulation/Data/extract_list_sim2.txt --make-bed --out Fine_Mapping_Simulation/Data/Sim2")

