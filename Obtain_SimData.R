rm(list = ls())

library(readxl)
Supplementary_Table <- read_excel("Fine_Mapping_Simulation/Data/Supplementary_Table.xlsx", sheet = "TS9", skip = 11)
FunctionalScores_Sim1 <- Supplementary_Table[Supplementary_Table$GWAS_Locus == "36_15q21.3",]
FunctionalScores_Sim2 <- Supplementary_Table[Supplementary_Table$GWAS_Locus == "11_6p21.33",]

FunctionalScores_Sim1$Chr <- unlist(lapply(strsplit(FunctionalScores_Sim1$`Variant_ID (hg19)`,":"),function(x){as.numeric(x[1])}))
FunctionalScores_Sim1$Position <- unlist(lapply(strsplit(FunctionalScores_Sim1$`Variant_ID (hg19)`,":"),function(x){as.numeric(x[2])}))

all_chr.tag <- read.delim("/data/DCEG_shared/data/simulated/simulated_multi_ances_genotype_600K/EUR/all_chr.tag.bim", header=FALSE)

bim_sim1 <- all_chr.tag[(all_chr.tag[,1] %in% FunctionalScores_Sim1$Chr) & (all_chr.tag[,4] %in% FunctionalScores_Sim1$Position),]
bim_sim2 <- all_chr.tag[(all_chr.tag[,1] %in% FunctionalScores_Sim2$Chr) & (all_chr.tag[,4] %in% FunctionalScores_Sim2$Position),]

write.table(bim_sim1[,2],file="Fine_Mapping_Simulation/Data/extract_list_sim1.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(bim_sim2[,2],file="Fine_Mapping_Simulation/Data/extract_list_sim2.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)

system("/data/williamsjacr/software/plink --bfile /data/DCEG_shared/data/simulated/simulated_multi_ances_genotype_600K/EUR/all_chr.tag --extract Fine_Mapping_Simulation/Data/extract_list_sim1.txt --make-bed --out Fine_Mapping_Simulation/Data/Sim1")
system("/data/williamsjacr/software/plink --bfile /data/DCEG_shared/data/simulated/simulated_multi_ances_genotype_600K/EUR/all_chr.tag --extract Fine_Mapping_Simulation/Data/extract_list_sim2.txt --make-bed --out Fine_Mapping_Simulation/Data/Sim2")

