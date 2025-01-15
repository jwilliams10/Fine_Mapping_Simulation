rm(list = ls())

library(susieR)
library(bigsnpr)

load("FineMapping_Application/Data/Sim1_Y.Rdata")

common_variants_Sim1 <- snp_attach("FineMapping_Application/Data/Sim1.rds")

Sim1_Bim <- read.delim("FineMapping_Application/Data/Sim1.bim", header=FALSE)

Sim1_Fam <- read.table("FineMapping_Application/Data/Sim1.fam", quote="\"", comment.char="")
Sim1_Fam[,1] <- 0
write.table(Sim1_Fam,"FineMapping_Application/Data/Sim1.fam", col.names = FALSE,row.names = FALSE)

G_Sim1 <- common_variants_Sim1$genotypes[,1:dim(common_variants_Sim1$genotypes)[2]]

load("FineMapping_Application/Data/Sim2_Y.Rdata")

common_variants_Sim2 <- snp_attach("FineMapping_Application/Data/Sim2.rds")

Sim2_Bim <- read.delim("FineMapping_Application/Data/Sim2.bim", header=FALSE)

G_Sim2 <- common_variants_Sim2$genotypes[,1:dim(common_variants_Sim2$genotypes)[2]]


results_susie_rss_pip_Sim1 <- list()
results_susie_rss_cs_Sim1 <- list()

results_FINEMAP_pip_Sim1 <- list()
results_FINEMAP_cs_Sim1 <- list()

results_susie_rss_pip_Sim2 <- list()
results_susie_rss_cs_Sim2 <- list()

results_FINEMAP_pip_Sim2 <- list()
results_FINEMAP_cs_Sim2 <- list()

i <- 1

for(i in 1:50){
  # Y_Sim1[[i]]$FID <- 0
  # Y_Sim1[[i]] <- Y_Sim1[[i]][,c(1,3,2)]
  # colnames(Y_Sim1[[i]]) <- c("IID","FID","Y")
  # write.table(Y_Sim1[[i]],file = "FineMapping_Application/Data/Y_Sim1.txt",sep = '\t',row.names = FALSE,quote = FALSE)
  # 
  # system("/data/williamsjacr/software/plink2 --bfile FineMapping_Application/Data/Sim1 --pheno FineMapping_Application/Data/Y_Sim1.txt --pheno-name Y --linear allow-no-covars --vif 999 --out FineMapping_Application/Data/Sim1_SumStats")
  # 
  # Sim1_SumStats <- read.delim("FineMapping_Application/Data/Sim1_SumStats.Y.glm.linear", header=FALSE, comment.char="#")
  # colnames(Sim1_SumStats) <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
  # 
  
  
  sumstats <- univariate_regression(G_Sim1, Y_Sim1[[i]]$Y)
  z_scores <- sumstats$betahat / sumstats$sebetahat
  R <- cor(G_Sim1)
  k <- susie_rss(z_scores,R, n=nrow(G_Sim1),L = 10,coverage = 0.95)
  output <- cbind(Sim1_Bim,k$pip)
  results_susie_rss_pip_Sim1[[i]] <- output
  results_susie_rss_cs_Sim1[[i]] <- Sim1_Bim[unlist(k$sets$cs),2]
  
  
  
  z <- data.frame(rsid = paste0("rs",1:ncol(G_Sim1)), chromosome = 1, position = 1:ncol(G_Sim1),allele1 = "A",allele2 = "G",maf = 0.1,beta = sumstats$betahat,se = sumstats$sebetahat)
  dir <- "FineMapping_Application/Data/"
  name_dir <- "Sim1_FINEMAP"
  
  write.table(R,file = paste0(dir,name_dir,".ld"),row.names = FALSE,col.names = FALSE,quote = FALSE)
  write.table(z,file = paste0(dir,name_dir,".z"),row.names = FALSE,col.names = TRUE,quote = FALSE)
  
  tmp1 <- c("z;ld;snp;config;cred;log;k;n_samples")
  tmp2 <- paste0(dir,name_dir,".z;",dir,name_dir,".ld;",dir,name_dir,".snp;",dir,name_dir,".config;",dir,name_dir,".cred;",dir,name_dir,".log;",dir,name_dir,".k;",nrow(G_Sim1))
  write.table(rbind(tmp1,tmp2),row.names = FALSE,col.names = FALSE,quote = FALSE,file = paste0(dir,name_dir))
  
  system(paste0("/data/williamsjacr/software/finemap_v1.4.2_x86_64/./finemap_v1.4.2_x86_64 --sss --in-files FineMapping_Application/Data/Sim1_FINEMAP"))
  
  output <- read.csv("FineMapping_Application/Data/Sim1_FINEMAP.snp", sep="")
  
  output <- cbind(Sim1_Bim,output$prob)
  results_FINEMAP_pip_Sim1[[i]] <- output
  
  output <- read.table("FineMapping_Application/Data/Sim1_FINEMAP.cred1", header=TRUE, quote="\"")
  system("rm FineMapping_Application/Data/Sim1_FINEMAP.cred1")
  cs <- as.numeric(gsub("rs","",output$cred1))
  if(file.exists("FineMapping_Application/Data/Sim1_FINEMAP.cred2")){
    output <- read.table("FineMapping_Application/Data/Sim1_FINEMAP.cred2", header=TRUE, quote="\"")
    system("rm FineMapping_Application/Data/Sim1_FINEMAP.cred2")
    cs <- c(cs,as.numeric(gsub("rs","",output$cred1)))
  }
  results_FINEMAP_cs_Sim1[[i]] <- Sim1_Bim[cs,2]
  
  
  
  sumstats <- univariate_regression(G_Sim2, Y_Sim2[[i]]$Y)
  z_scores <- sumstats$betahat / sumstats$sebetahat
  R <- cor(G_Sim2)
  k <- susie_rss(z_scores,R, n=nrow(G_Sim2),L = 10,coverage = 0.95)
  output <- cbind(Sim2_Bim,k$pip)
  results_susie_rss_pip_Sim2[[i]] <- output
  results_susie_rss_cs_Sim2[[i]] <- Sim2_Bim[unlist(k$sets$cs),2]
  
  
  
  z <- data.frame(rsid = paste0("rs",1:ncol(G_Sim2)), chromosome = 1, position = 1:ncol(G_Sim2),allele1 = "A",allele2 = "G",maf = 0.1,beta = sumstats$betahat,se = sumstats$sebetahat)
  dir <- "FineMapping_Application/Data/"
  name_dir <- "Sim2_FINEMAP"
  
  write.table(R,file = paste0(dir,name_dir,".ld"),row.names = FALSE,col.names = FALSE,quote = FALSE)
  write.table(z,file = paste0(dir,name_dir,".z"),row.names = FALSE,col.names = TRUE,quote = FALSE)
  
  tmp1 <- c("z;ld;snp;config;cred;log;k;n_samples")
  tmp2 <- paste0(dir,name_dir,".z;",dir,name_dir,".ld;",dir,name_dir,".snp;",dir,name_dir,".config;",dir,name_dir,".cred;",dir,name_dir,".log;",dir,name_dir,".k;",nrow(G_Sim2))
  write.table(rbind(tmp1,tmp2),row.names = FALSE,col.names = FALSE,quote = FALSE,file = paste0(dir,name_dir))
  
  system(paste0("/data/williamsjacr/software/finemap_v1.4.2_x86_64/./finemap_v1.4.2_x86_64 --sss --in-files FineMapping_Application/Data/Sim2_FINEMAP"))
  
  output <- read.csv("FineMapping_Application/Data/Sim2_FINEMAP.snp", sep="")
  
  output <- cbind(Sim2_Bim,output$prob)
  results_FINEMAP_pip_Sim2[[i]] <- output
  
  
  output <- read.table("FineMapping_Application/Data/Sim2_FINEMAP.cred2", header=TRUE, quote="\"")
  system("rm FineMapping_Application/Data/Sim2_FINEMAP.cred2")
  cs <- as.numeric(gsub("rs","",output$cred1))
  if(file.exists("FineMapping_Application/Data/Sim2_FINEMAP.cred1")){
    output <- read.table("FineMapping_Application/Data/Sim2_FINEMAP.cred1", header=TRUE, quote="\"")
    system("rm FineMapping_Application/Data/Sim2_FINEMAP.cred1")
    cs <- c(cs,as.numeric(gsub("rs","",output$cred1)))
  }
  if(file.exists("FineMapping_Application/Data/Sim2_FINEMAP.cred3")){
    output <- read.table("FineMapping_Application/Data/Sim2_FINEMAP.cred3", header=TRUE, quote="\"")
    system("rm FineMapping_Application/Data/Sim2_FINEMAP.cred3")
    cs <- c(cs,as.numeric(gsub("rs","",output$cred1)))
  }
  results_FINEMAP_cs_Sim2[[i]] <- Sim2_Bim[cs,2]
  
  print(i)
}

save(results_susie_rss_pip_Sim1,results_susie_rss_cs_Sim1,results_FINEMAP_pip_Sim1,results_FINEMAP_cs_Sim1,file = "FineMapping_Application/Data/Sim1_Result.RData")
save(results_susie_rss_pip_Sim2,results_susie_rss_cs_Sim2,results_FINEMAP_pip_Sim2,results_FINEMAP_cs_Sim2,file = "FineMapping_Application/Data/Sim2_Result.RData")
