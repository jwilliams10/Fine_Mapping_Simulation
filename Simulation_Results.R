rm(list = ls())

theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = 16),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold.italic", size =18),
            #legend.text = element_text(face ="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){   
  library(scales)   
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)    
} 

scale_colour_Publication <- function(...){   
  library(scales)   
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(cowplot)

load("Fine_Mapping_Simulation/Data/Sim1_Result.RData")

Sim1_Bim <- read.delim("Fine_Mapping_Simulation/Data/Sim1.bim", header=FALSE)

causal_SNPs_Sim1 <- which(Sim1_Bim[,4] == 56432339)

load("Fine_Mapping_Simulation/Data/Sim2_Result.RData")

Sim2_Bim <- read.delim("Fine_Mapping_Simulation/Data/Sim2.bim", header=FALSE)

causal_SNPs_Sim2 <- which(Sim2_Bim[,4] %in% c(167372978,167370971,167442994,167442994,167404689,167412048,167411788))


Sim1_Results <- matrix(NA,nrow = 50,ncol = 8)
Sim2_Results <- matrix(NA,nrow = 50,ncol = 8)

i <- 1
for(i in 1:50){
  identified_SNPs <- which(results_susie_rss_pip_Sim1[[i]]$`k$pip` > 0.95)
  Sim1_Results[i,1] <- sum(identified_SNPs %in% causal_SNPs_Sim1)
  Sim1_Results[i,2] <- sum(!(identified_SNPs %in% causal_SNPs_Sim1))
  
  Sim1_Results[i,3] <- sum(results_susie_rss_cs_Sim1[[i]] %in% Sim1_Bim[causal_SNPs_Sim1,2])
  Sim1_Results[i,4] <- sum(!(results_susie_rss_cs_Sim1[[i]] %in% Sim1_Bim[causal_SNPs_Sim1,2]))
  
  
  identified_SNPs <- which(results_susie_rss_pip_Sim2[[i]]$`k$pip` > 0.95)
  Sim2_Results[i,1] <- sum(identified_SNPs %in% causal_SNPs_Sim2)
  Sim2_Results[i,2] <- sum(!(identified_SNPs %in% causal_SNPs_Sim2))
  
  Sim2_Results[i,3] <- sum(results_susie_rss_cs_Sim2[[i]] %in% Sim2_Bim[causal_SNPs_Sim2,2])
  Sim2_Results[i,4] <- sum(!(results_susie_rss_cs_Sim2[[i]] %in% Sim2_Bim[causal_SNPs_Sim2,2]))
  
  
  
  
  
  identified_SNPs <- which(results_FINEMAP_pip_Sim1[[i]]$`output$prob` > 0.95)
  Sim1_Results[i,5] <- sum(identified_SNPs %in% causal_SNPs_Sim1)
  Sim1_Results[i,6] <- sum(!(identified_SNPs %in% causal_SNPs_Sim1))
  
  Sim1_Results[i,7] <- sum(unique(results_FINEMAP_cs_Sim1[[i]]) %in% Sim1_Bim[causal_SNPs_Sim1,2])
  Sim1_Results[i,8] <- sum(!(unique(results_FINEMAP_cs_Sim1[[i]]) %in% Sim1_Bim[causal_SNPs_Sim1,2]))
  
  
  identified_SNPs <- which(results_FINEMAP_pip_Sim2[[i]]$`output$prob` > 0.95)
  Sim2_Results[i,5] <- sum(identified_SNPs %in% causal_SNPs_Sim2)
  Sim2_Results[i,6] <- sum(!(identified_SNPs %in% causal_SNPs_Sim2))
  
  Sim2_Results[i,7] <- sum(unique(results_FINEMAP_cs_Sim2[[i]]) %in% Sim2_Bim[causal_SNPs_Sim2,2])
  Sim2_Results[i,8] <- sum(!(unique(results_FINEMAP_cs_Sim2[[i]]) %in% Sim2_Bim[causal_SNPs_Sim2,2]))
}

Sim1_Results <- as.data.frame(Sim1_Results)
colnames(Sim1_Results) <- c("SuSie_PIP_TP","SuSie_PIP_FP","SuSie_CS_TP","SuSie_CS_FP","FINEMAP_PIP_TP","FINEMAP_PIP_FP","FINEMAP_CS_TP","FINEMAP_CS_FP")

Sim1_Long <- data.frame(Method = rep(c("SuSie_PIP","SuSie_CS","FINEMAP_PIP","FINEMAP_CS"),each = 50),
                        TPs = c(Sim1_Results$SuSie_PIP_TP,Sim1_Results$SuSie_CS_TP,Sim1_Results$FINEMAP_PIP_TP,Sim1_Results$FINEMAP_CS_TP),
                        FPs = c(Sim1_Results$SuSie_PIP_FP,Sim1_Results$SuSie_CS_FP,Sim1_Results$FINEMAP_PIP_FP,Sim1_Results$FINEMAP_CS_FP))

Sim2_Results <- as.data.frame(Sim2_Results)
colnames(Sim2_Results) <- c("SuSie_PIP_TP","SuSie_PIP_FP","SuSie_CS_TP","SuSie_CS_FP","FINEMAP_PIP_TP","FINEMAP_PIP_FP","FINEMAP_CS_TP","FINEMAP_CS_FP")

Sim2_Long <- data.frame(Method = rep(c("SuSie_PIP","SuSie_CS","FINEMAP_PIP","FINEMAP_CS"),each = 50),
                        TPs = c(Sim2_Results$SuSie_PIP_TP,Sim2_Results$SuSie_CS_TP,Sim2_Results$FINEMAP_PIP_TP,Sim2_Results$FINEMAP_CS_TP),
                        FPs = c(Sim2_Results$SuSie_PIP_FP,Sim2_Results$SuSie_CS_FP,Sim2_Results$FINEMAP_PIP_FP,Sim2_Results$FINEMAP_CS_FP))

average_pips_Sim1 <- cbind(Sim1_Bim,rowMeans(do.call(cbind,lapply(results_susie_rss_pip_Sim1,function(x){x$`k$pip`}))),apply(do.call(cbind,lapply(results_susie_rss_pip_Sim1,function(x){x$`k$pip`})),1,function(x){sd(x)/sqrt(length(x))}),
                           rowMeans(do.call(cbind,lapply(results_FINEMAP_pip_Sim1,function(x){x$`output$prob`}))),apply(do.call(cbind,lapply(results_FINEMAP_pip_Sim1,function(x){x$`output$prob`})),1,function(x){sd(x)/sqrt(length(x))}))
average_pips_Sim2 <- cbind(Sim2_Bim,rowMeans(do.call(cbind,lapply(results_susie_rss_pip_Sim2,function(x){x$`k$pip`}))),apply(do.call(cbind,lapply(results_susie_rss_pip_Sim2,function(x){x$`k$pip`})),1,function(x){sd(x)/sqrt(length(x))}),
      rowMeans(do.call(cbind,lapply(results_FINEMAP_pip_Sim2,function(x){x$`output$prob`}))),apply(do.call(cbind,lapply(results_FINEMAP_pip_Sim2,function(x){x$`output$prob`})),1,function(x){sd(x)/sqrt(length(x))}))

### Simulation 1

Supplementary_Table <- read_excel("Fine_Mapping_Simulation/Data/Supplementary_Table.xlsx", sheet = "TS9", skip = 11)
FunctionalScores_Sim1 <- Supplementary_Table[Supplementary_Table$GWAS_Locus == "36_15q21.3",]
FunctionalScores_Sim2 <- Supplementary_Table[Supplementary_Table$GWAS_Locus == "11_6p21.33",]

FunctionalScores_Sim1$Chr <- unlist(lapply(strsplit(FunctionalScores_Sim1$`Variant_ID (hg19)`,":"),function(x){as.numeric(x[1])}))
FunctionalScores_Sim1$Position <- unlist(lapply(strsplit(FunctionalScores_Sim1$`Variant_ID (hg19)`,":"),function(x){as.numeric(x[2])}))

FunctionalScores_Sim1 <- FunctionalScores_Sim1[,c(17,18,2,15)]
colnames(FunctionalScores_Sim1) <- c("Chr","Position","RSID","Integrative_Score")

average_pips_Sim1 <- average_pips_Sim1[,c(1,2,4,7,8,9,10)]
colnames(average_pips_Sim1) <- c("Chr","SNPID","Position","Mean_SuSie_PIP","SE_SuSie_PIP","Mean_FINEMAP_PIP","SE_FINEMAP_PIP")

average_pips_Sim1 <- inner_join(average_pips_Sim1,FunctionalScores_Sim1)
average_pips_Sim1$Causal_SNPs <- "No"
average_pips_Sim1$Causal_SNPs[causal_SNPs_Sim1] <- "Yes"

average_pips_Sim1$Causal_SNPs <- factor(average_pips_Sim1$Causal_SNPs,levels = c("Yes","No"))

plot_dat_susie <- data.frame(SNP = c(average_pips_Sim1$RSID[average_pips_Sim1$Causal_SNPs == "Yes"],"All Other SNPs"),
           Average_PIP = c(mean(average_pips_Sim1$Mean_SuSie_PIP[average_pips_Sim1$Causal_SNPs == "Yes"]),mean(average_pips_Sim1$Mean_SuSie_PIP[average_pips_Sim1$Causal_SNPs == "No"])),
           SE_PIP = c(mean(average_pips_Sim1$SE_SuSie_PIP[average_pips_Sim1$Causal_SNPs == "Yes"]),mean(average_pips_Sim1$SE_SuSie_PIP[average_pips_Sim1$Causal_SNPs == "No"])),
           Average_Integrative_Score = c(mean(c(average_pips_Sim1$Integrative_Score[average_pips_Sim1$Causal_SNPs == "Yes"])),mean(c(average_pips_Sim1$Integrative_Score[average_pips_Sim1$Causal_SNPs == "No"]))),
           SE_Integrative_Score = c(rep(0,sum(average_pips_Sim1$Causal_SNPs == "Yes")),
                                    sd(c(average_pips_Sim1$Integrative_Score[average_pips_Sim1$Causal_SNPs == "No"]))/sqrt(length(c(average_pips_Sim1$Integrative_Score[average_pips_Sim1$Causal_SNPs == "No"])))),
           Causal_SNP = c("Yes","No"))

plot_dat_finemap <- data.frame(SNP = c(average_pips_Sim1$RSID[average_pips_Sim1$Causal_SNPs == "Yes"],"All Other SNPs"),
                             Average_PIP = c(mean(average_pips_Sim1$Mean_FINEMAP_PIP[average_pips_Sim1$Causal_SNPs == "Yes"]),mean(average_pips_Sim1$Mean_FINEMAP_PIP[average_pips_Sim1$Causal_SNPs == "No"])),
                             SE_PIP = c(mean(average_pips_Sim1$SE_FINEMAP_PIP[average_pips_Sim1$Causal_SNPs == "Yes"]),mean(average_pips_Sim1$SE_FINEMAP_PIP[average_pips_Sim1$Causal_SNPs == "No"])),
                             Average_Integrative_Score = c(mean(c(average_pips_Sim1$Integrative_Score[average_pips_Sim1$Causal_SNPs == "Yes"])),mean(c(average_pips_Sim1$Integrative_Score[average_pips_Sim1$Causal_SNPs == "No"]))),
                             SE_Integrative_Score = c(rep(0,sum(average_pips_Sim1$Causal_SNPs == "Yes")),
                                                      sd(c(average_pips_Sim1$Integrative_Score[average_pips_Sim1$Causal_SNPs == "No"]))/sqrt(length(c(average_pips_Sim1$Integrative_Score[average_pips_Sim1$Causal_SNPs == "No"])))),
                             Causal_SNP = c("Yes","No"))

plot1 <- ggplot(plot_dat_susie, aes(SNP,Average_PIP,fill = Causal_SNP)) + geom_col() + ylim(c(0,1)) + theme_Publication() + geom_hline(yintercept = 0.95,linetype=2,col = "red") + ylab("SuSie-RSS Average PIP") + scale_fill_Publication() + coord_flip()
plot2 <- ggplot(plot_dat_finemap, aes(SNP,Average_PIP,fill = Causal_SNP)) + geom_col() + ylim(c(0,1)) + theme_Publication() + geom_hline(yintercept = 0.95,linetype=2,col = "red") + ylab("FINEMAP Average PIP") + scale_fill_Publication() + coord_flip()
plot3 <- ggplot(plot_dat_susie, aes(SNP,Average_Integrative_Score,fill = Causal_SNP)) + geom_col() + theme_Publication()  + ylab("Average Integrative Score") + scale_y_continuous(breaks=c(0,2,4,6,8,10,12)) + scale_fill_Publication() + coord_flip()

legend_b <- get_plot_component(plot1 + guides(fill=guide_legend(title="Causal SNP?")), "guide-box-right")

# Remove PIP in TPs
# Remove TPs for Sim1
# Table TPs and Total Number of SNPs in Credible Sets, no PIP

table_dat1 <- Sim1_Long[Sim1_Long$Method %in% c("SuSie_CS","FINEMAP_CS"),]
aggregate(TPs/length(causal_SNPs_Sim1) ~ Method,data = table_dat1,mean)
aggregate(TPs + FPs ~ Method,data = table_dat1,mean)

### Simulation 2
FunctionalScores_Sim2 <- FunctionalScores_Sim2[,c(17,18,2,15)]
colnames(FunctionalScores_Sim2) <- c("Chr","Position","RSID","Integrative_Score")

average_pips_Sim2 <- average_pips_Sim2[,c(1,2,4,7,8,9,10)]
colnames(average_pips_Sim2) <- c("Chr","SNPID","Position","Mean_SuSie_PIP","SE_SuSie_PIP","Mean_FINEMAP_PIP","SE_FINEMAP_PIP")

average_pips_Sim2 <- inner_join(average_pips_Sim2,FunctionalScores_Sim2)
average_pips_Sim2$Causal_SNPs <- "No"
average_pips_Sim2$Causal_SNPs[causal_SNPs_Sim2] <- "Yes"

average_pips_Sim2$Causal_SNPs <- factor(average_pips_Sim2$Causal_SNPs,levels = c("Yes","No"))

plot_dat_susie <- data.frame(SNP = c(average_pips_Sim2$RSID[average_pips_Sim2$Causal_SNPs == "Yes"],"All Other SNPs"),
                       Average_PIP = c(average_pips_Sim2$Mean_SuSie_PIP[average_pips_Sim2$Causal_SNPs == "Yes"],mean(average_pips_Sim2$Mean_SuSie_PIP[average_pips_Sim2$Causal_SNPs == "No"])),
                       SE_PIP = c(average_pips_Sim2$SE_SuSie_PIP[average_pips_Sim2$Causal_SNPs == "Yes"],mean(average_pips_Sim2$SE_SuSie_PIP[average_pips_Sim2$Causal_SNPs == "No"])),
                       Average_Integrative_Score = c(average_pips_Sim2$Integrative_Score[average_pips_Sim2$Causal_SNPs == "Yes"],mean(average_pips_Sim2$Integrative_Score[average_pips_Sim2$Causal_SNPs == "No"])),
                       SE_Integrative_Score = c(rep(0,sum(average_pips_Sim2$Causal_SNPs == "Yes")),
                                                sd(c(average_pips_Sim2$Integrative_Score[average_pips_Sim2$Causal_SNPs == "No"]))/sqrt(length(c(average_pips_Sim2$Integrative_Score[average_pips_Sim2$Causal_SNPs == "No"])))),
                       Causal_SNP = c("Yes","Yes","Yes","Yes","Yes","Yes","No"))

plot_dat_finemap <- data.frame(SNP = c(average_pips_Sim2$RSID[average_pips_Sim2$Causal_SNPs == "Yes"],"All Other SNPs"),
                             Average_PIP = c(average_pips_Sim2$Mean_FINEMAP_PIP[average_pips_Sim2$Causal_SNPs == "Yes"],mean(average_pips_Sim2$Mean_FINEMAP_PIP[average_pips_Sim2$Causal_SNPs == "No"])),
                             SE_PIP = c(average_pips_Sim2$SE_FINEMAP_PIP[average_pips_Sim2$Causal_SNPs == "Yes"],mean(average_pips_Sim2$SE_FINEMAP_PIP[average_pips_Sim2$Causal_SNPs == "No"])),
                             Average_Integrative_Score = c(average_pips_Sim2$Integrative_Score[average_pips_Sim2$Causal_SNPs == "Yes"],mean(average_pips_Sim2$Integrative_Score[average_pips_Sim2$Causal_SNPs == "No"])),
                             SE_Integrative_Score = c(rep(0,sum(average_pips_Sim2$Causal_SNPs == "Yes")),
                                                      sd(c(average_pips_Sim2$Integrative_Score[average_pips_Sim2$Causal_SNPs == "No"]))/sqrt(length(c(average_pips_Sim2$Integrative_Score[average_pips_Sim2$Causal_SNPs == "No"])))),
                             Causal_SNP = c("Yes","Yes","Yes","Yes","Yes","Yes","No"))

plot4 <- ggplot(plot_dat_susie, aes(SNP,Average_PIP,fill = Causal_SNP)) + geom_col() + ylim(c(0,1)) + theme_Publication() + geom_hline(yintercept = 0.95,linetype=2,col = "red") + scale_fill_Publication() + ylab("SuSie-RSS Average PIP") + coord_flip()
plot5 <- ggplot(plot_dat_finemap, aes(SNP,Average_PIP,fill = Causal_SNP)) + geom_col() + ylim(c(0,1)) + theme_Publication() + geom_hline(yintercept = 0.95,linetype=2,col = "red") + scale_fill_Publication() + ylab("FINEMAP Average PIP") + coord_flip()
plot6 <- ggplot(plot_dat_susie, aes(SNP,Average_Integrative_Score,fill = Causal_SNP)) + geom_col() + theme_Publication()  + ylab("Average Integrative Score") + scale_y_continuous(breaks=c(0,2,4,6,8,10,12)) + scale_fill_Publication() + coord_flip()

pgrid <- plot_grid(
  plot1 + theme(legend.position="none"),
  plot2 + theme(legend.position="none"),
  plot3 + theme(legend.position="none"),
  plot4 + theme(legend.position="none"),
  plot5 + theme(legend.position="none"),
  plot6 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  ncol = 3,
  labels = c("A","B","C","D","E","F")
)

plot_grid(pgrid, legend_b, ncol = 2, rel_widths = c(1, .15))


table_dat2 <- Sim2_Long[Sim2_Long$Method %in% c("SuSie_CS","FINEMAP_CS"),]
aggregate(TPs/length(causal_SNPs_Sim2) ~ Method,data = table_dat2,mean)
aggregate(TPs + FPs ~ Method,data = table_dat2,mean)
