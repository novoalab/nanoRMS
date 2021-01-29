#####################################################################
######## SCRIPTS FOR PREDICTING Y FROM PAIRED SAMPLE #################
#####################################################################
###################### BY OGUZHAN BEGIK #############################
#####################################################################
#Rscript predict_twoSample_transcript_noknownsites.R test_data/ncRNA_normal_rep1_epinano.csv test_data/ncRNA_heatshock_rep1_epinano.csv test_data/ncRNA_normal_rep2_epinano.csv test_data/ncRNA_heatshock_rep2_epinano.csv


## Loading libraries
library(stringr)
library(dplyr)
library(VennDiagram)

# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)
#Arguments
input1 <- args[1]  #1st variable
input2 <- args[2]  #2nd variable
input3 <- args[3]  #3rd variable
input4 <- args[4]  #3rd variable



############################################
### PART1#### IMPORTING DATA 
############################################

### FILTERS DETERMINED
mis_thr <- 0.137
c_thr <- 0.578


#IMPORT ncRNA DATA
Normal_ncRNA_REP1 <- read.delim(input1 ,sep=",")
Stress_ncRNA_REP1 <- read.delim(input2 ,sep=",")
Normal_ncRNA_REP2 <- read.delim(input3 ,sep=",")
Stress_ncRNA_REP2 <- read.delim(input4,sep=",")


data_manipulation_forbases_ncrna <- function(data,label) {
	chrom <- str_split_fixed(data$X.Ref, n=2, pattern="_")
	data$RNA <- chrom[,2]
	data_s <- subset(data , RNA =="snRNA" | RNA == "snoRNA") 
	final<- vector()
	for (chro in unique(data_s$X.Ref)) {
		subs<- subset(data_s, X.Ref==chro)
		subs2<- subset(subs, pos >14)
		subs3<- subset(subs2, pos <(max(subs$pos)-14))
		final<- rbind(final, subs3)
	}
	data<- subset(final, cov > 10)
	data$chr_pos <- paste(data$X.Ref, data$pos, sep="_")
	bases <- str_split_fixed(data$ACGT_freq, n=4, pattern=":")
	colnames(bases) <- c("A", "C", "G", "T")
	data2 <- cbind(data, bases)
	data3 <- data2[,c("X.Ref", "pos", "chr_pos", "base", "cov", "q_mean", "q_median", "q_std", "ins", "del", "mis", "A", "T", "C", "G")]
	data3$A <- as.numeric(as.character(data3$A))
	data3$T <- as.numeric(as.character(data3$T))
	data3$G <- as.numeric(as.character(data3$G))
	data3$C <- as.numeric(as.character(data3$C))
	data_U <- subset(data3, base=="T")
	data_U$mismatches <- data_U$A+ data_U$C + data_U$G
 	data_U$count <- data_U$A+ data_U$C + data_U$G + data_U$T
	data_U$mis_freq <- data_U$mismatches/data_U$count
	data_U$A_freq <- data_U$A / data_U$mismatches
	data_U$C_freq <- data_U$C / data_U$mismatches
	data_U$G_freq <- data_U$G / data_U$mismatches
	data_U <- data_U %>% mutate(A_freq = coalesce(A_freq, 0))
	data_U <- data_U %>% mutate(C_freq = coalesce(C_freq, 0))
	data_U <- data_U %>% mutate(G_freq = coalesce(G_freq, 0))
	data_U$sample <- label
	data_final <- data_U[,c("X.Ref", "pos","chr_pos","cov", "mis", "C_freq")]
	return(data_final)
}

Normal_ncRNA_REP1_bases <- data_manipulation_forbases_ncrna(Normal_ncRNA_REP1, "Normal_ncRNA_Rep1")
Stress_ncRNA_REP1_bases <- data_manipulation_forbases_ncrna(Stress_ncRNA_REP1, "Stress_ncRNA_Rep1")
Normal_ncRNA_REP2_bases <- data_manipulation_forbases_ncrna(Normal_ncRNA_REP2, "Normal_ncRNA_Rep2")
Stress_ncRNA_REP2_bases <- data_manipulation_forbases_ncrna(Stress_ncRNA_REP2, "Stress_ncRNA_Rep2")




#############################################
# PART3: FILTERING FOR NORMAL vs HEAT #######
#############################################

predict_stress<- function(normal1,stress1, normal2,stress2) { 
	columns <- c("X.Ref", "pos","chr_pos")
	rep1 <- merge(normal1, stress1, by.x=columns, by.y=columns)
	colnames(rep1) <- c(columns, "normal1_cov", "normal1_mis", "normal1_c","stress1_cov", "stress1_mis", "stress1_c")
	rep2 <- merge(normal2, stress2, by.x=columns, by.y=columns)
	colnames(rep2) <- c(columns, "normal2_cov", "normal2_mis", "normal2_c","stress2_cov", "stress2_mis", "stress2_c")
	all <- merge(rep1, rep2,  by.x=columns, by.y=columns)
	all <- all[!duplicated(all$chr_pos.x), ]
	all$diff_mis1 <-  all$stress1_mis - all$normal1_mis
	all$diff_mis2 <-  all$stress2_mis - all$normal2_mis 
	#Coverage filtered
	cov_filter <- subset(all, normal1_cov > 30 & normal2_cov > 30 & stress1_cov > 30 & stress2_cov > 30 )
	print(paste("Positions passing coverage threshold of 30 = ", paste(nrow(cov_filter)), "sites"))
	#Alteration filtered list (Mismatch difference)
	alt_filter <- subset(cov_filter, diff_mis1 > 0.1 & diff_mis2 > 0.1)
	print(paste("Mismatch difference replicably changing (0.1) = ", paste(nrow(alt_filter)), "sites"))
	#Y-like filtered list (Mismatch and C_freq)
	novel_predicted <- subset(alt_filter, stress1_mis > mis_thr & stress1_c > c_thr & stress2_mis > mis_thr & stress2_c > c_thr )
	print(paste("Final Predicted Sites = ", paste(nrow(novel_predicted)), "sites"))
	#Predicted Known HEat
	#Bind all the positions, remove the duplicated ones from the table
	novel_predicted <- novel_predicted[!duplicated(novel_predicted$chr_pos), ]
	novel_predicted$diff_mis1 <- novel_predicted$stress1_mis - novel_predicted$normal1_mis 
	novel_predicted$diff_mis2 <- novel_predicted$stress2_mis - novel_predicted$normal2_mis 
	novel_predicted$mis_normal_merged <- rowMeans(novel_predicted[,c("normal1_mis", "normal2_mis")])
	novel_predicted$mis_stress_merged <- rowMeans(novel_predicted[,c("stress1_mis", "stress2_mis")])
	write.table(novel_predicted, file="predictions_ncRNA_normal_stress.tsv", sep="\t", quote=FALSE, row.names=FALSE)
	return(novel_predicted)
}

normal_vs_stress_predictions <- predict_stress(Normal_ncRNA_REP1_bases, Stress_ncRNA_REP1_bases, Normal_ncRNA_REP2_bases, Stress_ncRNA_REP2_bases)