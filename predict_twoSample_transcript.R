#####################################################################
######## SCRIPTS FOR PREDICTING Y FROM PAIRED SAMPLE #################
#####################################################################
###################### BY OGUZHAN BEGIK #############################
#####################################################################
#Rscript predict_twoSample_transcript.R ncRNA_normal_rep1_epinano.csv ncRNA_heatshock_rep1_epinano.csv ncRNA_normal_rep2_epinano.csv ncRNA_heatshock_rep2_epinano.csv


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
# Mod File
mod_ncRNA <- read.delim("supp/yeast_ncRNA_Y_Positions.tsv")

### FILTERS DETERMINED
mis_thr <- 0.137
c_thr <- 0.578


#IMPORT ncRNA DATA
Normal_ncRNA_REP1 <- read.delim(input1 ,sep=",")
Heat_ncRNA_REP1 <- read.delim(input2 ,sep=",")
Normal_ncRNA_REP2 <- read.delim(input3 ,sep=",")
Heat_ncRNA_REP2 <- read.delim(input4,sep=",")


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
	data_U2 <- merge(data_U, mod_ncRNA, all.x = TRUE, by.x=c("X.Ref", "pos"), by.y=c("chr", "gene_pos"))
	data_U3 <- data_U2[,c("X.Ref", "pos", "chr_pos.x", "base","sample", "Reference", "Heat", "cov", "q_mean", "q_median", "q_std", "ins", "del", "mis", "A", "T", "C", "G","mis_freq","A_freq","C_freq","G_freq")]
	data_U3$Heat <- as.character(data_U3$Heat)
	data_U3$Heat[is.na(data_U3$Heat)] <- "Unm"
	data_U3$Reference <- as.character(data_U3$Reference)
	data_U3$Reference[is.na(data_U3$Reference)] <- "Unm"
	data_final <- data_U3[,c("X.Ref", "pos","chr_pos.x", "Heat", "Reference","cov", "mis", "C_freq")]
	return(data_final)
}

Normal_ncRNA_REP1_bases <- data_manipulation_forbases_ncrna(Normal_ncRNA_REP1, "Normal_ncRNA_Rep1")
Heat_ncRNA_REP1_bases <- data_manipulation_forbases_ncrna(Heat_ncRNA_REP1, "Heat_ncRNA_Rep1")
Normal_ncRNA_REP2_bases <- data_manipulation_forbases_ncrna(Normal_ncRNA_REP2, "Normal_ncRNA_Rep2")
Heat_ncRNA_REP2_bases <- data_manipulation_forbases_ncrna(Heat_ncRNA_REP2, "Heat_ncRNA_Rep2")




#############################################
# PART3: FILTERING FOR NORMAL vs HEAT #######
#############################################

predict_heat<- function(normal1,heat1, normal2,heat2) { 
	columns <- c("X.Ref", "pos","chr_pos.x", "Heat", "Reference")
	rep1 <- merge(normal1, heat1, by.x=columns, by.y=columns)
	colnames(rep1) <- c(columns, "normal1_cov", "normal1_mis", "normal1_c","heat1_cov", "heat1_mis", "heat1_c")
	rep2 <- merge(normal2, heat2, by.x=columns, by.y=columns)
	colnames(rep2) <- c(columns, "normal2_cov", "normal2_mis", "normal2_c","heat2_cov", "heat2_mis", "heat2_c")
	all <- merge(rep1, rep2,  by.x=columns, by.y=columns)
	all <- all[!duplicated(all$chr_pos.x), ]
	all$diff_mis1 <-  all$heat1_mis - all$normal1_mis
	all$diff_mis2 <-  all$heat2_mis - all$normal2_mis 
	#Coverage filtered
	cov_filter <- subset(all, normal1_cov > 30 & normal2_cov > 30 & heat1_cov > 30 & heat2_cov > 30 )
	print(paste("Positions passing coverage threshold of 30 = ", paste(nrow(cov_filter)), "sites"))
	cov_filter_Y <- subset(cov_filter, Heat=="Yes" | Heat =="No")
	print(paste("Positions passing coverage threshold of 30 = ", paste(nrow(cov_filter_Y)), "Y sites"))
	#Alteration filtered list (Mismatch difference)
	alt_filter <- subset(cov_filter, diff_mis1 > 0.1 & diff_mis2 > 0.1)
	print(paste("Mismatch difference replicably changing (0.1) = ", paste(nrow(alt_filter)), "sites"))
	alt_filter_Y <- subset(alt_filter, Heat=="Yes" | Heat =="No")
	print(paste("Mismatch difference replicably changing (0.1) = ", paste(nrow(alt_filter_Y)), "Y sites"))
	#Y-like filtered list (Mismatch and C_freq)
	y_filter <- subset(alt_filter, heat1_mis > mis_thr & heat1_c > c_thr & heat2_mis > mis_thr & heat2_c > c_thr )
	print(paste("Final Predicted Sites = ", paste(nrow(y_filter)), "sites"))
	y_filter_y <- subset(y_filter, Heat=="Yes" | Heat =="No")
	print(paste("Final Predicted Heat Responsive Sites = ", paste(nrow(y_filter_y)), "Y sites"))
	#Y-like validated list
	print(paste("Final Predicted Validated Heat Responsive Sites = ", paste(nrow(y_filter_y)), "Y sites"))
	#FINAL CLASSIFICATION
	#Novel Predicted
	novel_predicted <- subset(y_filter ,Heat=="Unm")
	novel_predicted$New_Status <- "Novel"
	novel_predicted$Prediction <- "Predicted"
	#Predicted Known HEat
	knownwheat_predicted <- subset(y_filter ,Heat=="Yes")
	knownwheat_predicted$New_Status <- "Known Heat-Sensitive"
	knownwheat_predicted$Prediction <- "Predicted"
	#Predicted Known NonHeat
	knownwnonheat_predicted <- subset(y_filter ,Heat=="No")
	knownwnonheat_predicted$New_Status <- "Known NonHeat-Sensitive"
	knownwnonheat_predicted$Prediction <- "Predicted"
	#NonPredicted Known Heat
	heat_known <- subset(cov_filter,Heat=="Yes")
	heat_known$New_Status <- "Known Heat-Sensitive"
	heat_known$Prediction <- "Not Predicted"
	#NonPredicted Known NonHeat
	nonheat_known <- subset(cov_filter,Heat=="No")
	nonheat_known$New_Status <- "Known NonHeat-Sensitive"
	nonheat_known$Prediction <- "Not Predicted"
	#Bind all the positions, remove the duplicated ones from the table
	all_newstatus <- rbind(novel_predicted, knownwheat_predicted, knownwnonheat_predicted, heat_known, nonheat_known)
	all_newstatus2 <- all_newstatus[!duplicated(all_newstatus$chr_pos.x), ]
	all_newstatus2$diff_mis1 <- all_newstatus2$heat1_mis - all_newstatus2$normal1_mis 
	all_newstatus2$diff_mis2 <- all_newstatus2$heat2_mis - all_newstatus2$normal2_mis 
	all_newstatus2$mis_normal_merged <- rowMeans(all_newstatus2[,c("normal1_mis", "normal2_mis")])
	all_newstatus2$mis_heat_merged <- rowMeans(all_newstatus2[,c("heat1_mis", "heat2_mis")])
	write.table(all_newstatus2, file="predictions_ncRNA_normal_heat.tsv", sep="\t", quote=FALSE, row.names=FALSE)
	return(all_newstatus2)
}

normal_vs_heat_predictions <- predict_heat(Normal_ncRNA_REP1_bases, Heat_ncRNA_REP1_bases, Normal_ncRNA_REP2_bases, Heat_ncRNA_REP2_bases)




## COLOR SCHEME
color_known <- "#c9c9c9"
color_known_heat <- "#ff9a8c"
color_prediction <- "#53A045"

## Categories
known <- subset(normal_vs_heat_predictions, Heat== "No" | Heat == "Yes")
known_heat <- subset(normal_vs_heat_predictions, Heat== "Yes")
all_predictions <- subset(normal_vs_heat_predictions,Prediction =="Predicted")



List_Known_vs_Predicted_Heat  <- list(known[,3], known_heat[,3], all_predictions[,3])
#Reduce(intersect, List_Known_vs_Predicted_Heat)

venn.diagram(x =List_Known_vs_Predicted_Heat,
	 category.names = c("Reported (Schwartz&Carlile)" , "Reported Heat Sensitive(Schwartz)", "Predicted Heat Sensitive"),
	 filename = 'Predicted_KNOWN_ncRNA_VENN.tiff',
	 lwd = 0.1,
	 #lty = 'blank',
	 fill = c(color_known,  color_known_heat, color_prediction),
	 output=TRUE,
     # Output features
    imagetype="tiff" ,
    height = 6000 , 
    width = 6000 , 
    resolution = 6000,
    compression = "lzw",
   # Numbers
    cex = .5,
    fontface = "bold",
    fontfamily = "sans",
    # Set names
    cat.cex = 0.1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 25, 25 ),
    cat.col = c(color_known,  color_known_heat, color_prediction),
    #cat.dist = c(0.15, 0.15 ,-0.5),
    cat.fontfamily = "sans",
    #rotation = 1
)