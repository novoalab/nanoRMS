#####################################################################
######## SCRIPTS FOR PREDICTING Y FROM SINGLE SAMPLE  ###############
#####################################################################
###################### BY OGUZHAN BEGIK #############################
#####################################################################
#Rscript predict_singleSample.R wt_epinano.csv sn3_epinano.csv sn36_epinano.csv


## Loading libraries
library(stringr)
library(dplyr)
library(plyr)
library(VennDiagram)
library(RColorBrewer)
library(data.table)


# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)
#Arguments
input1 <- args[1]  #1st variable
input2 <- args[2]  #2nd variable
input3 <- args[3]  #3rd variable





############################################
### PART1 : Importing the files ############
############################################
### Mod File
mod_rRNA <- read.delim("yeast_all_rrna_mod_status.tsv")
### Importing the data



rRNA_REP1 <- read.delim(input1 ,sep=",")
rRNA_REP2 <- read.delim(input2 ,sep=",")
rRNA_REP3 <- read.delim(input3 ,sep=",")




### FILTERS DETERMINED
mis_thr <- 0.137
c_thr <- 0.578
#####


############################################
### PART2#### #PROCESS rRNA DATA
############################################
prediction <- function(data,label) {
	#Only Curlcake 3
	#Coverage filter
	data<- subset(data, cov > 30)
	data$chr_pos <- paste(data$X.Ref, data$pos, sep="_")
	bases <- str_split_fixed(data$ACGT_freq, n=4, pattern=":")
	colnames(bases) <- c("A", "C", "G", "T")
	data2 <- cbind(data, bases)
	data3 <- data2[,c("X.Ref", "pos", "chr_pos", "base", "cov", "q_mean", "q_median", "q_std", "ins", "del", "mis", "A", "T", "C", "G")]
	data3$A <- as.numeric(as.character(data3$A))
	data3$T <- as.numeric(as.character(data3$T))
	data3$G <- as.numeric(as.character(data3$G))
	data3$C <- as.numeric(as.character(data3$C))
	#Calculate mismatches
	mismatches <- c()  
	for (i in c(1:nrow(data3))){   
		base <- data3[i,c("base")]
		a <- sum(data3[i,c("A","T", "C", "G")])-data3[i,toString(base)]
		mismatches <- c(mismatches, a)
	}
	data3 <- cbind(data3, mismatches)
 	data3$count <- data3$A+ data3$C + data3$G + data3$T
	data3$mis_freq <- data3$mismatches/data3$count
	#Position filter
	final<- vector()
	for (chro in unique(data3$X.Ref)) {
		subs<- subset(data3, X.Ref==chro)
		subs2<- subset(subs, pos >30)
		subs3<- subset(subs2, pos <(max(subs$pos)-30))
		final<- rbind(final, subs3)
	}
	final_U <- subset(final, base=="T")
	final_U$A_freq <- final_U$A / final_U$mismatches
	final_U$C_freq <- final_U$C / final_U$mismatches
	final_U$G_freq <- final_U$G / final_U$mismatches
	final_U <- final_U %>% mutate(A_freq = coalesce(A_freq, 0))
	final_U <- final_U %>% mutate(C_freq = coalesce(C_freq, 0))
	final_U <- final_U %>% mutate(G_freq = coalesce(G_freq, 0))
	final_U$sample <- label
	final_U2 <- merge(final_U, mod_rRNA, by.x=c("X.Ref", "pos"), by.y=c("Chr", "Position"))
	final_U3 <- final_U2[,c("X.Ref", "pos", "chr_pos", "base","sample", "ModStatus", "Status", "cov", "q_mean", "q_median", "q_std", "ins", "del", "mis", "A", "T", "C", "G","mis_freq","A_freq","C_freq","G_freq")]
	final_U4 <- subset(final_U3, Status =="Unm" | Status =="Y")
	mitochondrial_rrna <- subset(final_U4, X.Ref == "15s" | X.Ref == "21s" )
	mitochondrial_rrna$Status <- factor(mitochondrial_rrna$Status)
	data_mod <- subset(mitochondrial_rrna, mis > mis_thr & C_freq >c_thr)
	data_mod$Prediction <- "Mod"
	data_unmod <- subset(mitochondrial_rrna, mis < mis_thr)
	data_unmod$Prediction <- "Unmod"
	data_final<- rbind(data_mod, data_unmod)
	table <- data.frame(table(data_final$Status,data_final$Prediction))
	colnames(table) <- c("Observation", "Prediction" ,"Freq")
	print(table)
	return(data_final)
}

rep1_predictions <- prediction(rRNA_REP1, "Rep1")
rep2_predictions <- prediction(rRNA_REP2, "Rep2")
rep3_predictions <- prediction(rRNA_REP3, "Rep3")

# Positions predicted as MODIFIED
rep1_pred_mod <- subset(rep1_predictions, Prediction =="Mod")
rep2_pred_mod <- subset(rep2_predictions, Prediction =="Mod")
rep3_pred_mod <- subset(rep3_predictions, Prediction =="Mod")

# REPLICABLE ONES
List_all <- list(rep1_pred_mod[,3], rep2_pred_mod[,3], rep3_pred_mod[,3])
reproducible_sites <-  Reduce(intersect, List_all)
write.table(reproducible_sites,"Yeast_MitochondrialrRNA_reproducible_predicted_sites.tsv", sep="\t", quote=FALSE, row.names=FALSE)

#VENN DIAGRAM
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(x =List_all,
	 category.names = c("Mit-rRNA Rep 1" , "Mit-rRNA Rep 2 " , "Mit-rRNA Rep 3"),
	 filename = 'WT_venn_diagramm_mitochondrial.tiff',
	 lwd = 2,
	 lty = 'blank',
	 fill = myCol,
	 output=TRUE,
     # Output features
    imagetype="tiff" ,
    height = 500 , 
    width = 500 , 
    resolution = 500,
    compression = "lzw",
   # Numbers
    cex = .3,
    fontface = "bold",
    fontfamily = "sans",
    # Set names
    cat.cex = 0.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
)