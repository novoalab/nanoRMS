#####################################################################
###  SCRIPTS FOR PREDICTING Y FROM PAIRED SAMPLE (Normal vs Heat)  ##
###################### BY OGUZHAN BEGIK #############################
#####################################################################
#WT VS KO, NORMAL VS STRESS
## Loading libraries
library(optparse)

parser <- OptionParser(usage = "%prog [options] -f <epinano_file1> -s <epinano_file2> ")
parser<- add_option(parser, c("-f", "--epinano1"), type="character", default=NULL, 
              help="first epinano file", metavar="character")
parser<- add_option(parser,c("-s", "--epinano2"), type="character", default=NULL, 
          help="second epinano file", metavar="character")
parser<- add_option(parser,c("-p", "--modpos"), type="character", default="positions/RNA_Mod_Positions_mRNAYeast_HeatSensitive.tsv", 
	      help="mod positions file [default= %default]", metavar="character")
parser <- add_option(parser, c("-m", "--misfreq"), type="numeric", default=0.137,
                help="Mismatch frequency threshold [default= 0.137]",
                metavar="number")
parser <- add_option(parser, c("-c", "--Cfreq"), type="numeric", default=0.578,
                help="C mismatch frequency threshold [default= 0.578]",
                metavar="number")
parser <- add_option(parser, c("-d", "--diff"), type="numeric", default=0.1,
                help="mismatch frequency difference threshold [default= 0.1]",
                metavar="number")

opt = parse_args(parser);



## Loading libraries
library(stringr)
library(dplyr)
library(plyr)
library(VennDiagram)
library(RColorBrewer)
library(data.table)

##FUNCTION DEFINED
processing <- function(data) {
	#Import the mod file
	mod_file <- read.delim(opt$modpos, sep="") #opt$modpos Mod positions file RNA_Mod_Positions_mRNAYeast_HeatSensitive.tsv
	#seperate one column into multiple columns
	columns <- str_split_fixed(data$V5, n=12,  pattern=",")
	#Add these columns to the original table
	data2 <- cbind(data, columns)
	#Select important columns
	data3 <- data2 [,c("1", "2","3","4", "5","V6","V7", "V8","V9", "V10" ,"6","7", "8","9", "10", "11","12")]
	#Rename the columns
	colnames(data3) <- c("Chr", "Pos", "Base", "Strand", "Coverage", "cds_chr", "cds_start", "cds_end", "cds_strand", "cds_name", "Q_Mean", "Q_Median", "Q_STD", "Mis", "Ins", "Del", "ACGT_Freq")
	#seperate one column into multple columns
	base_freq <-  str_split_fixed(data3$ACGT_Freq, n=4, pattern=":")
	#Rename columns
	colnames(base_freq) <- c("A", "C", "G","T")
	#Add these columns into the original table
	data4<- cbind(data3, base_freq)
	#Remove this column
	data4$ACGT_freq <- NULL
	#Make the counts numeric
	data4$A <- as.numeric(as.character(data4$A))
	data4$T <- as.numeric(as.character(data4$T))
	data4$G <- as.numeric(as.character(data4$G))
	data4$C <- as.numeric(as.character(data4$C))
	data4$Mis <- as.numeric(as.character(data4$Mis))
	data4$Pos <- as.numeric(as.character(data4$Pos))
	data4$Coverage <- as.numeric(as.character(data4$Coverage))
	data4$cds_start <- as.numeric(as.character(data4$cds_start))
	data4$cds_end <- as.numeric(as.character(data4$cds_end))
	#Coverage filter
	data5<- subset(data4, Coverage > 30)
	#Create a column with unique positions
	data5$Chr_Pos <- paste(data5$Chr, data5$Pos, sep="_")
	#Select important columns for analysis
	data6 <- data5[,c("Chr", "Pos", "Chr_Pos", "Base","Strand","cds_chr", "cds_start", "cds_end", "cds_strand","cds_name", "Coverage", "Mis", "A", "T", "C", "G")]
	#Treat each strand seperately
	#Positive Strand
	data_subset_positive <- subset(data6 , Strand=="+")
	data_subset_positive <- subset(data_subset_positive, Base=="T")
	data_subset_positive$rel_pos <- data_subset_positive$Pos-data_subset_positive$cds_start
	data_subset_positive$rel_pos <- data_subset_positive$rel_pos+1
	data_subset_positive$mismatches<- data_subset_positive$A+ data_subset_positive$C + data_subset_positive$G
	data_subset_positive$C_freq <-data_subset_positive$C/data_subset_positive$mismatches
	data_subset_positive <- data_subset_positive %>% mutate(C_freq = coalesce(C_freq, 0))
	data_subset_positive <- data_subset_positive[,c("Chr", "Pos", "Chr_Pos", "Base","Strand","cds_chr", "cds_start", "cds_end", "cds_strand","cds_name", "Coverage", "Mis", "C_freq")]
	#Negative Strand
	data_subset_negative <- subset(data6 , Strand=="-")
	data_subset_negative <- subset(data_subset_negative, Base=="A")
	data_subset_negative$rel_pos<- data_subset_negative$cds_end - data_subset_negative$Pos
	data_subset_negative$rel_pos<- data_subset_negative$rel_pos+1
	data_subset_negative$mismatches<- data_subset_negative$T+ data_subset_negative$C + data_subset_negative$G
	data_subset_negative$C_freq <-data_subset_negative$G/data_subset_negative$mismatches
	data_subset_negative <- data_subset_negative %>% mutate(C_freq = coalesce(C_freq, 0))
	data_subset_negative <- data_subset_negative[,c("Chr", "Pos", "Chr_Pos", "Base","Strand","cds_chr", "cds_start", "cds_end", "cds_strand","cds_name", "Coverage", "Mis", "C_freq")]
	#Merge both
	data_final <- rbind(data_subset_positive,data_subset_negative )
	data_final_mods <- merge(data_final, mod_file, all.x=TRUE, by.y=c("Chr", "Genome_Pos", "Strand"), by.x=c("Chr", "Pos", "Strand"))
	data_final_mods$Pseudouridine <- as.character(data_final_mods$Pseudouridine)
	data_final_mods$Pseudouridine[is.na(data_final_mods$Pseudouridine)] <- "No"
	data_final_mods$Gene_Name <- as.character(data_final_mods$Gene_Name)
	data_final_mods$Gene_Name[is.na(data_final_mods$Gene_Name)] <- "No"
	data_final_mods$Gene_Pos <- as.character(data_final_mods$Gene_Pos)
	data_final_mods$Gene_Pos[is.na(data_final_mods$Gene_Pos)] <- "No"
	return(data_final_mods)
}


difference <- function(data1, data2) {
	processed1 <- processing(data1)
	processed2 <- processing(data2)
	mis_thr <- opt$misfreq #
	c_thr <- opt$Cfreq #
	#PseudoU-like sites 
	#For the first sample
	pu_like1 <-  subset(processed1, Mis > mis_thr & C_freq> c_thr)
	pu_like2 <-  subset(processed2, Mis > mis_thr & C_freq> c_thr)
	#Colnames vector to use
	columns <- c("Chr", "Pos", "Chr_Pos", "Base", "Pseudouridine")
	pu_like_sites_merged<- rbind(pu_like1[,columns], pu_like2[, columns])
	pu_like_sites_merged<- pu_like_sites_merged[!duplicated(pu_like_sites_merged$Chr_Pos), ]
	merged1 <- merge(pu_like_sites_merged, processed1, by.x=columns, by.y=columns, all.x=TRUE)
	merged2 <- merge(merged1, processed2, by.x=columns, by.y=columns, suffixes = c(".data1",".data2"))
	#Remove NA 
	merged3 <- na.omit(merged2)
	#Taking the difference between mismatch frequency
	merged3$Mis.difference <- merged3$Mis.data1 -merged3$Mis.data2
	diff_thr <- opt$diff
	merged_different <- subset(merged3, abs(Mis.difference) > diff_thr)
	write.table(merged_different, file="Paired_comparison_altering_sites_pU_predictions.tsv", quote=FALSE, sep="\t", row.names=FALSE)
}





if (!is.null(opt$epinano1) && is.null(opt$epinano2)) {
	print_help(parser)
	stop("At least two epinano input must be supplied (input file).n", call.=FALSE)
} else if (!is.null(opt$epinano1) && !is.null(opt$epinano2)) {
	#Importing files
	data_rep1 <- read.delim(opt$epinano1,header=FALSE)
	data_rep2 <- read.delim(opt$epinano2,header=FALSE)
	difference(data_rep1, data_rep2)
} 
