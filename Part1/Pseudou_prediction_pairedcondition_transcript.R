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
parser<- add_option(parser,c("-p", "--modpos"), type="character", default="RNA_Mod_Positions_ncRNAYeast_HeatSensitive.tsv", 
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
	mod_file <- read.delim(opt$modpos) #opt$modpos Mod positions file RNA_Mod_Positions_ncRNAYeast_HeatSensitive.tsv
	#Rename the columns
	colnames(data) <- c("Chr", "Pos", "Base", "Strand", "Coverage", "Q_Mean", "Q_Median", "Q_STD", "Mis", "Ins", "Del", "ACGT_Freq")
	#Coverage filter
	data<- subset(data, Coverage > 30)
	#Create a column with unique positions
	data$Chr_Pos <- paste(data$Chr, data$Pos, sep="_")
	#split the base count column
	bases <- str_split_fixed(data$ACGT_Freq, n=4, pattern=":")
	#name columns
	colnames(bases) <- c("A", "C", "G", "T")
	#Add individual base counts to the data
	data2 <- cbind(data, bases)
	#Reorder the data
	data3 <- data2[,c("Chr", "Pos", "Chr_Pos", "Base", "Coverage", "Mis", "A", "T", "C", "G")]
	#Make the counts numeric
	data3$A <- as.numeric(as.character(data3$A))
	data3$T <- as.numeric(as.character(data3$T))
	data3$G <- as.numeric(as.character(data3$G))
	data3$C <- as.numeric(as.character(data3$C))
	#Extract only U positions
	data_U <- subset(data3, Base=="T")
	#Count mismatches
	data_U$mismatches <- data_U$A + data_U$C + data_U$G
	#C mismatch frequency
	data_U$C_freq <- data_U$C/data_U$mismatches
	#Convert NA into 0
	data_U <- data_U %>% mutate(C_freq = coalesce(C_freq, 0))
	#Position filter (Remove first and last 30 nucleotides)
	final<- vector()
	for (chro in unique(data_U$Chr)) {
		subs<- subset(data_U, Chr==chro)
		subs2<- subset(subs, Pos > 14)
		subs3<- subset(subs2, Pos < (max(subs$Pos) - 14))
		final<- rbind(final, subs3)
	}
	#Label the data
	#Merge the data with mod positions
	final2 <- merge(final, mod_file, all.x=TRUE, by.x=c("Chr", "Pos"), by.y=c("Chr", "Position"))
	final2$Pseudouridine <- as.character(final2$Pseudouridine)
	final2$Pseudouridine[is.na(final2$Pseudouridine)] <- "No"
	#Reorder the columns
	final3 <- final2[,c("Chr", "Pos", "Chr_Pos", "Base", "Pseudouridine", "Mis","C_freq")]
	return(final3)
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
	data_rep1 <- read.delim(opt$epinano1 ,sep=",")
	data_rep2 <- read.delim(opt$epinano2,sep=",")
	difference(data_rep1, data_rep2)
} 
