#####################################################################
### SCRIPTS FOR PREDICTING Y FROM SINGLE SAMPLE (Reproducibility)####
#####################################################################
###################### BY OGUZHAN BEGIK #############################
#####################################################################

## Loading libraries
library(optparse)



parser <- OptionParser(usage = "%prog [options] -f <epinano_file1> (-s <epinano_file2> -t <epinano_file3>)")
parser<- add_option(parser, c("-f", "--epinano1"), type="character", default=NULL, 
              help="first epinano file", metavar="character")
parser<- add_option(parser,c("-s", "--epinano2"), type="character", default=NULL, 
          help="second epinano file", metavar="character")
parser<- add_option(parser,c("-t", "--epinano3"), type="character", default=NULL, 
          help="third epinano file", metavar="character")
parser<- add_option(parser,c("-p", "--modpos"), type="character", default="positions/RNA_Mod_Positions_rRNAYeast.tsv", 
	      help="mod positions file [default= %default]", metavar="character")
parser <- add_option(parser, c("-m", "--misfreq"), type="numeric", default=0.137,
                help="Mismatch frequency threshold [default= 0.137]",
                metavar="number")
parser <- add_option(parser, c("-c", "--Cfreq"), type="numeric", default=0.578,
                help="C mismatch frequency threshold [default= 0.578]",
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
prediction <- function(data,label) {
	#Import the mod file
	mod_file <- read.delim(opt$modpos) #Mod positions file #RNA_Mod_Positions_rRNAYeast.tsv
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
		subs2<- subset(subs, Pos >30)
		subs3<- subset(subs2, Pos <(max(subs$Pos)-30))
		final<- rbind(final, subs3)
	}
	#Label the data
	final$sample <- label
	#Merge the data with mod positions
	final2 <- merge(final, mod_file, by.x=c("Chr", "Pos"), by.y=c("Chr", "Position"))
	#Reorder the columns
	final3 <- final2[,c("Chr", "Pos", "Chr_Pos", "Base","sample", "Mod","Neighbour",  "Coverage", "Mis","C_freq")]
	#Select only UNMODIFIED positions NOT neighbouring another modification
	final4 <- subset(final3, Mod =="Unm" & Neighbour=="No")
	#Remove other Levels
	final4$Mod <- factor(final4$Mod)
	#Passing the threshold
	mis_thr <- opt$misfreq #
	c_thr <- opt$Cfreq #
	data_mod <- subset(final4, Mis > mis_thr & C_freq > c_thr)
	return(data_mod)
}

one_input <- function(data,label) {
	predicted <- prediction(data, label)
	write.table(predicted, file=paste(label, "predicted_y_sites.tsv", sep="_"), sep="\t", quote=FALSE, row.names=FALSE)
}


two_input <- function(data1,label1,  data2, label2) {
	predicted1 <- prediction(data1, label1)
	predicted2 <- prediction(data2, label2)
	List_all <- list(predicted1[,3], predicted2[,3])
	reproducible_sites <-  Reduce(intersect, List_all)
	write.table(reproducible_sites,file="twosample_predicted_reproducible_y_sites.tsv", sep="\t", quote=FALSE, row.names=FALSE)
	#VENN DIAGRAM
	venn.diagram(x =List_all,
		 category.names = c(label1 , label2),
		 filename = 'Twosamples_venn_diagramm_mitochondrial.tiff',
		 lwd = 2,
		 lty = 'blank',
		 fill = c("red", "blue"),
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
	    #cat.pos = c(-27, 27),
	    #cat.dist = c(0.055, 0.055),
	    #cat.fontfamily = "sans",
	    #rotation = 1
	)
}


three_input <- function(data1,label1,  data2, label2, data3, label3) {
	predicted1 <- prediction(data1, label1)
	predicted2 <- prediction(data2, label2)
	predicted3 <- prediction(data3, label3)
	List_all <- list(predicted1[,3], predicted2[,3], predicted3[,3])
	reproducible_sites <-  Reduce(intersect, List_all)
	write.table(reproducible_sites,file="threesample_predicted_reproducible_y_sites.tsv", sep="\t", quote=FALSE, row.names=FALSE)
	#VENN DIAGRAM
	myCol <- brewer.pal(3, "Pastel2")
	venn.diagram(x =List_all,
		 category.names = c(label1 , label2, label3),
		 filename = 'Threesamples_venn_diagramm_mitochondrial.tiff',
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
}



if (is.null(opt$epinano1)) {
	print_help(parser)
	stop("At least one epinano input must be supplied (input file).n", call.=FALSE)
} else if (!is.null(opt$epinano1) && is.null(opt$epinano2) && is.null(opt$epinano3)) {
	#Importing files
	data_rep1 <- read.delim(opt$epinano1 ,sep=",")
	one_input(data_rep1, "Sample1")
} else if (!is.null(opt$epinano1) && !is.null(opt$epinano2) && is.null(opt$epinano3)) {
	#Importing files
	data_rep1 <- read.delim(opt$epinano1 ,sep=",")
	data_rep2 <- read.delim(opt$epinano2,sep=",")
	two_input(data_rep1, "Sample1", data_rep2, "Sample2")
} else if (!is.null(opt$epinano1) && !is.null(opt$epinano2) && !is.null(opt$epinano3)) {
	#Importing files
	data_rep1 <- read.delim(opt$epinano1 ,sep=",")
	data_rep2 <- read.delim(opt$epinano2,sep=",")
	data_rep3 <- read.delim(opt$epinano3,sep=",")
	three_input(data_rep1, "Sample1", data_rep2, "Sample2", data_rep3, "Sample3")
}


