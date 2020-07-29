#Scripts how to use
#Rscript --vanilla nanopolish_export_each_position.R data1 data2 data3 data4

#Libraries needed
library(dplyr)
library(ggplot2)
library(plyr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
# test the number of arguments
if (length(args) < 1) {
  stop("At least three arguments must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
   input1 <- args[1]
   len <- 1
} else if (length(args)==2) {
   input1 <- args[1]
   input2 <- args[2]
   len <- 2
} else if (length(args)==3) {
   input1 <- args[1]
   input2 <- args[2]
   input3 <- args[3]
   len <- 3
} else if (length(args)==4) {
   input1 <- args[1]
   input2 <- args[2]
   input3 <- args[3]
   input4 <- args[4]
   len <- 4
} else {
  stop("Incorrect number of arguments", call.=FALSE)
}




if (len == 1) {
	data1<- read.delim(input1)
	merged <- data1
} else if (len == 2) {
   	data1<- read.delim(input1)
	data2<- read.delim(input2)
	merged <- rbind(data1, data2)
} else if (len == 3) {
   	data1<- read.delim(input1)
	data2<- read.delim(input2)
	data3<- read.delim(input3)
	merged <- rbind(data1, data2, data3)
} else if (len == 4) {
   	data1<- read.delim(input1)
	data2<- read.delim(input2)
	data3<- read.delim(input3)
	data4<- read.delim(input4)
	merged <- rbind(data1, data2, data3, data4)
}


merged2<- merged[,c("modification","sample","read_index", "reference", "event_level_mean" )]
windows_casted<- dcast(merged2, modification+sample+read_index ~ reference )
windows_naomit<- na.omit(windows_casted)
colnames(windows_naomit)<- c("unique","Strain","read_index", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1", "2","3","4","5","6", "7")



for (strain in unique(windows_naomit$Strain)) {
  subs_strain <- subset(windows_naomit, Strain==strain)
    for (uniq in unique(windows_naomit$unique)) {
      subs_uniq<- subset(subs_strain, unique==uniq)
      write.table(subs_uniq, file=paste(uniq, strain, "15mer.perread.tsv", sep="_"), quote=FALSE, sep="\t")
    }
  }
