#!/usr/bin/env Rscript

# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Epinano output
input1 <- args[1] #1st variable
input2 <- as.character(args[2]) #1st variable


#Read input
data <- read.delim(input1,sep=",")
#add the sum column
data$sum <- data$mis + data$del + data$ins 
#Export the new table
write.table(data, file=paste(input2, "epinano_withsum.csv", sep="_"), quote=FALSE, sep=",")








