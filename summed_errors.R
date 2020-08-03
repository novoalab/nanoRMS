#!/usr/bin/env Rscript

#####################################################
##         RNA modification stoichiometry          ##
##  using direct RNA nanopore sequencing (NanoRMS) ##
#####################################################
## Epitranscriptomics and RNA Dynamics Lab  #########
## Center for Genomic Regulation (CRG)      #########
## License: MIT                             #########
## Author: Oguzhan Begik	                #########
#####################################################


# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Epinano output
input1 <- args[1] #1st variable


#Read input
data <- read.delim(input1,sep=",")
#add the sum column
data$sum <- data$mis + data$del + data$ins 
#Export the new table
write.table(data, file="output.txt", quote=FALSE, sep=",")








