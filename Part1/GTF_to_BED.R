#This script converts GTF to BED

#Libraries 
library(stringr)

#Arguments
args = commandArgs(trailingOnly=TRUE)

data <- read.delim(args[1], header=FALSE)
data2<- subset(data, V3=="CDS")
column<- str_split_fixed(data2$V9, ";",10)
gene_id <- column[,1]
gene_id2 <- str_split_fixed(gene_id, " ",2)
data2$gene_id<- gene_id2[,2]
data3 <- data2[,c("V1", "V4", "V5","V7", "gene_id")]
data3$chr<- "chr"
data3$ref<- paste(data3$chr, data3$V1, sep="")
data4 <- data3[,c("ref","V4" ,"V5", "V7" ,  "gene_id" )]
write.table(data4, file="Data.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)