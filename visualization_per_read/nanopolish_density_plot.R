#Scripts how to use
#Rscript --vanilla density_nanopolish.R data1 data2 data3 data4

#Libraries needed
library(ggplot2)

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

merged_2 <- subset(merged, reference=="0")

for (pos in unique(merged_2$Pos)) {
	subs<- subset(merged_2, Pos==pos)
	pdf(file=paste(pos, "density_plot.pdf", sep="_"),height=4,width=9,onefile=FALSE)
		print(ggplot(subs, aes(x= event_level_mean, fill=sample)) +
 		geom_density(alpha=0.3,adjust = 2)+
  		theme_bw()+
		ggtitle(paste(pos))+
		xlab("Current Intensity (pA) ")+
      	ylab("Density") +
		theme(axis.text=element_text(size=15),strip.text = element_text(size=13),
    		axis.title=element_text(size=17,face="bold"),
    		legend.title = element_text(size = 20),
    		plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
    		legend.text = element_text(color = "black", size=15)))
	dev.off()
}