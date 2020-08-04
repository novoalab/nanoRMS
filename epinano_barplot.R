#Delta barplot 
#Rscript epinano_barplot.R input1 label1 input2 label2 feature
#Libraries needed
library(plyr)
library(ggplot2)
library(ggrepel) 
library(MASS)
library(reshape2)

# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Arguments
input1 <- read.delim(args[1],sep=",")  #1st variable
label1 <- as.character(args[2])  #1st label
input2 <-read.delim(args[3],sep=",") #2nd variable
label2 <- as.character(args[4]) #2nd label
feature<- as.character(args[5]) #Feature



#Cleanup 
cleanup <- function(input, label) {
	#Filter low coverage reads
	input <- subset(input, cov>30)
	#Filter read starts
	input <- subset(input, pos>20)
	#For Stress
	#Add a column with position 
	input$position<- paste(input$X.Ref,input$pos)
	#Change column names 
	input <- input[, c("X.Ref","pos","position", "base", feature)]
	colnames(input)<- c("Chr","Position","chr_pos","base",feature )
	data_melted<- melt(data = input, id.vars = c("Chr", "Position", "chr_pos", "base"))
    colnames(data_melted)[which(names(data_melted) == "value")] <- paste(label, "value", sep="_")
	return(data_melted)
}


#Cleanup and process the data
data1 <- cleanup(input1, label1)
data2 <- cleanup(input2, label2)

merged <- join(data1,data2, by="chr_pos")
merged$Chr <- NULL
merged$Position <- NULL
merged$base <- NULL
merged$variable <- NULL
			

plot<- function(data)
for (chr in unique(data$Chr)) {
	subs <- subset(data,  Chr==chr)
	res<- rlm(subs[,c(paste(label1, "value", sep="_"))] ~ subs[,c(paste(label2, "value", sep="_"))]) #linear model  
	res_vec <- res$residuals#this contains residuals 
	threshold <-  5 * sd(res_vec) #The threshold
	subs$score<- abs(subs[,c(paste(label1, "value", sep="_"))] - subs[,c(paste(label2, "value", sep="_"))])
	pdf(file=paste(chr,feature, label1, label2, "barplot.pdf", sep="_"),height=5,width=20,onefile=FALSE)
	print(ggplot(subs, aes(x=Position, y=score)) +
      geom_bar(stat = "identity", width=1, fill="#2a7886") +
      geom_text_repel(data=subset(subs, score>threshold), aes(Position, score, label=Position,size=3, color="red"), segment.size  = 1,segment.color = "black")+
      ggtitle(paste(chr, feature, label1, label2, sep="_"))+
      xlab("Positions")+
      ylab("Delta Feature") +
      theme_bw()+
      theme(axis.text.x = element_text(face="bold", color="black",size=11),
        axis.text.y = element_text(face="bold", color="black", size=11),
        plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black", size=0.5)))
  dev.off()
}

plot(merged)