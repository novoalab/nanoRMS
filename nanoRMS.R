#!/usr/env/Rscript 

#####################################################
##         RNA modification stoichiometry          ##
##  using direct RNA nanopore sequencing (NanoRMS) ##
#####################################################
## Epitranscriptomics and RNA Dynamics Lab  #########
## Center for Genomic Regulation (CRG)      #########
## License: MIT                             #########
#####################################################


# Library requirements
library(stringr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(ggExtra)
library(caTools)

# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
	stop("At least two arguments must be supplied.\n--> Usage: R --vanilla < nanoRMS.R --args file1.tsv file2.tsv \n--> An optional third argument can be provided, being the clustering method to use: kmeans, kmeans_pca or knn\n", call.=FALSE)
}

########################################
# 1. READ DATA AND CLEAN UP 
########################################

# Test Read files
#args=c("LSR1_snRNA_35_normal_15mer.perread.tsv","LSR1_snRNA_35_warm_15mer.perread.tsv")


### Functions

# Clean files
clean_input<-function(input,label) {
	input$sample<- rep(label, nrow(input))
	input$read<- as.character(input$read_index)
	input<- input[!duplicated(input$read),]
	input_naomit<- na.omit(input)
	return(input_naomit)
}

# Merge files 
merge_with_wt_4samples<-function(wt,ko1,ko2,ko3) {
	merged<- rbind(wt,ko1,ko2,ko3)
	return(merged)
}

merge_with_wt_2samples<-function(wt,ko1) {
	merged<- rbind(wt,ko1)
	return(merged)
}

### Run

colnames<-c("unique","Strain","read_index",paste("pos",rep(1:15),sep=""))

if (length(args)==2) {
	wt_input<- read.delim(args[1], col.names=colnames)
	ko1_input<- read.delim(args[2], col.names=colnames)
	
	wt<-clean_input(wt_input,"sample1")
	ko1<-clean_input(ko1_input,"sample2")
	
	dat<-merge_with_wt_2samples(wt,ko1)
	
} else if (length(args)==3) {
	wt_input<- read.delim(args[1], col.names=colnames)
	ko1_input<- read.delim(args[2], col.names=colnames)
	method<-args[3]
	
	wt<-clean_input(wt_input,"sample1")
	ko1<-clean_input(ko1_input,"sample2")
	
	dat<-merge_with_wt_2samples(wt,ko1)
} else {
	stop("At least two arguments must be supplied.\n--> Usage: R --vanilla < nanoRMS.R --args file1.tsv file2.tsv \n--> An optional third argument can be provided, being the clustering method to use: kmeans, kmeans_pca or knn\n", call.=FALSE)
}


########################################
# 2. STOICHIOMETRY 
########################################
	
# a) K-means PCAed data
########################

do_kmeans_pca<-function(dat) {
	x<-scale(as.matrix(dat[,4:18]), center=TRUE) # data without labels (only 15 positions)
	pca<-prcomp(x)	
	
	pcadat<-cbind(as.data.frame(pca$x),dat$sample)
	colnames(pcadat)<-c(paste("PC",rep(1:15),sep=""),"sample")
	require(ggplot2)
	ggplot(pcadat, aes(PC1,PC2,color=sample))+geom_point(alpha=0.1)
	
	kmeans_fit<-kmeans(pca$x[,1:2],2)
	#kmeans_fit<-kmeans(pca$x,2)

	kmeans_clusters<-cbind(as.data.frame(kmeans_fit$cluster),dat$sample)
	print(table(kmeans_clusters))
	
	pcadat2<-cbind(as.data.frame(pca$x),dat$sample, kmeans_fit$cluster)
	colnames(pcadat2)<-c(paste("PC",rep(1:15),sep=""),"sample","cluster")
	ggplot(pcadat2, aes(PC1,PC2,color=cluster))+geom_point(alpha=0.1)
	
	return(kmeans_fit)
}
		
# b) K-means SCALED data
#########################

do_kmeans<-function(dat) {
	x<-scale(as.matrix(dat[,4:18]), center=TRUE) # data without labels (only 15 positions)
	
	pcadat<-cbind(as.data.frame(x),dat$sample)
	colnames(pcadat)<-c(paste("PC",rep(1:15),sep=""),"sample")
	require(ggplot2)
	ggplot(pcadat, aes(PC1,PC2,color=sample))+geom_point(alpha=0.1)
	
	kmeans_fit<-kmeans(x,2)
	
	kmeans_clusters<-cbind(as.data.frame(kmeans_fit$cluster),dat$sample)
	print(table(kmeans_clusters))
	
	pcadat2<-cbind(as.data.frame(x),dat$sample, kmeans_fit$cluster)
	colnames(pcadat2)<-c(paste("PC",rep(1:15),sep=""),"sample","cluster")
	ggplot(pcadat2, aes(PC1,PC2,color=cluster))+geom_point(alpha=0.1)
	
	return(kmeans_fit)
}
		

# c) KNN 2 samples
#########################

subdivide_training_testing<-function(dat, type, split_ratio) {
	# Get reproducible "sampling"	
	set.seed(101) 
	print(head(dat))
	sample_case = sample.split(dat$sample, SplitRatio = split_ratio)	
	
	# Test or train
	if (type=="train") {
		train = subset(dat, sample_case == TRUE)
		print (dim(train))
		return(train)
		
	}
	if (type=="test") {
		test  = subset(dat, sample_case == FALSE)		
		print(dim(test))
		return(test)
	}	
}

get_knn<-function(dat.2129,label_vector,x) { 
	
	# subdivide into train and test
	dat<-dat.2129[dat.2129$sample %in% label_vector,]
	train<-subdivide_training_testing(dat,"train",0.5) # 50% training
	test<-subdivide_training_testing(dat,"test",0.5)  # 50% testing
	
	# knn
	library(class)
	knn <- knn(train=train[,4:18], test=test[,4:18], cl=as.factor(train$sample), k=x)
	print(table(knn,test$sample))
	#accuracy<-mean(knn == test$sample)
	#print(paste("Accuracy:", accuracy))
	return(knn)
}


# d) KNN 4 samples (2 for training, 2 for testing)
####################################################

get_knn_4samples<-function(dat.2129,label_vector,x) {
	
	# subdivide into train and test
	dat<-dat.2129[dat.2129$sample %in% label_vector,]
	train<-subdivide_training_testing(dat,"train",0.5) # 50% training
	test_part1<-subdivide_training_testing(dat,"test",0.5)  # 50% testing
	
	'%!in%' <- function(x,y)!('%in%'(x,y))
	
	test_part2<-dat.2129[dat.2129$Strain %!in% label_vector,]
	test<-rbind(test_part1,test_part2)

	# knn
	library(class)
	knn <- knn(train=train[,4:18], test=test[,4:18], cl=as.factor(train$sample), k=x)
	print(table(knn,test$sample))
	#accuracy<-mean(knn == test$sample)
	#print(paste("Accuracy:", accuracy))
	return(knn)
}


########################################
# 3. RUN 
########################################

## Choose default method
if (exists("method")=="FALSE") {
	method="kmeans"
	#print ("Using kmeans")
}

## Run prediction based on choice:
set.seed(101)
if (method=="kmeans_pca") {	
	dat.kmeans_pca<-do_kmeans_pca(dat)
} else if  (method=="kmeans") {	
	dat.kmeans<-do_kmeans(dat)
} else if  (method=="knn") {	
	knn_model<-get_knn(dat,c("sample1","sample2"),7)
#} else if (method=="knn_4samples") {	
#	knn_model_4samples<-get_knn_4samples(dat,c("wt","ko1"),7)
} else {
	print ("ERROR: Please choose among the following methods: kmeans_pca, kmeans, knn")
}
# Have a nice day! :)

