#Scripts how to use
# Rscript --vanilla nanopolish_window.R positions_file input_table label

args = commandArgs(trailingOnly=TRUE)
argument1<- args[1]
argument2<- args[2]
argument3<- args[3]

# test the number of arguments
if (length(args) < 3) {
  stop("At least three arguments must be supplied (input file).n", call.=FALSE)
} else if (length(args)==3) {
   input1 <- args[2]
   label1<- args[3]
} else {
  stop("Incorrect number of arguments", call.=FALSE)
}


# Import the positions to be looked into
positions <- read.delim(argument1)
input1 <- read.csv(argument2, sep="\t", header=TRUE)
label1<- argument3

# Process the data input
process <- function(data, label) { 
  data$position <- data$position+3 # Add 3 nt to each position 
  data$Pos<- paste(data$contig, data$position, sep="_") # Unique column
  data$sample <- label #Add label
  windows <- vector() # Create an empty vector
  for (rown in 1:nrow(positions)){ #Create windows file
    chr <- positions[rown, 1]
    subs <- subset(data, contig==as.character(chr))
    mod <- positions[rown, 2]
    window<- c(mod-7,mod-6, mod-5, mod-4, mod-3, mod-2, mod-1,mod, mod+1, mod+2, mod+3, mod+4,mod+5, mod+6, mod+7)
      for(wind in window[1:length(window)]) {
        subs2<- subset(subs, position==wind)
        subs2$modification<- paste(chr, mod, sep="_")
        subs2$reference<- paste(wind-mod)
        windows<-rbind(windows, subs2)
      }
  }
return(windows)
}

final_window <- process(input1, label1)

write.table(final_window, file=paste(label1, "window_file.tsv", sep="_"), quote=FALSE, sep="\t")
