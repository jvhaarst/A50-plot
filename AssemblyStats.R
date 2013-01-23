# R script to:
# Read in multiple genome assemblies
# Calculate metrics on scaffolds
# Generate table of metrics of all assemblies
# Generate graph to show assembly growth

max_count <- 20000

# Load the needed R libraries (from Bioconductor)
require("Biostrings")
require("IRanges")

# Create list to hold assembly info from csv file
readSeq<-function(filename){
  sortedSizes <- sort(width(readDNAStringSet(filename)), decreasing=TRUE)
  return(sortedSizes)
}

input<-read.csv("input.csv", as.is=TRUE, header=TRUE, blank.lines.skip=TRUE, comment.char="#")
attach(input)
inputCols<-ncol(input)
inputRows<-nrow(input)
inputValid<-0

N <- list()

for (rowNum in 0:inputRows) {
  seqName<-Seq_Name[ rowNum ]
  seqFile<-Seq_File[ rowNum ]
    
  if ( length(seqName) > 0 && nchar(seqName) > 0 ) {
    print(paste("Row", rowNum))
    print(paste("  seq name: '", seqName, "'", sep=""))
    print(paste("  seq file: '", seqFile, "'", sep=""))

    N[[ seqName ]] <- readSeq(seqFile)

    inputValid<-inputValid + 1
  }
}

print(paste("TOTAL VALID SEQUENCES: ", inputValid, sep=""))

# Get the maximum contig count from the list
#max_count <- max(unlist(lapply(N,length)))

# Read external script for functions
source('./contigStats.R')

# Create plot and statistics
reflength <- sapply(N, sum)
print(reflength)
stats<-contigStatsFlipped(style="data",N=N, reflength=reflength)
print(stats)
cat("\"Name\"\t" , file="Rplots_stats.csv",           append=FALSE)
write.table(stats, file="Rplots_stats.csv", sep="\t", append=TRUE )
contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlim=c(0,max_count),
                    trimSize=max_count,
                    xlab="Number of contigs",
                    ylab="Cumulative contig length",
                    main="Cumulative Plot of N Statistic"
)
