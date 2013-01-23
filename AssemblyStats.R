# R script to:
# Read in multiple genome assemblies
# Calculate metrics on scaffolds
# Generate table of metrics of all assemblies
# Generate graph to show assembly growth

# max contig count. Any assembly with more than this number of contigs will be trimmed
# TODO: read from setup file
max_count <- 20000



# Load the needed R library (from Bioconductor)
require("Biostrings")
require("IRanges")

library("Biostrings")
library("IRanges")

source('./contigStats.R')



# Function to read fasta files returning a size histogram
readSeq<-function(filename){
  #sortedSizes <- sort(width(readDNAStringSet(filename)), decreasing=TRUE) # Rstudio
  sortedSizes <- sort(width(read.DNAStringSet(filename)), decreasing=TRUE) # Command line
  return(sortedSizes)
}



# Load CSV
input<-read.csv("input.csv", as.is=TRUE, header=TRUE, blank.lines.skip=TRUE, comment.char="#")
attach(input)
inputCols<-ncol(input)
inputRows<-nrow(input)
inputValid<-0



# Load Fasta
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



# Use maximal reference length for N50
print("Use maximal reference length for N50")
reflength <- sapply(N, sum)
max_ref <- as.numeric(max(reflength))

# Get table and save it
stats<-contigStatsFlipped(style="data",N=N, reflength=reflength)
print(stats)
cat("\"Name\"\t" , file="Rplots_stats.csv",           append=FALSE)
write.table(stats, file="Rplots_stats.csv", sep="\t", append=TRUE )

# Generate graphic
contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlim=c(0,max_count),
                    trimSize=max_count,
                    xlab="Number of contigs",
                    ylab="Cumulative contig length",
                    main="Cumulative Plot of A50 Statistic"
)

quit()


########### TRASH BIN ###########################



# Create list to hold assembly info
# assemblies <- list("/home/assembly/dev_150/assemblies/allpaths_lg_sample_heinz_raw/sl/data/run/ASSEMBLIES/test/final.assembly.fasta",
#                 "/home/assembly/dev_150/assemblies/allpaths_lg_sample_heinz_raw_with454/sl/data/run/ASSEMBLIES/test/final.assembly.fasta",
#                 "/home/assembly/dev_150/assemblies/clc-default/clc_contigs.fa",
#                 "/home/assembly/dev_150/assemblies/S_lycopersicum_scaffolds.2.40.fa",
#                 "/home/assembly/progs/fermi/heinz/fmdef.p4.fa"
#                 )
# assemblies <- list("velvet-SRR001665.interleaved.fasta","Galaxy143-[Contigs].fasta")
# )
# Create list of contig lenghts from assemblies
# assembly_contigs <- list()
# for(assembly in assemblies){
#   print(assembly)
#   contigs <-read.DNAStringSet(assembly, "fasta",nrec=10000)
#   assign(paste(as.name(assembly)),width(contigs))
#   assembly_contigs <- c(assembly_contigs,paste(as.name(assembly)))
#   #assembly_contigs <- c(assembly_contigs,paste(as.name(assembly))=width(contigs))
# }
# N <- list(assemblies)
# reflength <- sapply(N, sum)
# max_ref <- max(reflength)
# reflength <- sapply(N, function(x) x <- max_ref)


## Use own reference length for N50
#print("Use own reference length for N50")
#reflength <- sapply(N, sum)
#max_ref <- as.numeric(max(reflength))
#print(reflength)
#print(contigStatsFlipped(style="data",N=N, reflength=reflength))
#contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlim=c(0,max_count),
#                   xlab="Number of contigs",
#                   ylab="Cumulative contig length",
#                   main="Cumulative Plot of N Statistic (own)"
#)

# Use Heinz reference length for N50
#print("Use Heinz reference length for N50")
#reflength <- sapply(N, function(x) x <-as.numeric(reflength["heinz reference (2.40)"]))
#print(reflength)
#print(contigStatsFlipped(style="data",N=N, reflength=reflength))
#contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlim=c(0,max_count),
#                   xlab="Number of contigs",
#                   ylab="Cumulative contig length",
#                   main="Cumulative Plot of N Statistic (Heinz as reference length)"
#)




# Use maximal reference length for N50
#print("Use maximal reference length for N50")
#reflength <- sapply(N, sum)
#max_ref <- as.numeric(max(reflength))
#print(reflength)
#print(contigStatsFlipped(style="data",N=N, reflength=reflength))
#contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlim=c(0,max_count),
#                   xlab="Number of contigs",
#                   ylab="Cumulative contig length",
#                   main="Cumulative Plot of N Statistic (longest)"
#)
