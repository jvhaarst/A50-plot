# R script to:
# Read in multiple genome assemblies
# Calculate metrics on scaffolds
# Generate table of metrics of all assemblies
# Generate graph to show assembly growth

# Load the needed R library (from Bioconductor)
#install.packages("gdata")
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#biocLite("IRanges")

require("Biostrings")
require("IRanges")

library("Biostrings")
library("IRanges")


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

readSeq<-function(filename){
  sortedSizes <- sort(width(read.DNAStringSet(filename)), decreasing=TRUE)
  return(sortedSizes)
}

#heinz_allpaths                                       <- readSeq("/home/assembly/dev_150/assemblies/allpaths_lg_sample_heinz_raw/sl/data/run/ASSEMBLIES/test/final.assembly.fasta")
#allpaths_454                                         <- readSeq("/home/assembly/dev_150/assemblies/allpaths_lg_sample_heinz_raw_with454/sl/data/run/ASSEMBLIES/test/final.assembly.fasta")
heinz_reference                                      <- readSeq("/home/assembly/dev_150/assemblies/S_lycopersicum_scaffolds.2.40.fa")


#habrochaites_aplg_ill_170_2k                         <- readSeq("/home/assembly/dev_150/assemblies/allpaths_lg_habrochaites_raw/sh/data/run/ASSEMBLIES/test/final.assembly.fasta")
#habrochaites_aplg_ill_170_500_2k                     <- readSeq("/home/assembly/dev_150/assemblies/allpaths_lg_habrochaites_raw/sh/data/run/ASSEMBLIES/test2/final.assembly.fasta")
#habrochaites_aplg_ill_170_2k_opera_454_8k_20k        <- readSeq("/home/aflit001/temptive/opera/habrochaites/output1/scaffoldSeqFmt.fasta")
habrochaites_aplg_ill_170_500_2k_opera_454_8k_20k    <- readSeq("/home/aflit001/temptive/opera/habrochaites/output/scaffoldSeqFmt.fasta")


#pennellii_aplg_ill_170_2k                            <- readSeq("/home/assembly/dev_150/assemblies/allpaths_lg_pennellii_raw/sp/data/run/ASSEMBLIES/test/final.assembly.fasta")
#pennellii_aplg_ill_170_500_2k                        <- readSeq("/home/assembly/dev_150/assemblies/allpaths_lg_pennellii_raw/sp/data/run/ASSEMBLIES/test2/final.assembly.fasta")
#pennellii_aplg_ill_170_2k_opera_454_3k_8k_20k        <- readSeq("/home/aflit001/temptive/opera/pennellii/output1/scaffoldSeqFmt.fasta")
pennellii_aplg_ill_170_500_2k_opera_454_3k_8k_20k    <- readSeq("/home/aflit001/temptive/opera/pennellii/output/scaffoldSeqFmt.fasta")


#clc                                                  <- readSeq("/home/assembly/dev_150/assemblies/clc-default/clc_contigs.fa")
#arcanum_clc_default                                  <- readSeq("/home/assembly/dev_150/assemblies/clc_arcanum/CLC-780MB-tryout1.fa")
arcanum_clc_nondefault                               <- readSeq("/home/assembly/dev_150/assemblies/clc_arcanum/CLC-830MB-tryout2.fa")


# Create named list of contig lengths. Could do that in one go, but this is a bit more flexible, and lets me change visual names to something more friendly
N <- list(
  "heinz reference (2.40)"                                = heinz_reference,
  #"heinz_allpaths"                                        = heinz_allpaths,
  
  #"habrochaites_allpaths"                                 = habrochaites_allpaths,
  #"habrochaites_opera_scaf"                               = habrochaites_opera_scaf,
  #"pennellii_allpaths"                                    = pennellii_allpaths,
  #"pennellii_opera_scaf"                                  = pennellii_opera_scaf,
  
  #"habrochaites (aplg ill 170+2k)"                        = habrochaites_aplg_ill_170_2k,
  #"habrochaites (aplg ill 170+500+2k)"                    = habrochaites_aplg_ill_170_500_2k,
  #"habrochaites (aplg ill 170+2k)+(opera 454 8k+20k)"     = habrochaites_aplg_ill_170_2k_opera_454_8k_20k,
  "habrochaites (aplg ill 170+500+2k)+(opera 454 8k+20k)" = habrochaites_aplg_ill_170_500_2k_opera_454_8k_20k,

  #"pennellii (aplg ill 170+2k)"                           = pennellii_aplg_ill_170_2k,
  #"pennellii (aplg ill 170+500+2k)"                       = pennellii_aplg_ill_170_500_2k,
  #"pennellii (aplg ill 170+2k)+(opera 454 3k+8k+20k)"     = pennellii_aplg_ill_170_2k_opera_454_3k_8k_20k,
  "pennellii (aplg ill 170+500+2k)+(opera 454 3k+8k+20k)" = pennellii_aplg_ill_170_500_2k_opera_454_3k_8k_20k,

  #arcanum_clc_default                                     = arcanum_clc_default,
  "arcanum (clc ill 170+2k 454 8k+20k)"                   = arcanum_clc_nondefault
  )

# Get the maximum contig count from the list
max_count <- max(unlist(lapply(N,length)))
#max_count <- 30000
source('./contigStats.R')

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

max_count <- 20000

# Use maximal reference length for N50
print("Use maximal reference length for N50")
reflength <- sapply(N, sum)
max_ref <- as.numeric(max(reflength))
print(reflength)
stats<-contigStatsFlipped(style="data",N=N, reflength=reflength)
print(stats)
write.table(stats, file="Rplots.txt", sep="\t", append=FALSE)
contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlim=c(0,max_count),
                    trimSize=max_count,
                    xlab="Number of contigs",
                    ylab="Cumulative contig length",
                    main="Cumulative Plot of N Statistic (longest)"
)
