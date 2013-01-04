# R script to: 
# Read in multiple genome assemblies
# Calculate metrics on scaffolds
# Generate table of metrics of all assemblies
# Generate graph to show assembly growth

# Load the needed R library (from Bioconductor)
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#biocLite("IRanges")
require("Biostrings")
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
#   contigs <-readDNAStringSet(assembly, "fasta",nrec=10000)
#   assign(paste(as.name(assembly)),width(contigs))
#   assembly_contigs <- c(assembly_contigs,paste(as.name(assembly)))
#   #assembly_contigs <- c(assembly_contigs,paste(as.name(assembly))=width(contigs))
# }
# N <- list(assemblies)
# reflength <- sapply(N, sum)
# max_ref <- max(reflength)
# reflength <- sapply(N, function(x) x <- max_ref)

heinz_allpaths           <- width(readDNAStringSet("/home/assembly/dev_150/assemblies/allpaths_lg_sample_heinz_raw/sl/data/run/ASSEMBLIES/test/final.assembly.fasta"))
#allpaths_454            <- width(readDNAStringSet("/home/assembly/dev_150/assemblies/allpaths_lg_sample_heinz_raw_with454/sl/data/run/ASSEMBLIES/test/final.assembly.fasta"))
#clc                     <- width(readDNAStringSet("/home/assembly/dev_150/assemblies/clc-default/clc_contigs.fa"))
arcanum_clc_default      <- width(readDNAStringSet("/home/assembly/dev_150/assemblies/clc_arcanum/CLC-780MB-tryout1.fa")) 
arcanum_clc_nondefault   <- width(readDNAStringSet("/home/assembly/dev_150/assemblies/clc_arcanum/CLC-830MB-tryout2.fa"))
heinz_reference          <- width(readDNAStringSet("/home/assembly/dev_150/assemblies/S_lycopersicum_scaffolds.2.40.fa"))
#fermi                   <- width(readDNAStringSet("/home/assembly/progs/fermi/heinz/fmdef.p4.fa"))
habrochaites_allpaths    <- width(readDNAStringSet("/home/assembly/dev_150/assemblies/allpaths_lg_habrochaites_raw/sh/data/run/ASSEMBLIES/test/final.assembly.fasta"))
habrochaites_opera_scaf  <- width(readDNAStringSet("/home/aflit001/temptive/opera/habrochaites/output/scaffoldSeq.fasta"))
pennellii_allpaths       <- width(readDNAStringSet("/home/assembly/dev_150/assemblies/allpaths_lg_pennellii_raw/sp/data/run/ASSEMBLIES/test/final.assembly.fasta"))
pennellii_opera_scaf     <- width(readDNAStringSet("/home/aflit001/temptive/opera/pennellii/output/scaffoldSeq.fasta"))

N <- list(
  heinz_reference=heinz_reference,
  heinz_allpaths=heinz_allpaths,
  habrochaites_allpaths=habrochaites_allpaths,
  habrochaites_opera_scaf=habrochaites_opera_scaf,
  pennellii_allpaths=pennellii_allpaths,
  pennellii_opera_scaf=pennellii_opera_scaf,
  arcanum_clc_default=arcanum_clc_default,
  arcanum_clc_nondefault=arcanum_clc_nondefault
  )

source('~/code/assemblystats/contigStats.R')
# Use own reference length for N50
print("Use own reference length for N50")
reflength <- sapply(N, sum)
max_ref <- as.numeric(max(reflength))
print(reflength)
print(contigStatsFlipped(style="data",N=N, reflength=reflength))
#contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlab="Percentage of Assembly Covered by Contigs of Size >=Y", ylab="Contig Size", main="Cumulative Plot of N Statistic")

contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlim=c(0,30000),
                   xlab="Number of contigs", 
                   ylab="Cumulative contig length", 
                   main="Cumulative Plot of N Statistic (own)"
)
# Use Heinz reference length for N50
print("Use Heinz reference length for N50")
reflength <- sapply(N, function(x) x <-as.numeric(reflength["heinz_reference"]))
print(reflength)
print(contigStatsFlipped(style="data",N=N, reflength=reflength))
#contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlab="Percentage of Assembly Covered by Contigs of Size >=Y", ylab="Contig Size", main="Cumulative Plot of N Statistic")

contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlim=c(0,30000),
                   xlab="Number of contigs", 
                   ylab="Cumulative contig length", 
                   main="Cumulative Plot of N Statistic (heinz)"
)

# Use maximal reference length for N50
print("Use maximal reference length for N50")
reflength <- sapply(N, function(x) x <- max_ref)
print(reflength)
print(contigStatsFlipped(style="data",N=N, reflength=reflength))

contigStatsFlipped(style="base",N=N, reflength=reflength, pch=20, xlim=c(0,30000),
                   xlab="Number of contigs", 
                   ylab="Cumulative contig length", 
                   main="Cumulative Plot of N Statistic (longest)"
)
