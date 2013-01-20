# R script to: 
# Read in multiple genome assemblies
# Split scaffolds into contigs, based on N distribution
# Calculate metrics on scaffolds and contigs
# Generate table of metrics of all assemblies
# Generate graph to show assembly growth

# Load the needed R library (from Bioconductor)
require("Biostrings")
#require("ggplot2")
#require("Cairo")
# Read in assembly in FASTA format 
# Read filename from commandline
#args=(commandArgs(TRUE))
#input=args[1]
#print(input)
input <- "/home/jvh/data/degradome_tomato/SL2.40ch.fasta"
testdata<-read.DNAStringSet(input, "fasta")
# Calculate N50 (TODO put this in function)
assemblyStats <- function(input=input, name=name){
        # Contig length data
        N <- width(input)
        # Total length
        reflength <- sum(N)
        # Created a sorted,cumulative list.
        N_sorted <- sort(N,decreasing=TRUE)
        N_sorted_cumulative <- cumsum(N_sorted)
        # Find the length at which the cumulative length exceeds half of the total length
        # This works like this :
        # sapply(seq_along(N), function(x))                     : loop from 1 to length of N, so number of contigs
        # which(N_sorted_cumulative[[x]] - reflength/2 >= 0)[1]	: return NA or 1, depending on whether the result is true or false
        # Which then gives(if true)				                      : N_sorted[704][1]
        # After that we need to get the first item that is not NA
        N90 <- sapply(seq_along(N), function(x) N_sorted[x][ which(N_sorted_cumulative[x] - reflength/(100/90) >= 0)[1] ] )
        N_90_index = which(N90 != 'NA')[1]
        N_90 <- N90[N_90_index]
        rm(N90)
        
        N50 <- sapply(seq_along(N), function(x) N_sorted[x][ which(N_sorted_cumulative[x] - reflength/(100/50) >= 0)[1] ] )
        N_50_index = which(N50 != 'NA')[1]
        N_50 <- N50[N_50_index]
        rm(N50)
        
        N25 <- sapply(seq_along(N), function(x) N_sorted[x][ which(N_sorted_cumulative[x] - reflength/(100/25) >= 0)[1] ] )
        N_25_index = which(N25 != 'NA')[1]
        N_25 <- N_sorted[N_25_index]
        rm(N25)
        
        total_contig_length =  sum(width(testdata))
        Max_contig_length = max(width(testdata))
        Number_of_contigs_larger_than_100 = length(which(sapply(seq_along(N), function(x) N[x]>100 )))
        Number_of_bases_in_contigs_over_100 = sum(na.omit(sapply(seq_along(N), function(x) N_sorted[x][which(N_sorted[x]>=100)[1] ] )))
        Number_of_contigs_larger_than_1000 = length(which(sapply(seq_along(N), function(x) N[x]>=1000 )))
        Number_of_bases_in_contigs_over_1000 = sum(na.omit(sapply(seq_along(N), function(x) N_sorted[x][which(N_sorted[x]>=1000)[1] ] )))
        Percentage_of_Ns_in_consensus = 100*(alphabetFrequency(testdata, collapse=TRUE)['N']/total_contig_length)
        
        Mean_contig_length = mean(width(testdata))
        Standard_deviation_of_contig_length = sd(width(testdata))
        Median_contig_length = median(width(testdata))
        Total_known_base_number=alphabetFrequency(testdata, collapse=TRUE)['A']+
        			alphabetFrequency(testdata, collapse=TRUE)['C']+
        			alphabetFrequency(testdata, collapse=TRUE)['G']+
        			alphabetFrequency(testdata, collapse=TRUE)['T']
        
        GC_Content_of_contigs=100*(alphabetFrequency(testdata, collapse=TRUE)['G']+alphabetFrequency(testdata, collapse=TRUE)['C'])/Total_known_base_number
        #print(alphabetFrequency(testdata, baseOnly=TRUE,collapse=TRUE))
        
        stats <- cbind(name,reflength,
                       N_25,N_25_index,N_50,N_50_index,N_90,N_90_index,
                       total_contig_length,Max_contig_length,
                       Number_of_contigs_larger_than_100,Number_of_bases_in_contigs_over_100,
                       Number_of_contigs_larger_than_1000,Number_of_bases_in_contigs_over_1000,
                       Percentage_of_Ns_in_consensus,
                       Mean_contig_length,Standard_deviation_of_contig_length,Median_contig_length,
                       GC_Content_of_contigs)
        return( list( array(N_sorted_cumulative), array(N_sorted), stats ) )
}

stats <- assemblyStats(input=testdata,name=input)
print(summary(unlist(stats[1])))
plot(unlist(stats[1]))
stop("Stop")
#sequences<-as.character(testdata, use.names=FALSE)
#gregexpr("N+",sequences)

# Plot the data nicely with ggplot2

# Plot the cumulative length, with N50/N25/N90 lines 
#ggplot(plotdf, aes(perc, length, color=Samples)) + 
#geom_point() + 
#scale_x_continuous(xlab) + 
#scale_y_continuous(ylab) + 
#opts(title = main) +
#opts(plot.title = theme_text(size = sizetitle)) +
#opts(axis.title.x = theme_text(size = sizex)) + 
#opts(axis.title.y = theme_text(size = sizey, angle = 90)) +
#opts(legend.text = theme_text(size = sizelegend)) +
#opts(legend.title = theme_text(size = 0))
#CairoX11(display=Sys.getenv("DISPLAY"), width = 10, height = 10)
#CairoPDF("plot.pdf", width = 6, height = 6, onefile = TRUE, bg="white")
# Multiple plots are multiple pages
plot(N_sorted_cumulative,main="Cumulative length" )
# Put data in same plot with :
#par(new=T)
#plot(cumsum(na.omit(sapply(seq_along(N), function(x) N_sorted[x][which(N_sorted[x]>=1000)[1] ] ))))


plotdata <- data.frame(N_sorted_cumulative)
plotdata$counter<-row(plotdata)

#hist(width(testdata),breaks="scott")

# Read in contigStats
#source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/contigStats.R")
#contigStats(N=N, reflength=reflength, style="ggplot2")
#stats <- contigStats(N=N, reflength=reflength, style="data")
#stats[["Contig_Stats"]]

#temp<-gregexpr("A+T*|A*T+",a[[7]])[[1]]
#
#temp1<-data.frame(chr=names(a)[7],pos=as.numeric(temp),l=attr(temp,which="match.length"))
#
###Split
#a<-"AAAAAAAAAAANNNNNNNNNNNTTTTTTTTT"
#strsplit(a,"NNNNN")
