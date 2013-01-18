# Simple R script to calculate N50 of assembled contigs
# also plots sum(contig_len) y against descending sorted list of contig sizes.
#
# Inputfile:
#	 inputfile is derived from FASTA headers of ABYSS >contigID length kmercov by parsing
#    contiglen in col 2 by default (else reconfig below)
#    has headers (i.e. [contigN] [len] [abs_kmer_coverage])
#
# info: alex.bossers@wur.nl

#config
lencol <- c(2) # column to use for len of contigs
file <- c('Cb14160-002_abyss4-contigs.lengthAlx')

#main
contigdata <- read.table(file, header=T, sep="\t")
sorted.data <- contigdata[order(contigdata[lencol],decreasing=T),]
total <- c(0)
maxreached <- c(0)
sumtotal <- sum(sorted.data$len)
numcontigs <- nrow(sorted.data)
#dim(contigplot) <- c(0,0)
contigplot <- data.frame(c_len=c(0),c_sum=c(0))

for (i in 1:numcontigs) {
	total <- total+sorted.data[i,lencol]
	
	#store values for plot
	contigplot[i,1] <- sorted.data[i,lencol]
	contigplot[i,2] <- total
	
	#get contig len that pushes the sum over half total sum
	if (!maxreached & total >= (sumtotal/2)) {
		cat("number of contigs    :",numcontigs,"\n")
		cat("sum total contig len :",sumtotal,"\n")
		cat("half sum contigs     :",sumtotal/2,"\n")
		cat("at index",i,"contig len",sorted.data[i,2],"summs up to",total,"greater than half_sumcontigs\n")
		cat("==> N50 determined   :",sorted.data[i,2],"bp\n")
		maxreached <- c(1)
	}
}
#finally plot the data x=contigs y=sumcontigs
plot(contigplot$c_sum,type="l",col="blue")