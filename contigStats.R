###########################################################################################
## Compute N50 values for Assembly Results and Plot their Cumulative Length Distribution ##
###########################################################################################
## Author: Thomas Girke
## Last update: January 28, 2010
## Utilities: Compute N50 and plot cumulative length distribution.
## More details can be found here: http://manuals.bioinformatics.ucr.edu/home/ht-seq#TOC-Analyzing-Assembly-Results
## Originally at http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/contigStats.R
## Changed by : jan.vanhaarst@wur.nl

## Overview:
## The N50 is a weighted median statistic such that 50% of the entire assembly
## is contained in contigs equal to or larger than this value. It can be calculated
## as follows: store the lengths of all contigs in a vector L of positive integers.
## Create a vector Lr where each integer is repeated as many times as defined by
## the sum its values in L. After this the N50 value is computed as the median of Lr.
## For example: if L = c(2, 2, 3, 3, 5, 7), then Lr consists of four 2's, six 3's,
## five 5's, and seven 7's; the N50 of L is the median of Lr, which is 5.
## In R the N50 according to this definition can be computed like this:
##      L <- c(2, 2, 3, 3, 5, 7)
##      Lr <- unlist(tapply(L, L, function(x) rep(x[1], sum(x))))
##      median(Lr)
##
## With very large contig sizes the above implementation for computing the N50
## value can become inefficient. A much more efficient approach is to sort L
## decreasingly, compute for it the cumulative sum and then identify the length
## value for which the cumulative sum covers at least 50% of the combined length
## of all contigs or the length of a known target/reference sequence. Note,
## comparisons among N50 values from different assemblies are only meaningful
## when using for their calculation the same combined length value. Thus, a known
## target length value can be often a good solution for comparing assembly results.

## The following contigStats function calculates the N50 values according to the more
## efficient second method. In addition, it creates a distribution plot of the
## cumulative contig lengths, which is a very effective way for comparing assembly
## results.

#saulo
require("gdata")
library("gdata")
source('./rbind.na.R')
#oluas

contigStats <- function(N=N, reflength, style="ggplot2", pch=20, xlab="Percentage of Assembly Covered by Contigs of Size >=Y", ylab="Contig Size [bp]", main="Cumulative Length of Contigs", sizetitle=14, sizex=12, sizey=12, sizelegend=9, xlim, ylim) {
        ## Compute cumulative length vectors for contig sets
        Nl <- lapply(names(N), function(x) rev(sort(N[[x]]))); names(Nl) <- names(N)
        Nlcum <- lapply(names(Nl), function(x) cumsum(Nl[[x]])); names(Nlcum) <- names(Nl)

        ## Compute N50 values
        N50 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x]/2 >= 0)[1]])
        names(N50) <- names(N)

        ## Return only data (no plot)
        if(style=="data") {
                N90 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.90 >= 0)[1]]); names(N50) <- names(N)
                N75 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.75 >= 0)[1]]); names(N50) <- names(N)
                N25 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.25 >= 0)[1]]); names(N50) <- names(N)
                stats <- cbind(N25, N50, N75, N90, Longest=sapply(N, max), Mean=sapply(N, mean), Median=sapply(N, median), Shortest=sapply(N, min), N_Contigs=sapply(N, length))
                return(c(Nlcum, Contig_Stats=list(stats)))
        }
        ## Plot cumulative contig length with base graphics, only necessary when ggplot is unavailable
	if(style=="base") {
            if(missing(xlim)) xlim <- c(0, 100)
            if(missing(ylim)) ylim <- c(0, max(unlist(N)))
            split.screen(c(1,1))
            for(i in seq(along=Nl)) {
                    if(i==1) {
                            plot(Nlcum[[i]]/reflength[[i]] * 100, Nl[[i]], col=i, pch=pch, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
                    }
                    screen(1, new=FALSE)
                    plot(Nlcum[[i]]/reflength[[i]] * 100, Nl[[i]], col=i, pch=pch, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n")
            }
            legend("bottomright", legend=paste(names(N50), ": N50 = ", N50, sep=""), cex=1.2, bty="n", pch=15, pt.cex=1.8, col=seq(along=Nl))
            close.screen(all=TRUE)
        }
        ## Plot cumulative contig length with ggplot2
        ## Note: ggplot2 plotting options can be looked up with theme_get()
	if(style=="ggplot2") {
                require("ggplot2")
                plotdf <- data.frame(Samples=rep(names(Nlcum), sapply(Nlcum, length)), length=unlist(Nl), perc=unlist(lapply(names(Nlcum), function(x) Nlcum[[x]]/reflength[[x]]*100)))
                counts <- table(plotdf[,1]); counts <- counts[names(N50)]
                N50rep <- paste(plotdf[,1], ": N50=", unlist(lapply(as.character(unique(plotdf[,1])), function(x) rep(N50[x], counts[names(N50[x])]))), sep="")
                plotdf[,1] <- N50rep
                ggplot(plotdf, aes(perc, length, color=Samples)) +
                geom_point() +
                geom_line() +
                scale_x_continuous(xlab) +
                scale_y_continuous(ylab) +
                labs(title = main) +
                theme(plot.title = element_text(size = sizetitle)) +
                theme(axis.title.x = element_text(size = sizex)) +
                theme(axis.title.y = element_text(size = sizey, angle = 90)) +
                theme(legend.text = element_text(size = sizelegend)) +
                theme(legend.title = element_text(size = 0))
        }
}
contigStatsFlipped <- function(N=N, reflength, style="ggplot2", pch=20, xlab="Percentage of Assembly Covered by Contigs of Size >=Y", ylab="Contig Size [bp]", main="Cumulative Length of Contigs", sizetitle=14, sizex=12, sizey=12, sizelegend=9, trimSize=25000, xlim, ylim) {
        ## Compute cumulative length vectors for contig sets
		
        Nl    <- lapply(names(N ), function(x) {
				nN<-rev(sort(N[[x]]))

				trimSizel<-trimSize
				if ( length(nN) < trimSize ) {
						print(paste("reducing trim size to", length(nN)))
						trimSizel<-length(nN)-1
				} else {
						print(paste("keeping trim size. size ", length(nN)))
				}
				nN <- trimSum(nN, trimSizel, right=TRUE, na.rm=FALSE)

				return(nN);
		}); names(Nl)    <- names(N)
        Nlcum <- lapply(names(Nl), function(x) cumsum(Nl[[x]]));   names(Nlcum) <- names(Nl)


        ## Compute N50 values for use on graph
        N50 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x]/2 >= 0)[1]])
        Ns  <- sapply(seq(along=N), function(x) length(N[[x]]))
        names(N50) <- names(N)
        names(Ns ) <- names(N)
		#print(paste("Ns", Ns))

        ## Return only data (no plot)
        if(style=="data") {
                N90 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.90 >= 0)[1]]); names(N50) <- names(N)
                N75 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.75 >= 0)[1]]); names(N50) <- names(N)
                N25 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.25 >= 0)[1]]); names(N50) <- names(N)
                stats <- cbind(
                		N25, N50, N75, N90, 
                		Longest=sapply(N, max), Shortest=sapply(N, min), 
                		Mean=round(sapply(N, mean)), Median=round(sapply(N, median)),
                		N_Contigs=sapply(N, length), Total_length=sapply(N, sum)
                		)

				alldata<-do.call(rbind.na, Nlcum)
				#print(names(alldata)); quit()

				cat("Contig Number\t" , file="Rplots_cumm.csv", append=FALSE)
				cat(names(Ns)         , file="Rplots_cumm.csv", append=TRUE , sep="\t" )
				cat("\n"              , file="Rplots_cumm.csv", append=TRUE )
				write.table(t(alldata), file="Rplots_cumm.csv", sep="\t", na="\"\"", col.names=FALSE, append=TRUE )

                return(Contig_Stats=list(stats))
        }
        ## Plot cumulative contig length with base graphics, only necessary when ggplot is unavailable
	if(style=="base") {
            if(missing(xlim)) xlim <- c(0, max(unlist(N)))
            #if(missing(ylim)) ylim <- c(0, 100)
            if(missing(ylim)) ylim <- c(0, max(unlist(Nlcum)))
            split.screen(c(1,1))
            for(i in seq(along=Nl)) {
                    if(i==1) {
                    	#plot(y=Nlcum[[i]]/reflength[[i]] * 100, x=seq_along(Nlcum[[i]]),col=i, pch=pch, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
                    	plot(y=Nlcum[[i]], x=seq_along(Nlcum[[i]]),col=i, pch=pch, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
                    }
                    screen(1, new=FALSE)
                    #plot(y=Nlcum[[i]]/reflength[[i]] * 100, x=seq_along(Nlcum[[i]]),  col=i, pch=pch, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n")
                    plot(y=Nlcum[[i]], x=seq_along(Nlcum[[i]]),  col=i, pch=pch, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n")
            }
            legend("bottomright", legend=paste(names(N50), ": N50 = ", N50, " Size = ", Ns, sep=""), cex=0.6, bty="n", pch=15, pt.cex=0.8, col=seq(along=Nl),
                   xjust=1
                   )
            close.screen(all=TRUE)
        }
        ## Plot cumulative contig length with ggplot2
        ## Note: ggplot2 plotting options can be looked up with theme_get()
	if(style=="baseperc") {
            if(missing(xlim)) xlim <- c(0, max(unlist(N)))
            if(missing(ylim)) ylim <- c(0, 100)
            split.screen(c(1,1))
            for(i in seq(along=Nl)) {
                    if(i==1) {
                    	plot(y=Nlcum[[i]]/reflength[[i]] * 100, x=seq_along(Nlcum[[i]]),col=i, pch=pch, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
                    }
                    screen(1, new=FALSE)
                    plot(y=Nlcum[[i]]/reflength[[i]] * 100, x=seq_along(Nlcum[[i]]),  col=i, pch=pch, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n")
            }
            legend("bottomright", legend=paste(names(N50), ": N50 = ", N50, " Size = ", Ns, sep=""), cex=0.6, bty="n", pch=15, pt.cex=0.8, col=seq(along=Nl),
                   xjust=1
                   )
            close.screen(all=TRUE)
        }
        ## Plot cumulative contig length with ggplot2
        ## Note: ggplot2 plotting options can be looked up with theme_get()
	if(style=="ggplot2") {
                require("ggplot2")
                plotdf <- data.frame(Samples=rep(names(Nlcum), sapply(Nlcum, length)), length=unlist(Nl), perc=unlist(lapply(names(Nlcum), function(x) Nlcum[[x]]/reflength[[x]]*100)))
                counts <- table(plotdf[,1]); counts <- counts[names(N50)]
                N50rep <- paste(plotdf[,1], ": N50=", unlist(lapply(as.character(unique(plotdf[,1])), function(x) rep(N50[x], counts[names(N50[x])]))), sep="")
                plotdf[,1] <- N50rep
                ggplot(plotdf, aes(perc, length, color=Samples)) +
                geom_point() +
                geom_line() +
                scale_x_continuous(xlab) +
                scale_y_continuous(ylab) +
                labs(title = main) +
                theme(plot.title = element_text(size = sizetitle)) +
                theme(axis.title.x = element_text(size = sizex)) +
                theme(axis.title.y = element_text(size = sizey, angle = 90)) +
                theme(legend.text = element_text(size = sizelegend)) +
                theme(legend.title = element_text(size = 0))
        }
}
## Usage of contigStats function
## Test sample of simulated contig length values provided in list
# 
#N1 <- sample(1:500, 100)
# 
#N2 <- sample(1:480, 80)
# 
#N <- list(N1=N1, N2=N2)

## Run contigStats function
# 
#reflength <- sapply(N, sum)
#max_ref <- max(reflength)
#reflength <- sapply(N, function(x) x <- max_ref)
# 
#N50 <- contigStats(N=N, reflength=reflength, pch=20, xlab="Percentage of Assembly Covered by Contigs of Size >=Y", ylab="Contig Size", main="Cumulative Plot of N Statistic")
#N50 <- contigStats(style="data",N=N, reflength=reflength, pch=20, xlab="Percentage of Assembly Covered by Contigs of Size >=Y", ylab="Contig Size", main="Cumulative Plot of N Statistic")


