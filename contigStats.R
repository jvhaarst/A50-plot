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

require("gdata"   , warn.conflicts=FALSE)
require("sitools" , warn.conflicts=FALSE)
require("plotmath", warn.conflicts=FALSE)
source('./rbind.na.R')

#contigStats <- function(N=N, reflength, style="ggplot2", pch=20, xlab="Percentage of Assembly Covered by Contigs of Size >=Y", ylab="Contig Size [bp]", main="Cumulative Length of Contigs", sizetitle=14, sizex=12, sizey=12, sizelegend=9, xlim, ylim) {
#        ## Compute cumulative length vectors for contig sets
#        Nl <- lapply(names(N), function(x) rev(sort(N[[x]]))); names(Nl) <- names(N)
#        Nlcum <- lapply(names(Nl), function(x) cumsum(Nl[[x]])); names(Nlcum) <- names(Nl)
#
#        ## Compute N50 values
#        N50 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x]/2 >= 0)[1]])
#        names(N50) <- names(N)
#
#        ## Return only data (no plot)
#        if(style=="data") {
#                N90 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.90 >= 0)[1]]); names(N50) <- names(N)
#                N75 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.75 >= 0)[1]]); names(N50) <- names(N)
#                N25 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.25 >= 0)[1]]); names(N50) <- names(N)
#                stats <- cbind(N25, N50, N75, N90, Longest=sapply(N, max), Mean=sapply(N, mean), Median=sapply(N, median), Shortest=sapply(N, min), N_Contigs=sapply(N, length))
#                return(c(Nlcum, Contig_Stats=list(stats)))
#        }
#        ## Plot cumulative contig length with base graphics, only necessary when ggplot is unavailable
#	if(style=="base") {
#            if(missing(xlim)) xlim <- c(0, 100)
#            if(missing(ylim)) ylim <- c(0, max(unlist(N)))
#            split.screen(c(1,1))
#            for(i in seq(along=Nl)) {
#                    if(i==1) {
#                            plot(Nlcum[[i]]/reflength[[i]] * 100, Nl[[i]], col=i, pch=pch, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
#                    }
#                    screen(1, new=FALSE)
#                    plot(Nlcum[[i]]/reflength[[i]] * 100, Nl[[i]], col=i, pch=pch, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n")
#            }
#            legend("bottomright", legend=paste(names(N50), ": N50 = ", N50, sep=""), cex=1.2, bty="n", pch=15, pt.cex=1.8, col=seq(along=Nl))
#            close.screen(all=TRUE)
#        }
#        ## Plot cumulative contig length with ggplot2
#        ## Note: ggplot2 plotting options can be looked up with theme_get()
#	if(style=="ggplot2") {
#                require("ggplot2")
#                plotdf <- data.frame(Samples=rep(names(Nlcum), sapply(Nlcum, length)), length=unlist(Nl), perc=unlist(lapply(names(Nlcum), function(x) Nlcum[[x]]/reflength[[x]]*100)))
#                counts <- table(plotdf[,1]); counts <- counts[names(N50)]
#                N50rep <- paste(plotdf[,1], ": N50=", unlist(lapply(as.character(unique(plotdf[,1])), function(x) rep(N50[x], counts[names(N50[x])]))), sep="")
#                plotdf[,1] <- N50rep
#                ggplot(plotdf, aes(perc, length, color=Samples)) +
#                geom_point() +
#                geom_line() +
#                scale_x_continuous(xlab) +
#                scale_y_continuous(ylab) +
#                labs(title = main) +
#                theme(plot.title = element_text(size = sizetitle)) +
#                theme(axis.title.x = element_text(size = sizex)) +
#                theme(axis.title.y = element_text(size = sizey, angle = 90)) +
#                theme(legend.text = element_text(size = sizelegend)) +
#                theme(legend.title = element_text(size = 0))
#        }
#}


fmtnum<-function(num){
    return(formatC(num, format="d", big.mark=","))
}

contigStatsFlipped <- function(N=N, reflength, style="ggplot2", pch=20, xlab="Percentage of Assembly Covered by Contigs of Size >=Y", ylab="Contig Size [bp]", main="Cumulative Length of Contigs", sizetitle=14, sizex=12, sizey=12, sizelegend=9, trimSize=25000, xunity=1e3, doLookup=FALSE, yunity=1e7, xlim, ylim, outBaseName="Rplots") {
        ## Compute cumulative length vectors for contig sets, trimming at TRIMSIZE
        cat("Trimming histograms\n")

		NlTrim<- lapply(names(N ), function(x) {
				nN<-rev(sort(N[[x]]))

				trimSizel<-trimSize
				if ( length(nN) < trimSize ) {
						cat(paste("reducing trim size to"  , length(nN), "\n"))
						trimSizel<-length(nN)-1
				} else {
						cat(paste("keeping trim size. size", length(nN), "\n"))
				}
                    
				nN <- trimSum(nN, trimSizel, right=TRUE, na.rm=FALSE)

				return(nN);
		}); names(NlTrim)    <- names(N)

        numEls<-length(names(N))

        cat(paste("number of elements:", numEls, "\n"))
        
		Nl        <- lapply(names(N), function(x) rev(sort(N[[x]]))  ); names(Nl       ) <- names(N)
        Nlcum     <- lapply(names(N), function(x) cumsum(Nl[[x]])    ); names(Nlcum    ) <- names(N)
        NlTrimCum <- lapply(names(N), function(x) cumsum(NlTrim[[x]])); names(NlTrimCum) <- names(N)

        ## Compute N50 values for use on graph
        N50 <- sapply(seq(along=N), function(x) {Nl[[x]][which(Nlcum[[x]] - reflength[x]/2 >= 0)[1]]}); names(N50) <- names(N)
        I50 <- sapply(seq(along=N), function(x) {        which(Nlcum[[x]] - reflength[x]/2 >= 0)[1] }); names(I50) <- names(N)
        Ns  <- sapply(seq(along=N), function(x) length(N[[x]]))                                       ; names(Ns ) <- names(N)
        
		cat("N50\n")
		cat(paste(N50, "\n"))
		cat("I50\n")
		cat(paste(I50, "\n"))
		#print(paste("Ns", Ns))

        #add tick values every 5th part of the X length
		maxx<-trimSize
		maxx<-round(maxx/xunity)*xunity
        if(missing(xlim)) xlim <- c(0, max(unlist(N    )))
		tickvaluesX<-seq(0, maxx, by=maxx/5)
		cat("Tick values X\n"    )
		cat(paste(     tickvaluesX , "\n"))
		cat(paste(f2si(tickvaluesX), "\n"))

        lookup_max<-unlist(lapply(seq(along=Nlcum), function(x) Nlcum[[x]][which.max(abs(Nlcum[[x]]))]))
        #cat(paste("qmax ", qmax, "\n"))
        
        lookup_mean<-floor(1/mean(1/lookup_max))
        cat(paste("Nlookup:  ", lookup_mean  , "\n"))
        lookup_mean75<-floor(lookup_mean * .75)
        cat(paste("Nlookup75:", lookup_mean75, "\n"))
        
        Nlookup_mean75     <- sapply(seq(along=N), function(x) {Nl[[x]][which(Nlcum[[x]] - lookup_mean75 >= 0)[1]]}); names(Nlookup_mean75) <- names(N)
        Ilookup_mean75     <- sapply(seq(along=N), function(x) {        which(Nlcum[[x]] - lookup_mean75 >= 0)[1] }); names(Ilookup_mean75) <- names(N)
        Ilookup_mean75Trim <- sapply(seq(along=N), function(x) {if ( Ilookup_mean75[x] > trimSize ) { return(trimSize); } else { return(Ilookup_mean75[x]); }} ); ; names(Ilookup_mean75Trim) <- names(N)
        cat(paste(names(N), "Nlookup_mean75", Nlookup_mean75, "Ilookup_mean75", Ilookup_mean75, "\n"))
        print(Ilookup_mean75Trim)
        #qmaxmax<-max(qmax)
        #qmaxmin<-min(qmax)
        #cat(paste("qmaxmax ", qmaxmax, "\n"))
        #cat(paste("qmaxmin ", qmaxmin, "\n"))
        #
        #qmaxmax50<-ceiling( qmaxmax/2 )
        #qmaxmin50<-ceiling( qmaxmin/2 )
        #cat(paste("qmaxmax50 ", qmaxmax50, "\n"))
        #cat(paste("qmaxmin50 ", qmaxmin50, "\n"))
        #
        #
        #Nqmaxmax50 <- sapply(seq(along=N), function(x) {Nl[[x]][which(Nlcum[[x]] - qmaxmax50 >= 0)[1]]}); names(Nqmaxmax50) <- names(N)
        #Iqmaxmax50 <- sapply(seq(along=N), function(x) {        which(Nlcum[[x]] - qmaxmax50 >= 0)[1] }); names(Iqmaxmax50) <- names(N)
        #
        #Nqmaxmin50 <- sapply(seq(along=N), function(x) {Nl[[x]][which(Nlcum[[x]] - qmaxmin50 >= 0)[1]]}); names(Nqmaxmin50) <- names(N)
        #Iqmaxmin50 <- sapply(seq(along=N), function(x) {        which(Nlcum[[x]] - qmaxmin50 >= 0)[1] }); names(Iqmaxmin50) <- names(N)
        #
        #cat(paste(names(N), "Nqmaxmax50 ", Nqmaxmax50, " Iqmaxmax50 ", Iqmaxmax50, "\n"))
        #cat(paste(names(N), "Nqmaxmin50 ", Nqmaxmin50, " Iqmaxmin50 ", Iqmaxmin50, "\n"))
        
        ## Return only data (no plot)
        if(style=="data") {
            outTsvCumm = paste(outBaseName, ".cumm.tsv", sep="")
            outTsvStat = paste(outBaseName, ".stat.tsv", sep="")

            
            N90 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.90 >= 0)[1]]); names(N90) <- names(N)
            I90 <- sapply(seq(along=N), function(x)         which(Nlcum[[x]] - reflength[x] * 0.90 >= 0)[1] ); names(I90) <- names(N)

            N75 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.75 >= 0)[1]]); names(N75) <- names(N)
            I75 <- sapply(seq(along=N), function(x)         which(Nlcum[[x]] - reflength[x] * 0.75 >= 0)[1] ); names(I75) <- names(N)

            N25 <- sapply(seq(along=N), function(x) Nl[[x]][which(Nlcum[[x]] - reflength[x] * 0.25 >= 0)[1]]); names(N25) <- names(N)
            I25 <- sapply(seq(along=N), function(x)         which(Nlcum[[x]] - reflength[x] * 0.25 >= 0)[1] ); names(I25) <- names(N)
            
            stats <- cbind(
                    #N25, I25, N50, I50, N75, I75, N90, I90, qmaxmax, qmaxmin, qmaxmax50, qmaxmin50, Nqmaxmax50, Iqmaxmax50, Nqmaxmin50, Iqmaxmin50,
                    N25, I25, N50, I50, N75, I75, N90, I90, lookup_mean, lookup_mean75, Nlookup_mean75, Ilookup_mean75,
                    Longest=sapply(N, max), Shortest=sapply(N, min), 
                    Mean=round(sapply(N, mean)), Median=round(sapply(N, median)),
                    N_Contigs=sapply(N, length), Total_length=sapply(N, sum)
                    )

            alldata<-do.call(rbind.na, NlTrimCum)
            #print(names(alldata)); quit()

            cat("Contig Number\t" , file=outTsvCumm, append=FALSE)
            cat(names(Ns)         , file=outTsvCumm, append=TRUE , sep="\t" )
            cat("\n"              , file=outTsvCumm, append=TRUE )
            write.table(t(alldata), file=outTsvCumm, sep="\t", na="\"\"", col.names=FALSE, append=TRUE )

            cat("\"Name\"\t",   file=outTsvStat,           append=FALSE)
            write.table(stats , file=outTsvStat, sep="\t", append=TRUE )
            
            return(Contig_Stats=list(stats))
        }
        ## Plot cumulative contig length with base graphics, only necessary when ggplot is unavailable
        if(style=="base") {
            outPng     = paste(outBaseName, ".png"     , sep="")

            png(filename=outPng, width=2048, height=2048, bg="transparent", type="cairo", pointsize=48)

		    #add tick values every 5th part of the X length
            maxy<-max(unlist(NlTrimCum))
			maxy<-round(maxy/yunity)*yunity

            if(missing(ylim)) ylim <- c(0, max(unlist(NlTrimCum)))
			
			tickvaluesY<-seq(0, maxy, by=maxy/5)
			
            cat("Tick values Y\n"    )
            cat(paste(     tickvaluesY , "\n"))
            cat(paste(f2si(tickvaluesY), "\n"))
            par(mar=c(4,4,5,1)+.1)
            
            palette(rainbow(numEls))
            split.screen(c(1,1))
            for(i in seq(along=NlTrim)) {
                # TODO: IF I == J50, CHANGE SYMBOL?
                if(i==1) {
                    plot(y=NlTrimCum[[i]], x=seq_along(NlTrimCum[[i]]),col=i, pch=pch, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", xlab=xlab, ylab=ylab, main=main, type='o', bty="n")
                    box(lwd=6)
                }
                screen(1, new=FALSE)

                plot(y=NlTrimCum[[i]], x=seq_along(NlTrimCum[[i]]),  col=i, pch=pch, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n", ann=FALSE, type='o')
            }
		    axis(1, at=tickvaluesX, labels=f2si(tickvaluesX))
		    axis(2, at=tickvaluesY, labels=f2si(tickvaluesY))


            legendTxt = paste(names(N50), ":",
                            " N50=" , fmtnum(N50),
                            " I50=" , fmtnum(I50),
                            " Size=", fmtnum(Ns ), sep="")
                            
            if ( doLookup ) {
                axis(3, at=tickvaluesX, labels=f2si(tickvaluesX))
    
                maxylim=max(ylim)
                minylim=min(ylim)
                maxxlim=max(xlim)
                yblock =maxylim/80
                xblock =maxxlim/350
                screen(1, new=FALSE)
                for (nPos in seq(1:length(Ilookup_mean75))) {
                    #print(paste("pos ", nPos, " I50 ", Iqmaxmax50[nPos], " ylim ", max(ylim), " I50 ", Iqmaxmax50[nPos], " qmax ", qmaxmax50))
                    #segments(Iqmaxmax50[nPos], max(ylim), Iqmaxmax50[nPos], qmaxmax50, col=nPos, lty=5, lwd=3)
                    #segments(Iqmaxmin50[nPos], qmaxmin50, Iqmaxmin50[nPos],         0, col=nPos, lty=5, lwd=3)
                    #rect(Iqmaxmax50[nPos]-xblock, maxylim+(yblock*2), Iqmaxmax50[nPos]+xblock, maxylim+(yblock*3), col=nPos, border=NA)
                    #rect(Iqmaxmin50[nPos]-xblock, minylim-(yblock*3), Iqmaxmin50[nPos]+xblock, minylim-(yblock*2), col=nPos, border=NA)
                    rect(Ilookup_mean75Trim[nPos]-xblock, maxylim+(yblock*1), Ilookup_mean75Trim[nPos]+xblock, maxylim+(yblock*3), col=nPos, border=NA)
                }


                #qmaxtxt = bquote(" Mmax50 = " ^ .(qmaxmax50))
                #qmaxtxt = paste(" Mmax50=", fmtnum(qmaxmax50))
                #text(trimSize *.85, qmaxmax50 * 1.05, qmaxtxt)
                #abline(h=qmaxmax50, lty=5, lwd=3)
                #
                #qmintxt = paste(" Mmin50=", fmtnum(qmaxmin50))
                #text(trimSize *.85, qmaxmin50 * 0.95, qmintxt)
                #abline(h=qmaxmin50, lty=5, lwd=3)

                qmaxtxt = paste("Lookup75=", fmtnum(lookup_mean75), sep="")
                text(trimSize *.85, lookup_mean75 * 1.05, qmaxtxt)
                abline(h=lookup_mean75, lty=5, lwd=8)

                
                #for (seqname in names(N)) {
                #    print(paste("adding square to name:", seqname))
                #    q50maxM = Nqmaxmax50[seqname]
                #    q50maxJ = Iqmaxmax50[seqname]
                #    print(paste("  M:", q50maxM, " J:",q50maxJ))
                #   
                #    rect(q50maxJ, maxy * .20, q50maxJ*1.01, maxy * .19)
                #}
                
                legendTxt = paste(names(N50), ":",
                            " N50=", fmtnum(N50),
                            " I50=", fmtnum(I50),
                            #" Jmin=", fmtnum(Iqmaxmin50), " Jmax=", fmtnum(Iqmaxmax50),
                            " Nl=", fmtnum(Nlookup_mean75),
                            " Il=", fmtnum(Ilookup_mean75),
                            " Size=", fmtnum(Ns), sep="")
            }

            legend("bottomright", legend=legendTxt, cex=0.75, bty="n", pch=15, pt.cex=1.5, col=seq(along=Nl), xjust=1 )

            close.screen(all=TRUE)
        }
        ## Plot cumulative contig length with ggplot2
        ## Note: ggplot2 plotting options can be looked up with theme_get()
	if(style=="baseperc") {
		    maxy<-100
			
			if(missing(ylim)) ylim <- c(0, 100)
            
			split.screen(c(1,1))
            for(i in seq(along=NlTrim)) {
                    if(i==1) {
                    	plot(y=NlTrimCum[[i]]/reflength[[i]] * 100, x=seq_along(NlTrimCum[[i]]),col=i, pch=pch, xlim=xlim, xaxt="n", yaxt="n", ylim=ylim, xlab=xlab, ylab=ylab, main=main, type='o')
                    }
                    screen(1, new=FALSE)
                    plot(y=NlTrimCum[[i]]/reflength[[i]] * 100, x=seq_along(NlTrimCum[[i]]),  col=i, pch=pch, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n", ann=FALSE, type='o')
					axis(1, at=tickvalues, labels=f2si(tickvalues))
					axis(2, at=tickvalues, labels=f2si(tickvalues))
            }
            legend("bottomright", legend=paste(names(N50), ": N50 = ", N50, " Size = ", Ns, sep=""), cex=0.6, bty="n", pch=15, pt.cex=0.8, col=seq(along=Nl), xjust=1 )
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


