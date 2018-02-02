# plot functions, similar to one in molecular nets

#' Map x into a set of ranges
mapRange = function(x, ranges){
	if(is.na(x))
		x = 0
	for(i in 1:(length(ranges)-1)){
		if(x >= ranges[i] & x < ranges[i+1])
			return(i)
	}
	if( x <= ranges[1] )
		return(1)
	if (x >= ranges[length(ranges)])
		return(length(ranges))
	return(NA)
}

myPalette = function(){
	palette.pos = apply(colorRamp(c("white", brewer.pal(9, "OrRd")))((1:25)/25), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
	palette.neg = apply(colorRamp(c("white", brewer.pal(9, "PuBu")))((1:25)/25), 1, function(x) rgb(x[1]/255, x[2]/255, x[3]/255))
    palette = c(rev(palette.neg[1:length(palette.neg)]), palette.pos)
    
    palette
}

#' Make a plot of error distribution
#' 
#' @param err the matrix of errors
#' @param main the main tile
#' @param cex.main scaling for main title
#' @param cex.axis scaling for the labels
#' @param cex.text scaling for the numbers in the table
#' @param val.col the colour matrix, if NULL will be generated based on matrix values
#' @param round.num if to round numbers to two decimals
#' @param ranges the ranges of colour-coding. Based on this and using mapRange() the colour-coding is calculated. It is a vector of numbers.
#' @param rev.pal if to reverse the default pallette
#' @param gray.diag if to gray out the diagonals
#' @param bidir the 2-column matrix of bidirectional edges (same format as GBNet), used to gray them out if neccessary
#' @param nocol.last.row if to not colour the last row (which are wildtype values)
#' @param suffix suffix to add to printer numbers, e.g. "%"
#' @param right.labels the labels to show at the right side of the plot (axis=4)
#' @param show.text if to show the textual representation of the values
#' @param las.top the las style of the annotations at the top of the plot
#' @param top.side the side of the top labels
#' @param row.round.num the per-row round.num
#' @param row.nocol rows that should not be coloured
#' @param rotate.top if to rotate the top axis labels by 45 degress
#' @param palette a custom palette to use
#' @param plot if to actually plot something
#' @param row.labels rowname labels to use if rownames() is meaningless
plot.err = function(err, main="", cex.main=1.2, cex.axis=0.9, cex.text=0.9, val.col=NULL, round.num=2, ranges=NULL, rev.pal=FALSE, gray.diag=TRUE,
bidir = NULL, nocol.last.row=FALSE, suffix="", right.labels=NULL, show.text=TRUE, las.top=0, top.side=3, row.round.num=NULL, row.nocol=NULL, rotate.top=FALSE,
palette=NULL, plot=TRUE, row.labels=NULL){
	#palette.pos = c("#FFFFFF", brewer.pal(9, "OrRd"))[1:7]
    #palette.neg = c("#FFFFFF", brewer.pal(9, "PuBu"))[1:7]	
    if(is.null(palette))
		palette = myPalette()
    
    sizex = ncol(err)
    sizey = nrow(err)
    
    if(rev.pal)
    	palette = rev(palette)
    
    if( is.null(val.col) ){
    	if(is.null(ranges))
			ranges = seq(-0.35, 0.35, length.out=length(palette))
			
		val.col = matrix("", ncol=sizex, nrow=sizey)
		for(i in 1:sizey){
			for(j in 1:sizex){
				val.col[i,j] = palette[ mapRange( err[i,j], ranges ) ]
			}
		}
		
		if(gray.diag)
			diag(val.col) = "lightgray"
			
		if(!is.null(bidir)){
			if(nrow(bidir)>0){
				for(i in 1:nrow(bidir)){
					val.col[bidir[i,1], bidir[i,2]] = "lightgray"
					val.col[bidir[i,2], bidir[i,1]] = "lightgray"
				}
			}
		}
    }
    
    if(nocol.last.row){
    	val.col[nrow(val.col),] = "white"
    }
    
    if(!is.null(row.nocol)){
    	val.col[row.nocol,] = "white"
    
    }
    
    if(plot){
    	if(is.null(row.labels))
    		row.labels = rownames(err)
    
        plot(NULL, type = "n", yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(0-0.01,sizex+0.01), ylim=c(0-0.01,sizey+0.01), xaxs="i", yaxs="i",
                    main=main, cex.main=cex.main)
        rasterImage(val.col, 0, 0, sizex, sizey, interpolate=F)
        axis(2, at=1:sizey-0.5, labels=rev(row.labels), tick=F, mgp=c(3,0.2,0), cex.axis=cex.axis, las=1)
        if(rotate.top){
        	text((1:sizex)-0.5, par("usr")[4] + 0.15, srt = 45, adj = 0, labels = colnames(err), xpd = TRUE, cex=cex.axis)
        } else {
        	axis(top.side, at=1:sizex-0.5, labels=colnames(err), tick=F, mgp=c(3,0.2,0), cex.axis=cex.axis, las=las.top)
       	}
        if(!is.null(right.labels))
        	axis(4, at=1:sizey-0.5, labels=rev(right.labels), tick=F, mgp=c(3,0.2,0), cex.axis=cex.axis, las=1)
        
        if(show.text){
		    for(i in 1:sizex){
			    for(j in 1:sizey){
				    if(is.na(err[j,i])){
					    label = ""
				    } else if(round.num>=0){
					    if(j == nrow(val.col) & nocol.last.row){
						    label = paste(err[j,i])
					    } else {
						    if(!is.null(row.round.num)){
							    label = paste(sprintf(paste("%.",row.round.num[j],"f",sep=""), round(err[j,i],row.round.num[j])), suffix, sep="")
						    } else {
							    label = paste(sprintf(paste("%.",round.num,"f",sep=""), round(err[j,i],round.num)), suffix, sep="")
						    }
					    }
				    } else {
					    label = paste(err[j,i])
				    }
		            text(i-0.5, sizey-j+0.5, label, cex=cex.text)
			    }
		    }
        }
    }
    
	invisible(val.col)

}

