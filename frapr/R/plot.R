# Plotting functions


#' Helper functions for plotFL to map to coordinates from image()
#'
#' @param x x coordinate (columns)
#' @param nc number of columns in the simulation matrix
.mapX = function(x, nc){
	width = 1 + 1/(nc-1)
	start = 0.5/(nc-1)
	
	(x/nc - start) * width
}

#' Helper functions for plotFL to map to coordinates from image()
#'
#' @param y y coordinate (rows)
#' @param nr number of rows in the simulation matrix
.mapY = function(y, nr){
	height = 1 + 1/(nr-1)
	start = 0.5/(nr-1)
	
	(y/nr - start) * height
}


#' Plot a rectangular fluorescent simulation area
#'
#' @param x a matrix of fluorescence values
#' @param max.value the brighest green value
#' @param digits digits to show for time
#' @param bleachMask the frap mask to visualise
#' @param bleach.col colour of the frap bleach mask
#' @param bleach.lty line type of the frap bleach mask area
plotFL = function(x, t=0, max.value=max(x), digits=2, bleachMask=NULL, bleach.col="white", bleach.lty=3){
	col.map = sapply(1:255, function(x) rgb(0, x/255, 0))
	col.breaks = seq(0, max.value, length.out=256)

	if(!is.null(t)){
	    par(mar=c(4,4,4,4))
	} else {
	    par(mar=c(0,0,0,0))
	}
	# set black background
	par(bg="black")
	# plot image and time
	showImage(x, col=col.map, breaks=col.breaks, useRaster=TRUE, xaxt="n", yaxt="n")
	if(!is.null(t))
	    mtext(sprintf(paste("t = %.", digits, "f s", sep=""), t), side=1, line=1, at=1, adj=c(1,1), col="white")
	
	# plot the frap (bleaching) mask
	if(!is.null(bleachMask)){
		#bleachMask = t(bleachMask) # to match what we output with image()
		bleachMask = bleachMask[nrow(bleachMask):1,]
		# normalisation to drawing area
		nc = ncol(bleachMask)
		nr = nrow(bleachMask)
		for(i in 1:nrow(bleachMask)){
			for(j in 1:ncol(bleachMask)){
				# transition from mask -> no mask 
				if(bleachMask[i,j] && (j == ncol(bleachMask) || !bleachMask[i,j+1])){
					lines(.mapX(c(j,j),nc), .mapY(c(i-1,i),nr), col=bleach.col, lwd=2, lty=bleach.lty)
				# no mask -> mask
				} else if (!bleachMask[i,j] && j != ncol(bleachMask) && bleachMask[i,j+1]){
					lines(.mapX(c(j,j),nc), .mapY(c(i-1,i),nr), col=bleach.col, lwd=2, lty=bleach.lty)
				# first	
				} else if(j == 1 && bleachMask[i,j]){
					lines(.mapX(c(j-1,j-1),nc), .mapY(c(i-1,i),nr), col=bleach.col, lwd=2, lty=bleach.lty)
				}
				
				# same for horizontal lines... 
				if(bleachMask[i,j] && (i == nrow(bleachMask) || !bleachMask[i+1,j])){
					lines(.mapX(c(j-1,j),nc), .mapY(c(i,i),nr), col=bleach.col, lwd=2, lty=bleach.lty)
				} else if(!bleachMask[i,j] && i != nrow(bleachMask) && bleachMask[i+1,j]){
					lines(.mapX(c(j-1,j),nc), .mapY(c(i,i),nr), col=bleach.col, lwd=2, lty=bleach.lty)
				} else if(i==1 && bleachMask[i,j]){
					lines(.mapX(c(j-1,j),nc), .mapY(c(i-1,i-1),nr), col=bleach.col, lwd=2, lty=bleach.lty)
				}
			}
		}
	}
}

#' Save simulation results as an animation
#' 
#' @param sim the return value of gillespieSim()
#' @param outfile output filename
#' @param start.time start time to use (inclusive), defaults to 0
#' @param end.time end time to use (inclusive), default to the whole simulation
#' @param time.speed how fast should simulation time pass with respect to real time, defaults to 1.
#'        A value of e.g. 0.5 would slow down the simulation 2x time. This is achieved by showing
#'        the frames from the simulation at a slower rate.  
#' @param max.fl maximal fluorescene value (if not specified will use max(sim$FL))
saveAnimation = function(sim, outfile, start.time=0, end.time=max(sim$t), time.speed=1, max.fl=max(sim$FL)){
	p = sim$params
	
	# set the paramters for speed of animation
	# HACK: note we are injecting -dither none not to lose in quality ('animation' doesn't support this argument)
	ani.options(interval=1/p$fps * time.speed, loop="0 -dither none")
	# turn off automatic play
	ani.options(autobrowse=FALSE, autoplay=FALSE)

	if(file.exists("animation.gif"))
		file.remove("animation.gif")

	saveGIF({
		time.sel = which(sim$t >= start.time & sim$t <= end.time)
		for(i in time.sel){
			plotFL(sim$FL[,,i], sim$t[i], max.fl, bleachMask=p$bleachMask)
		}
	})
	
	file.copy("animation.gif", outfile, overwrite=TRUE)
	file.remove("animation.gif")
}

#' Moving average
#' @param x the vector of values
#' @param r range of values
movingAvg = function(x, r){
	y = x
	for(i in 1:length(x)){
		r1 = i - r
		r2 = i + r
		
		if(r1 < 1)
			r1 = 1
		if(r2 > length(x))
			r2 = length(x)
			
		y[i] = mean(x[r1:r2])
	}
	
	y
}

#' Plot the recovery curves
#'
#' @param img set of images (of real data)
#' @param p simulation parameters
#' @param res.time residence time
#' @param Free proportion of free molecules
#' @param D diffusion constant
#' @param experiment experiment name
#' @param replicate the replicate number
#' @param mfrow layout of the graphs
#' @param smooth smoothing radius to use
#' @param cex.text scaling factor for the final inferred parameters panel
#' @param sim.mean simulation results if already available
#' @param xlim custom x limit for the axis
plotParamCurves = function(img, p, res.time=NULL, Free=NULL, D=NULL, experiment="", replicate="", mfrow=c(2,3), smooth=0, 
	cex.text=1, sim.mean=NULL, xlim=NULL){
	if(!is.null(res.time) && !is.null(Free) && !is.null(D)){
		p = initRates(p, residence.time=res.time, Free=Free, D=D)
	}
	if(is.null(sim.mean))
		sim.mean = simulateMeans2D(p)
	
	ylim = range(na.omit(c(sim.mean$curve, img$stats$rec)))
	ylim[1] = 0
	
	#main = paste(experiment, "D =",  p$D/1e-12, "um2/s, res.time =", res.time, "s, Free =", Free)
	par(mfrow=mfrow)
	for(i in 1:ncol(img$stats$rec)){		
		y = img$stats$rec[,i]
		if(smooth != 0){
			y = movingAvg(y, smooth)
		}
		
		col = "blue"
		if(i == 1)
			col = "black"
		
		if(is.null(xlim))
			xlim = range(img$t)
		
		plot(img$t, y, main=colnames(sim.mean$curve)[i], xlim=xlim, ylim=ylim, 
			xlab="Time (s)", ylab="% recovery", col="darkgray")
		lines(img$t, y, col="darkgray")
		lines(sim.mean$t + p$bleach.time, sim.mean$curve[,i], col=col, lwd=2)
		
		abline(v=p$bleach.time, lty=2)
		
		if(i != 1)
			lines(lines(sim.mean$t + p$bleach.time, sim.mean$curve[,1], col="black", lwd=1, lty=2))
		
		abline(h=1, lty=2, col="darkgray")		
	}
	
	# plot the parameters
	plot(NULL, xlim=c(0, 1), ylim=c(0,1), xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	text(0.05, 0.95, experiment, cex=cex.text, adj=c(0,1))
	text(0.05, 0.85, paste("Replicate", replicate), cex=cex.text, adj=c(0,1))
	text(0.05, 0.52, paste("D =", D/1e-12, "um2/s"), cex=cex.text, adj=c(0,1))
	text(0.05, 0.52-0.12, paste("Free molecules =", round(Free*100), "%"), cex=cex.text, adj=c(0,1))
	text(0.05, 0.52-0.24, paste("Residence time =", paste(res.time, collapse="s, "), "s"), cex=cex.text, adj=c(0,1))	
}

#' Plot the error rate in a consistent format
#'
#' NOTE: this function has an external dependency to plot.err! 
#'
#' @param err a 3D array of error values
#' @param mfrow how to layout the graphs, default is automatic
#' @param ... other parameters to plot.err
plotErr = function(err, mfrow=NULL, ...){
	# automatic mfrow selection
	if(is.null(mfrow)){
		d = dim(err)[3]
		if(d == 1){
			mfrow=c(1,1)
		} else if(d>1 & d<=2){
			mfrow=c(1,2)
		} else if(d>2 & d<=4){
			mfrow=c(2,2)
		} else if(d>4 & d<=6){
			mfrow=c(2,3)
		} else if(d>6 & d<=9){
			mfrow=c(3,3)
		} else if(d>9 & d<=12){
			mfrow=c(3,4)
		} else {
			mfrow=c(4,4)
		}
	}
	
	# range for colour-coding
	r = range(-log(err))
	r = seq(r[1], r[2], length.out=50)
		
	# plot stuffs
	par(mfrow=mfrow, mar=c(5,4,6,2))
	for(i in 1:dim(err)[3]){
		plot.err(-log(err[,,i]), gray.diag=FALSE, ranges=r, main=paste("D =",  as.numeric(dimnames(err)[[3]][i])/1e-12, "um2/s"), ...)
		mtext(side=2, line=1.2, "Residence time (s)", cex=0.6)
		mtext(side=3, line=1.2, "Proportion free", cex=0.6)
	}
}

