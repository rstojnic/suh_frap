# image properties file

#' Load the post-bleaching frame
#' If post.avg > 1, this will actually be an average of multiple frames
#' @param img the img list with the image properties
loadPostFrame = function(img){
    post.avg = img$post.avg
    if(post.avg == 1){
        post.img = loadTiff(paste(img$base.path, "/", img$frames[img$post.frame], sep=""))
        list(t=img$t[img$post.frame], t.inx=img$post.frame, post.img=post.img)
    } else {
        # load images and average
        l = list()
        for(i in 1:post.avg){
            l[[i]] = loadTiff(paste(img$base.path, "/", img$frames[img$post.frame+i-1], sep=""))
        }
        
        # average over images
        l.sum = l[[1]]
        for(i in 2:length(l))
            l.sum = l.sum + l[[i]]
        l.mean = l.sum / length(l)
        
        # figure out time
        t.inx = img$post.frame:(img$post.frame+post.avg-1)
        list(t=mean(img$t[t.inx]), t.inx=t.inx, post.img=l.mean)
    }
}

#' Make the image properties object
#'
#' @param experiment experiment name
#' @param replicate replicate number
#' @param base.path path where the TIFF files are
#' @param frame.pat file pattern for frames, needs to contain %F as a placeholder for the frame number
#' @param bleach.frames frames that correspond to bleaching
#' @param time.inc time increment per frame (in seconds)
#' @param real.len.x real length (in meters) along the x axis
#' @param real.len.y real length (in meters) along the y axis
#' @param remove.frames frames to remove
#' @param roi.file filename of the ImageJ ROIs
#' @param roi.total names of the total area ROI
#' @param roi.bg.pat name pattern for backgrounds
#' @param roi.exclude name of the ROI exclusion pattern
#' @param post.avg how many post-bleaching frames to average, if 1 then don't average just take the one frame
#' @param roi.total.file if the total ROI is in a different file... 
#' @param bleach.duration the time in second that the bleaching took place
#' @param bleach.norm.sd the standard deviation of the normal distribution used to model the laser bleaching profile (default: 4 microns)
#' @param bleach.delay the delay between the end of the bleaching, and beginning of imaging (default: 0.5s)
#' @param bleach.norm.sd.range the range of values to try for bleach.norm.sd, defaults to a single value
makeImageProps = function(experiment, replicate, base.path, frame.pat, bleach.frames, time.inc, real.len.x, real.len.y, remove.frames=NULL, roi.file="rois.zip", 
	roi.total="total", roi.bg.pat="background", roi.exclude="", post.avg=3, roi.total.file=NULL, bleach.duration=1, bleach.norm.sd=1.7*1e-6,
	bleach.delay=0.5, bleach.norm.sd.range=bleach.norm.sd){
	img = list()
	
	img$experiment = experiment
	img$replicate = replicate
	img$base.path = base.path
	img$frame.pat = frame.pat
	img$time.inc = time.inc
	img$real.len.x = real.len.x
	img$real.len.y = real.len.y
	img$post.avg = post.avg
	img$bleach.duration = bleach.duration
	img$bleach.norm.sd = bleach.norm.sd
	img$bleach.norm.sd.range = bleach.norm.sd.range
	img$bleach.delay = bleach.delay
	
	# double-check the arguments
	if(!str_detect(frame.pat, "%F"))
		stop("Variable 'frame.pat' needs to contain '%F' to designate where the frame number is located in the file name.")
		
	if(!file.exists(base.path))
		stop("Path in 'base.path' does not exist.")
		
	if(length(bleach.frames) == 0)
		stop("Vector of bleach frames 'bleach.frames' need to have at least one value")
		
	# connect the files to frames... 
	frames = dir(base.path, pattern=glob2rx(gsub("%F", "*", frame.pat, fixed=TRUE)))
	frames.inx = as.numeric(unlist(strapply(frames, gsub("%F", "([0-9]+)", frame.pat), identity)))
	o = order(frames.inx)
	
	if(length(frames) == 0)
		stop("No frames found. Please check your input parameters.")
	
	# order frames by their number
	frames = frames[o]
	frames.inx = frames.inx[o]
	
	if(!is.null(remove.frames)){
		to.remove = frames.inx %in% remove.frames
		frames = frames[!to.remove]
		frames.inx = frames.inx[!to.remove]
	}
	
	if(!(all(frames.inx[-1] - frames.inx[-length(frames.inx)] == 1))){
		stop("Some frames are missing (frame numbers are non-consecutive)!")
	}

	# frames and time
	img$frames = frames
	img$t = frames.inx * time.inc
		
	# convert the bleach frame indicies from original... 
	img$bleach.frames = which(frames.inx %in% bleach.frames)	
	# record which frames are pre and post-bleaching
	img$pre.frames = 1:(min(img$bleach.frames)-1)
	img$post.frame = max(img$bleach.frames)+1
	
	pf = loadPostFrame(img)
	img$post.img = pf$post.img
	if(post.avg > 1){
	    # Need to modify time and frames
	    # remove the times that were averaged over
        img$t = c(img$t[1:(pf$t.inx[1]-1)], pf$t, img$t[-(1:pf$t.inx[length(pf$t.inx)])])	    
        # remove the frames that were averaged over
        img$frames = c(img$frames[1:(pf$t.inx[1]-1)], "should never tried to load this image", img$frames[-(1:pf$t.inx[length(pf$t.inx)])])	    
	}

	# get the image properties from the first frame
	f = readTIFF(paste(img$base.path, "/", img$frames[1], sep=""))
	img$size.x = dim(f)[2]
	img$size.y = dim(f)[1]
	
	# now load the ROIs
	rois.raw = read.ijzip(paste(base.path, "/", roi.file, sep=""))
	
	# Sometimes the Total area is in a different file.. 
	if(!is.null(roi.total.file)){
	    rois.total.raw = read.ijzip(paste(base.path, "/", roi.total.file, sep=""))
	} else {
	    rois.total.raw = rois.raw
	}
	
	# extract ROIs
	if(!(roi.total %in% names(rois.total.raw))){
		stop(paste("The specified 'roi.total' value of '", roi.total, "' is not present in the loaded ROIs", sep=""))
	}
	
	src.size = c(img$size.x, img$size.y)
	img$rois = lapply(rois.raw, roiToMask, src.size)
	rois.total.converted = lapply(rois.total.raw, roiToMask, src.size)
	
	img$roi.total = rois.total.converted[[roi.total]]
	# subtract the exclusion ROI if present
	if(roi.exclude %in% names(img$rois)){
		img$roi.total = img$roi.total & !img$rois[[roi.exclude]]
	}
	if(roi.exclude %in% names(rois.total.converted)){
	    img$roi.total = img$roi.total & !rois.total.converted[[roi.exclude]]
	}
	
	
	img$roi.background = img$rois[names(rois.raw)[grep(roi.bg.pat, names(rois.raw))]]
	
	# adjust the ROIs so there is no interpolation error when the images is caled down
	
	img

}

#' Draw a line into a boolean mask matrix
#'
#' @param m a boolean matrix
#' @param x0 starting x coordinate
#' @param y0 starting y coordinate
#' @param x1 ending x coordinate
#' @param y1 ending y coordinate
.drawLine = function(m, x0, y0, x1, y1){
	# make sure these are all integers
	x0 = as.integer(x0)
	x1 = as.integer(x1)
	y0 = as.integer(y0)
	y1 = as.integer(y1)

	dx = abs(x1 - x0)
	dy = abs(y1 - y0)

	# special case: horizontal line
	if(dy == 0){
		if(y0 >= 1 && y0 <= nrow(m))
			m[y0, max(1,x0):min(x1,ncol(m))] = TRUE
		return(m)
	} else if(dx == 0){
		# special case: vertical line
		if(x0 >= 1 && x0 <= ncol(m))
			m[max(1,y0):min(y1,nrow(m)), x0] = TRUE
		return(m)
	}

	sx = ifelse(x0 < x1, 1L, -1L)
	sy = ifelse(y0 < y1, 1L, -1L)
	
	err = dx - dy
	
	while(TRUE){
		# only draw inside the bounds
		if(y0 >= 1 && y0 <= nrow(m) && x0 >= 1 && x0 <= ncol(m))
			m[y0, x0] = TRUE
		
		if(x0 == x1 && y0 == y1)
			break
		
		e2 = 2*err
		if(e2 > -dy){
			err = err - dy
			x0 = x0 + sx
		}
		if(e2 < dx){
			err = err + dx
			y0 = y0 + sy
		}
	}
		
	return(m)
}

#' Convert a ROI into a binary mask (old version!)
#'
#' @param roi a ROI object from RImageJROI
#' @param src.size a vector of size 2 specifying the x and y size of the source ROI
#' @param dest.size a vector of size 2 specifying how big the destination size should be
roiToMask = function(roi, src.size, dest.size=src.size) {
	must_be("ijroi", roi)
	must_be("numeric", src.size, dest.size)
	must_have("length", 2, src.size, dest.size)

	# output mask
	mask = matrix(FALSE, nrow=dest.size[2], ncol=dest.size[1])

	# draw a polygon
	if(roi$strType %in% c("polygon", "traced") || (roi$strType == "freehand" && !("strSubtype" %in% names(roi)))){
		cds = roi$coords
		
		# rescale for the output coordinates (which are originally 0-based)
		cds[,"x"] = cds[,"x"] / (src.size[1]-1) * (dest.size[1]-1) + 1
		cds[,"y"] = cds[,"y"] / (src.size[2]-1) * (dest.size[2]-1) + 1
		
		# draw the individual lines
		for(i in 1:nrow(cds)){
			c1 = cds[i,]
			if(i == nrow(cds)){
				c2 = cds[1,]
			} else {
				c2 = cds[i+1,]
			}
			
			mask = .drawLine(mask, c1[1], c1[2], c2[1], c2[2])			
		}		
		
		# fill the contour
		outmask = mask
		for(i in 1:nrow(mask)){
			inside = FALSE
			j.vals = which(mask[i,])
			if(length(j.vals) > 0){
				for(j in seq(min(j.vals), max(j.vals))){
					if(mask[i,j] && (j==1 || !mask[i,j-1]))
						inside = !inside
		
					if(inside)
						outmask[i,j] = TRUE
				}
			}
		}

		mask = outmask
		
	} else if(roi$strType == "oval"){
		cds = roi$coords

		# rescale for the output coordinates (which are originally 0-based)
		cds[,"x"] = cds[,"x"] / src.size[1] * dest.size[1] + 1
		cds[,"y"] = cds[,"y"] / src.size[2] * dest.size[2] + 1

		height = abs(cds["right","y"]  - cds["left","y"])/2
		width = abs(cds["right","x"]  - cds["left","x"])/2

		origin.x = (cds["right","x"]  + cds["left","x"])/2
		origin.y = (cds["right","y"]  + cds["left","y"])/2
		
		
		# this value might overflow!!!
		hhww = height*height*width*width
		h2 = height*height
		w2 = width*width
		
		if(is.infinite(hhww)){
			stop("Numerical overflow when drawing the oval ROI, please use a smaller image!")
		}
		
		# the algorithm from http://stackoverflow.com/questions/10322341/simple-algorithm-for-drawing-filled-ellipse-in-c-c
		for(y in seq(-height, height)) {
			for(x in seq(-width, width)) {
				if(x*x*h2+y*y*w2 <= hhww){
					# limit the pixel to the edge
					ny = origin.y + y
					nx = origin.x + x
					
					if(ny <= 0)
						ny = 1
					if(ny > nrow(mask))
						ny = nrow(mask)
					if(nx <= 0)
						nx = 1
					if(nx > ncol(mask))
						nx = ncol(mask)
						
				    mask[ny, nx] = TRUE
				}
			}
		}
	} else if(roi$strType == "freehand" && roi$strSubtype == "ELLIPSE"){
		# convert coordinates
		x1 = roi$x1 / src.size[1] * dest.size[1] + 1
		x2 = roi$x2 / src.size[1] * dest.size[1] + 1
		y1 = roi$y1 / src.size[2] * dest.size[2] + 1
		y2 = roi$y2 / src.size[2] * dest.size[2] + 1
	
		# code from plot.ijroi
		centerX = (x1 + x2)/2
		centerY = (y1 + y2)/2
		theta = seq(0, 2*pi, len=360)
		dx = x2 - x1
		dy = y2 - y1
		major = sqrt(dx^2 + dy^2)
		minor = major*roi$aspectRatio
		a = major/2
		b = minor/2
		phi = atan2(dy, dx)
		ellipX = centerX + a*cos(theta)*cos(phi) - b*sin(theta)*sin(phi)
		ellipY = centerY + a*cos(theta)*sin(phi) + b*sin(theta)*cos(phi)
		
		cds = round(cbind(ellipX, ellipY))
		
		# draw the individual lines
		for(i in 1:nrow(cds)){
			c1 = cds[i,]
			if(i == nrow(cds)){
				c2 = cds[1,]
			} else {
				c2 = cds[i+1,]
			}
			
			mask = .drawLine(mask, c1[1], c1[2], c2[1], c2[2])			
		}		
		
		# fill the contour
		outmask = mask
		for(i in 1:nrow(mask)){
			inside = FALSE
			j.vals = which(mask[i,])
			if(length(j.vals) > 0){
				for(j in seq(min(j.vals), max(j.vals))){
					if(mask[i,j] && (j==1 || !mask[i,j-1]))
						inside = !inside
		
					if(inside)
						outmask[i,j] = TRUE
				}
			}
		}

		mask = outmask
	} else if(FALSE && roi$strType == "rect"){
		# NOTE: this code is disabled!
		cds = roi$coords
		
		# convert coordinates
		cds[,"x"] = cds[,"x"] / (src.size[1]-1) * (dest.size[1]-1) + 1
		cds[,"y"] = cds[,"y"] / (src.size[2]-1) * (dest.size[2]-1) + 1

		left = unlist(cds["left",])
		right = unlist(cds["right",])
		
		for(x in left[1]:right[1]){
			for(y in left[2]:right[2]){
				mask[y,x] = TRUE
			}
		}
		
	} else {
		stop(paste("ROI type '", roi$strType, "' not supported", sep=""))
	}
	
	if(!any(mask)){
	    stop("Loading one of the ROIs failed")
	}
	
	return(mask)
}

loadTiff = function(infile, blur.sd=5){
    raw = readTIFF(infile)
    
    # these should be single-channel so merge channels
    if(length(dim(raw)) == 3){
        img.raw = raw[,,1] + raw[,,2] + raw[,,3]
    } else {
        img.raw = raw
    }
    
    if(!is.null(blur.sd)){
        #if(require("EBImage")){
        #    img.raw = gblur(img.raw, blur.sd)
        #} else {
            # fall back to the slower library here... 
            img.raw = as.matrix(blur(as.im(img.raw), blur.sd))
        #}
    }
    img.raw
}

#' Extract the recovery curve for a ROI and a set of image
#' 
#' @param pr image property list
#' @param roi.fg names of foreground ROIs, others (total and bg) should be in the property list
#' @param custom.bg if to use custom BG values
#' @param scale.by if to scale down the image
imageStats = function(img, roi.fg, custom.bg=NULL, scale.by=1){
	frames = img$frames
	
	# make the ROIs for total, background and foregrounds
	src.size = c(img$size.x, img$size.y)
	roi = list()
	roi$total = img$roi.total
	roi$bg = img$roi.background
	roi$fg = img$rois[roi.fg]
	
	# scale down ROIs
	for(j in 1:length(roi)){
	    if(is.matrix(roi[[j]])){
	        roi[[j]] = .scaleDownROI(roi[[j]], scale.by)
	    } else {
	        if(length(roi[[j]]) != 0){
	            for(k in 1:length(roi[[j]])){
	                roi[[j]][[k]] = .scaleDownROI(roi[[j]][[k]], scale.by)
	            }
	        }
	    }
	}

	# raw fluorescence from different regions
	fl = list()
	fl$total = rep(NA, length(frames))
	fl$bg = matrix(NA, nrow=length(frames), ncol=length(roi$bg))
	fl$fg = matrix(NA, nrow=length(frames), ncol=length(roi$fg))
	
	colnames(fl$bg) = names(img$roi.background)
	colnames(fl$fg) = roi.fg

	for(i in 1:length(frames)){
		#cat("Processing image", i, "/", length(frames), "\n")
		
		# skip bleached frames
		if(i %in% img$bleach.frames)
			next
		
	    if(i == img$post.frame){
	        img.raw = img$post.img
	    } else {
		    img.raw = loadTiff(paste(img$base.path, "/", frames[i], sep=""))
	    }

		img.raw = .scaleDownImage(img.raw, scale.by)
				# get the raw values for the different components
		for(j in 1:length(roi)){
			if(is.matrix(roi[[j]])){
				fl[[j]][i] = mean(img.raw[ roi[[j]] ])
			} else {
				if(length(roi[[j]]) != 0){
					for(k in 1:length(roi[[j]])){
						fl[[j]][i,k] = mean(img.raw[ roi[[j]][[k]] ])
					}
				}
			}
		}
	}
	
	# now do the double normalisation
	if(!is.null(custom.bg)){
		bg = custom.bg
		if(!is.matrix(bg)){
			fl$bg = matrix(bg, ncol=1, nrow=length(fl$total))
		} else {
			if(nrow(bg) < length(fl$total))
				stop("Matrix background is not long enough. It should have at least as many rows as there is frames.")
			fl$bg = bg[1:nrow(fl$fg),,drop=FALSE]
			bg = rowMeans(fl$bg)
		}
	} else {
		bg = rowMeans(fl$bg)
	}
	pre.sel = img$pre.frames
	
	# record values
	norm = list()
	norm$fg = fl$fg - bg
	norm$total = fl$total - bg
	norm$fg.pre = colMeans(norm$fg[pre.sel,,drop=FALSE])
	norm$total.pre = mean(norm$total[pre.sel])
	
	recovery = t(t(norm$fg) / norm$fg.pre) / (norm$total / norm$total.pre)
	
	img$stats = list(raw = fl, norm=norm, rec=recovery)	
	img
}

#' Use image() to correctly draw an image... 
#'
#' @param img a matrix to draw
#' @param ... other parameters to pass to image()
showImage = function(img, ...){
	img = t(img)
	img = img[, ncol(img):1]
	image(img, ...)
}

#' Extend the selectors to be exactly divisible by scale.by
#' 
#' @param sel a vector of selectors
#' @param scale.by bu which number is should be exactly devisible
.extendArea = function(sel, scale.by){
	r = range(sel)

	r.min = r[1]
	r.max = r[2]
	
	# if not divisible, go to a lower value
	if((r.min-1) %% scale.by != 0){
	    r.min = (r.min %/% scale.by) * scale.by + 1
	}
	
	# if not divisible, go the a higher value
	if(r.max %% scale.by != 0){
	    r.max = (r.max %/% scale.by + 1) * scale.by
	}
		
	r.min:r.max
}

#' Scale down ROIs by a certain factor, use majority rule
#' 
#' @param roi a boolean matrix of the ROI
#' @param scale.by by how many pixels to scale down
.scaleDownROI = function(roi, scale.by){
    if(scale.by == 1)
        return(roi)
    
	out = matrix(FALSE, nrow=nrow(roi)/scale.by, ncol=ncol(roi)/scale.by)
	
	for(i in 1:nrow(out)){
		for(j in 1:ncol(out)){
			sel.i = ((i-1)*scale.by) : (i*scale.by-1) + 1
			sel.j = ((j-1)*scale.by) : (j*scale.by-1) + 1
			
			if(mean(roi[sel.i, sel.j]) >= 0.5)
				out[i,j] = TRUE
		}
	}
	
	out
}

#' Scale a ROI up by a fixed scaling factor
#'
#' @param roi scaled down ROI
#' @param scale.by integer scaling factor
.scaleUpROI = function(roi, scale.by){
    assert_that(is.matrix(roi))
    assert_that(is.numeric(scale.by), scale.by == as.integer(scale.by))
    
    out.roi = matrix(FALSE, nrow=nrow(roi)*scale.by, ncol=ncol(roi)*scale.by)
    inx = 1
    for(i in 1:ncol(roi)){
        out.col = rep(roi[,i], each=scale.by)
        # copy over the output column 'scale.by' times
        for(j in 1:scale.by){
            out.roi[,inx] = out.col
            inx = inx + 1
        }
    }
    out.roi
}

#' Scale down image by a certain factor, use majority rule
#' 
#' @param img a real values matrix of the image intensity
#' @param scale.by by how many pixels to scale down
.scaleDownImage = function(img, scale.by){
    if(scale.by == 1)
        return(img)
    
	out = matrix(0, nrow=nrow(img)/scale.by, ncol=ncol(img)/scale.by)
	
	for(i in 1:nrow(out)){
		for(j in 1:ncol(out)){
			sel.i = ((i-1)*scale.by) : (i*scale.by-1) + 1
			sel.j = ((j-1)*scale.by) : (j*scale.by-1) + 1
			
			out[i,j] = mean(img[sel.i, sel.j])
		}
	}
	
	out
}

#' Delete any negative values
#'
#' @param x input matrix
.posOnly = function(x) {
	if(any(x < 0))
		x[x<0] = 0
	
	x
}

#' Init a simulation from a set of images
#'
#' @param img image parameters and stats
#' @param scale.by by how many pixels to rescale the image, default is 10, so that a 10x10 pixel area
#'                 will become a single pixel. 
initFromImage = function(img, scale.by=10, model="full"){
	must_be_in(c("diffusion", "full"), model)
	
	stats = img$stats
	
	# get the total area
	roi.total = img$roi.total
	area.sel.x = which(colSums(roi.total)>0)
	area.sel.y = which(rowSums(roi.total)>0)
	
	# numbers of pixels to add so we have exactly divisible pixels
	area.sel.x = .extendArea(area.sel.x, scale.by)
	area.sel.y = .extendArea(area.sel.y, scale.by)
	
	# load the first frame after bleaching
	post.frame = img$post.frame
	post.img = img$post.img
	post.img = post.img[area.sel.y, area.sel.x]
	
	# now create the shape matrix and bleach area by scaling down ROIs
	shapeMat = .scaleDownROI(img$roi.total[area.sel.y, area.sel.x], scale.by)
	# NOTE: by convention the first in fg is the bleaching area!!!
	bleach.roi = colnames(img$stats$rec)[1]
	other.roi = colnames(img$stats$rec)[2]
	bleachMask = .scaleDownROI(img$rois[[bleach.roi]][area.sel.y, area.sel.x], scale.by)
	otherMask = .scaleDownROI(img$rois[[other.roi]][area.sel.y, area.sel.x], scale.by)
	
	originalBleachMask = img$rois[[bleach.roi]][area.sel.y, area.sel.x]
	originalScaleBy = scale.by
	
	# real lengths in meters
	real.x = img$real.len.x / img$size.x * length(area.sel.x)
	real.y = img$real.len.y / img$size.y * length(area.sel.y)
	
	# trim the bleach mask to fit into the shape... 
	bleachMask = bleachMask & shapeMat
	otherMask = otherMask & shapeMat

	# scale image! <- these are the values after bleaching
	img.scaled = .scaleDownImage(post.img, scale.by)
	
	pre = matrix(0, nrow=nrow(img.scaled), ncol=ncol(img.scaled))
	pre.mean = rep(0, length(img$pre.frames)) # this value is for debug only
	pre.total = rep(0, length(img$pre.frames)) # total fluorescene per frame
	# reload the images to scale the abundance with pre-bleach intensity
	for(f.inx in img$pre.frames){
		frame.img = loadTiff(paste(img$base.path, "/", img$frames[f.inx], sep=""))
		frame.img = frame.img[area.sel.y, area.sel.x]
		frame.scaled = .scaleDownImage(frame.img, scale.by)
		
		# normalised fluorescence (to TOTAL) for each frame
		frame.bg = mean(stats$raw$bg[f.inx,])
		frame.pre.norm = .posOnly(frame.scaled - frame.bg) / (mean(frame.scaled[shapeMat]) - frame.bg)
		pre = pre + frame.pre.norm
		pre.mean[f.inx] = mean(frame.pre.norm[bleachMask])
		pre.total[f.inx] = sum(frame.scaled[shapeMat])
	}
	
	pre = pre / length(img$pre.frames)
	
	# total values of fluorescence before and after (unnormalised)
	total.pre = mean(pre.total)
	total.post = sum(img.scaled[shapeMat])
	
	# corrected post-bleach image
	bg.post = mean(stats$raw$bg[post.frame,])
	# NOTE: dividing by mean(pre) assumes there is no receptor blinding!
	if(model == "diffusion"){
		img.post = .posOnly((img.scaled - bg.post) / (mean(img.scaled[shapeMat]) - bg.post) / pre)
		#img.post = .posOnly(img.scaled - bg.post) / mean(pre[shapeMat])
	} else {
		img.post = .posOnly((img.scaled - bg.post) / (mean(img.scaled[shapeMat]) - bg.post) / mean(pre[shapeMat]))
	    #img.post = .posOnly((img.scaled - bg.post) / (mean(img.scaled[shapeMat]) - bg.post) / pre)
	}
	
	#img.post = img.scaled #- bg.post

	# clean up	
	img.post[!shapeMat] = 0

	# create a list of parameters
	p = list()
	p$grid.x = grid.x = ncol(img.post)		
	p$grid.y = grid.y = nrow(img.post)
	p$area.sel.x = area.sel.x
	p$area.sel.y = area.sel.y
	if(model == "diffusion"){	
		F.init = img.post
		F.init[!shapeMat] = 0
		
		# assume diffusion only... 
		p$F.init = F.init
		p$C.init = matrix(0, nrow=grid.y, ncol=grid.x)
	} else if(model == "full"){
		# Calculate the I.prop parameter
		I.prop = matrix(0, nrow=grid.y, ncol=grid.x)		
		I.prop[shapeMat] = pre[shapeMat] / mean(pre[shapeMat])
		
		# proportion value
		p$I.prop = I.prop
		
		# initial fluorescence vaue
		I.init = img.post	
		I.init[!shapeMat] = 0	
		p$I.init = I.init
		
		p$pre.image.sim = matrix(0, nrow=grid.y, ncol=grid.x)
		p$pre.image.sim[shapeMat] = pre[shapeMat]
		
		p$post.image.sim = img.post
		p$total.pre = total.pre
		p$total.post = total.post
		# total amount of removed fluorescence
		p$bleach.fl.prop = p$total.post / p$total.pre 
	}
	
	# no reaction, these need to be set separately
	p$k.off = 0
	p$k.on = 0
	# simulation time: 5s
	p$max.time = img$time.inc * length(img$frames)
	# total length : 20um
	p$L.x = real.x
	p$L.y = real.y
	# diffusion: 30um^2/s
	p$D = 30 * (10^-6)^2
	# shape of the area
	p$shapeMat = shapeMat
	# bleaching stuffs
	p$bleachMask = bleachMask
	p$fps = 1/img$time.inc
	p$bleach.duration = img$bleach.duration
	
	# other parameters
	p$L = 0
	p$bleach.time = img$t[img$post.frame]
	p$bleach.depth = stats$rec[post.frame,1]
	p$mol = 0
	
	p$originalBleachMask = originalBleachMask
	p$originalScaleBy = scale.by
	p$bleach.norm.sd = img$bleach.norm.sd
	p$bleach.norm.sd.range = img$bleach.norm.sd.range
	p$bleach.delay = img$bleach.delay
	
	# have all the ROIs here
	p$rois = lapply(img$rois[colnames(stats$rec)], function(x) {
		m = .scaleDownROI(x[area.sel.y, area.sel.x], scale.by)
		m & shapeMat
	})
	
	# use the time from img
	p$t = img$t[(max(img$bleach.frames)+1):length(img$t)]
	p$t = p$t - min(p$t)

	p	
}

#' Initialize k.on and k.off rates based on number of residence time 
#' and number of free molecules
#'
#' @param p simulation parameters
#' @param residence.time 1/k.off in seconds. If a single value is supplied
#'        that is used as the global residence time. Alternatively, 
#'        multiple values can be provided where names correspond to
#'        different ROIs in which k_off is set. Parameter with no
#'        name is set to be the total, e.g. c(20, "band"=40)
#' @param Free proportion of total free molecules
#' @param D diffusion constant
initRates = function(p, residence.time, Free, D){
	if(Free == 1){
		# diffusion only model
		p$Free = Free
		p$D = D
		p$k.off = 0
		p$k.on = 0
		p$F.init = p$shapeMat + 0
		p$C.init = p$shapeMat * 0
		p$F.init[p$bleachMask] = p$F.init[p$bleachMask] * p$bleach.depth
	} else {
		if(length(residence.time) == 1){
			p$k.off = 1/residence.time
		} else {
			# multiple values
			target.rois = names(residence.time)
			total.inx = which(target.rois == "")
			if(length(total.inx) != 1){
				stop("Residence time needs to contain exactly one value for global residence time")
			}
			k.off = (p$I.prop > 0) * (1/residence.time[total.inx])
			# now add the exceptions
			for(roi.name in setdiff(target.rois, "")){
				if(!(roi.name %in% names(p$rois)))
					stop(paste("ROI with name '", roi.name, "' not found", sep=""))
				k.off[ p$rois[[roi.name]] ] = 1/residence.time[roi.name]
			}
		
			p$k.off = k.off
		}
		p$Free = Free
		p$D = D
	
		# varying ON rates
		p$k.on = p$k.off / p$Free * (p$I.prop - p$Free)
		if(any(p$k.on < 0))
			p$k.on[p$k.on < 0] = 0

		# K_d
		p$kd = p$k.off / p$k.on

		# initial number of Free molecules
		F.init = p$I.init * p$kd / (1+ p$kd)
		if(any(!is.finite(F.init)))
			F.init[!is.finite(F.init)] = 0

		# at these locations all the protein is free and there is no DNA!
		if(any(F.init[p$shapeMat] == 0)){
			F.init[F.init == 0 & p$shapeMat] = p$I.init[F.init == 0 & p$shapeMat]
		}	
		C.init = p$I.init - F.init
		
		# Sep 2016 change: use the steady states to initialise the non-bleached area!!!!
		# this is to minimise the amount of influence of noise on the results
		ss = steadyState(sum(p$I.init[p$shapeMat]), p$shapeMat, p$k.on, p$k.off)
		p$F.init = ss$F
		p$C.init = ss$C
		
		# only copy over the bleaced pixels
		#p$F.init[p$bleachMask] = F.init[p$bleachMask]
		#p$C.init[p$bleachMask] = C.init[p$bleachMask]
		# Copy over only the bleach depth
		p$F.init[p$bleachMask] = p$F.init[p$bleachMask] * p$bleach.depth
		p$C.init[p$bleachMask] = p$C.init[p$bleachMask] * p$bleach.depth
		
	}
	
	p
}

#' Inner function that does the actual fitting of k.on/k.off, F.init, C.init
#' @param p the vector of parameters
#' @param scale.k.on the scaling parameter (if to apply)
initRatesOldFit = function(p, scale.k.on=1){
    # varying ON rates
    p$k.on = p$k.off / p$Free * (p$I.prop - p$Free)
    if(any(p$k.on < 0))
        p$k.on[p$k.on < 0] = 0
    
    p$k.on = p$k.on * scale.k.on
    
    # K_d
    p$kd = p$k.off / p$k.on
    
    # initial number of Free molecules
    p$F.init = 	p$I.init * p$kd / (1+ p$kd)
    if(any(!is.finite(p$F.init)))
        p$F.init[!is.finite(p$F.init)] = 0
    
    # at these locations all the protein is free and there is no DNA!
    if(any(p$F.init[p$shapeMat] == 0)){
        p$F.init[p$F.init == 0 & p$shapeMat] = p$I.init[p$F.init == 0 & p$shapeMat]
    }	
    
    p$C.init = p$I.init - p$F.init
    
    p
}

#' Old version of initRates that is using the post-bleach image initialisation
#' 
#' This is the version we are actually using in the paper... 
#' 
#' @param p parameters list
#' @param residence.time
#' @param Free
#' @param D
#' @param exact.free if to perform the correction to exactly estimate the proportion of free molecules
initRatesOld = function(p, residence.time, Free, D, exact.free=TRUE){
    cat("Using OLD initRates()\n")
    if(Free == 1){
        # diffusion only model
        p$Free = Free
        p$D = D
        p$k.off = 0
        p$k.on = 0
    } else {
        if(length(residence.time) == 1){
            p$k.off = 1/residence.time
        } else {
            # multiple values
            target.rois = names(residence.time)
            total.inx = which(target.rois == "")
            if(length(total.inx) != 1){
                stop("Residence time needs to contain exactly one value for global residence time")
            }
            k.off = (p$I.prop > 0) * (1/residence.time[total.inx])
            # now add the exceptions
            for(roi.name in setdiff(target.rois, "")){
                if(!(roi.name %in% names(p$rois)))
                    stop(paste("ROI with name '", roi.name, "' not found", sep=""))
                k.off[ p$rois[[roi.name]] ] = 1/residence.time[roi.name]
            }
            
            p$k.off = k.off
        }
        p$Free = Free
        p$D = D
        
        p = initRatesOldFit(p)
        
        if(exact.free){
            cat("Performing the Free rate correction\n")
            # This assumes `p` is available in the parent env (ie this function)
            opt = optimize(function(scale.k.on){
                p1 = initRatesOldFit(p, scale.k.on)
                
                F.fitted = sum(p1$F.init) / sum(p1$F.init+p1$C.init)
                
                (F.fitted - p1$Free)^2
                
            }, c(1e-5, 2), tol=1e-4)
            
            p = initRatesOldFit(p, opt$minimum)
        }

    }
    
    p
}

#' This version uses initialisation with two species of bound molecules
#' based on Old version of initRates that is using the post-bleach image initialisation
#'
#' @param p the global parameters for fitting
#' @param residence.time the main residence time, can have multiple values within different ROIs
#' @param Free proportion of free molecules
#' @param residence.time.c2 the residence time of second species of molecules
#' @param prop.c2 the proportion (relative to total) of number of molecules in C2 state
initRatesC2 = function(p, residence.time, Free, D, residence.time.c2, prop.c2){
    cat("Using initRatesC2()\n")
    if(Free == 1){
        stop("This function is not appopriate for initialising diffusion-only models. Please use initRatesOld()")
    } else {
        if(length(residence.time) == 1){
            p$k.off.c1 = 1/residence.time
        } else {
            # multiple values
            target.rois = names(residence.time)
            total.inx = which(target.rois == "")
            if(length(total.inx) != 1){
                stop("Residence time needs to contain exactly one value for global residence time")
            }
            k.off = (p$I.prop > 0) * (1/residence.time[total.inx])
            # now add the exceptions
            for(roi.name in setdiff(target.rois, "")){
                if(!(roi.name %in% names(p$rois)))
                    stop(paste("ROI with name '", roi.name, "' not found", sep=""))
                k.off[ p$rois[[roi.name]] ] = 1/residence.time[roi.name]
            }
            
            p$k.off.c1 = k.off
        }
        p$Free = Free
        p$D = D
        
        # second residence time
        p$k.off.c2 = 1/residence.time.c2
        # proportions of various molecule species
        p$prop.c2 = prop.c2
        p$prop.c1 = 1 - Free - prop.c2
        
        # relative proportion of the two populations, needed for correct initialisation of bound molecules
        p$prop.c1c2 = p$prop.c1 / (p$prop.c1 + p$prop.c2)
        
        # varying ON rates
        p$k.on.c1 = p$k.off.c1 / p$Free * (p$I.prop - p$Free) * p$prop.c1c2
        if(any(p$k.on.c1 < 0))
            p$k.on.c1[p$k.on.c1 < 0] = 0
        
        p$k.on.c2 = p$k.off.c2 / p$Free * (p$I.prop - p$Free) * (1-p$prop.c1c2)
        if(any(p$k.on.c2 < 0))
            p$k.on.c2[p$k.on.c2  < 0] = 0
        
        
        # K_b (nod K_d)
        p$kb.c1 = p$k.on.c1 / p$k.off.c1
        p$kb.c2 = p$k.on.c2 / p$k.off.c2
        
        # initial number of Free molecules
        p$F.init = 	p$I.init / (1 + p$kb.c1 + p$kb.c2)
        p$F.init[!p$shapeMat] = 0
        if(any(!is.finite(p$F.init)))
            p$F.init[!is.finite(p$F.init)] = 0
        
        # at these locations all the protein is free and there is no DNA!
        if(any(p$F.init[p$shapeMat] == 0)){
            p$F.init[p$F.init == 0 & p$shapeMat] = p$I.init[p$F.init == 0 & p$shapeMat]
        }	
        
        # directly follows for C1 + C2 + F = I.init
        p$C1.init = p$I.init - p$F.init - p$kb.c2 * p$F.init
        p$C2.init = p$I.init - p$F.init - p$kb.c1 * p$F.init
        
        if(any(p$C1.init < 0))
            p$C1.init[p$C1.init<0] = 0
        
        if(any(p$C2.init < 0))
            p$C2.init[p$C2.init<0] = 0
    }
    
    p
}

#' Initialise the parameters to be uniform over the whole nucleus
#' 
#' @param p the parameters list
#' @param residence.time the res.time to use
#' @param Free percent free
#' @param D diffusion
uniformInit = function(p, residence.time, Free, D){
    cat("Using uniformInit()\n")
    if(Free == 1){
        stop("This function is not appopriate for initialising diffusion-only models. Please use initRatesOld()")
    }
    
    # everything is uniform
    p$D = D
    p$Free = Free
    p$k.off = (1 / residence.time) * p$shapeMat
    p$k.on = p$k.off * (1-Free)/Free
    
    # bleaching based on actual pre/post images
    bleach = p$post.image.sim / p$pre.image.sim
    bleach[!p$shapeMat] = 0
    
    total.sum = sum(p$shapeMat)
    
    ss = steadyState(total.sum, p$shapeMat, p$k.on, p$k.off)
    
    # init value based on the image
    p$F.init = ss$F * bleach
    p$C.init = ss$C * bleach
    
    p$F.init[!p$shapeMat] = 0
    p$C.init[!p$shapeMat] = 0
    
    # rescale so the total is sum(p$shapeMat) as used in the steady state calculation
    scale = total.sum / sum(p$F.init + p$C.init)
    
    p$F.init = p$F.init * scale
    p$C.init = p$C.init * scale
    
    stopifnot(abs(sum(p$F.init+p$C.init) - total.sum) < 1e-8)
    
    p
}

#' Initialise the parameters to be uniform over the whole nucleus
#' 
#' C2 implementation with two species of bound molecules 
#' 
#' @param p the global parameters for fitting
#' @param residence.time the main residence time, can have multiple values within different ROIs
#' @param Free proportion of free molecules
#' @peram D diffusion contact
#' @param residence.time.c2 the residence time of second species of molecules
#' @param prop.c2 the proportion (relative to total) of number of molecules in C2 state
uniformInitC2 = function(p, residence.time, Free, D, residence.time.c2, prop.c2){
    cat("Using uniformInitC2()\n")
    if(Free == 1){
        stop("This function is not appopriate for initialising diffusion-only models. Please use initRatesOld()")
    }
    
    # everything is uniform
    p$D = D
    p$Free = Free
    p$k.off.c1 = (1 / residence.time) * p$shapeMat
    p$k.off.c2 = (1 / residence.time.c2) * p$shapeMat
    p$k.on.c1 = p$k.off.c1 * (1-Free-prop.c2)/Free
    p$k.on.c2 = p$k.off.c2 * prop.c2/Free
    
    # bleaching based on actual pre/post images
    bleach = p$post.image.sim / p$pre.image.sim
    bleach[!p$shapeMat] = 0
    
    total.sum = sum(p$shapeMat)
    
    ss = steadyStateC2(total.sum, p$shapeMat, p$k.on.c1, p$k.on.c2, p$k.off.c1, p$k.off.c2)
    
    # init value based on the image
    p$F.init = ss$F * bleach
    p$C1.init = ss$C1 * bleach
    p$C2.init = ss$C2 * bleach
    
    p$F.init[!p$shapeMat] = 0
    p$C1.init[!p$shapeMat] = 0
    p$C2.init[!p$shapeMat] = 0
    
    # rescale so the total is sum(p$shapeMat) as used in the steady state calculation
    scale = total.sum / sum(p$F.init + p$C1.init + p$C2.init)
    
    p$F.init = p$F.init * scale
    p$C1.init = p$C1.init * scale
    p$C2.init = p$C2.init * scale
    
    stopifnot(abs(sum(p$F.init+p$C1.init+p$C2.init) - total.sum) < 1e-8)
    
    p
}



