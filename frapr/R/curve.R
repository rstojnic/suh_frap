# Recovery curve extraction and fitting

#' Calculate the steady state with a possibly uneven distribution of k_on and k_off
#' @param total Total fluorescence (number)
#' @param shapeMat a TRUE/FALSE matrix showing which pixels are inside the nucleus
#' @param k.on the K_on rate, either a single value or a matrix
#' @param k.off the K_off rate, either a signel value of a matrix
#' @return a list of:
#'         - FL - predicted steady state fluorescence
#'         - F - free molecules
#'         - C - bound molecules
steadyState = function(total, shapeMat, k.on, k.off){
    if(length(k.on)==1 && k.on == 0 && length(k.off)==1 && k.off == 0){
        # Diffusion only model
        F.ss = shapeMat * total/sum(shapeMat)
        return(list("FL"=F.ss, "F"=F.ss, "C"=shapeMat*0))
    }
    
    if(!is.matrix(k.on)){
        k.on = shapeMat * k.on
    }
    if(!is.matrix(k.off)){
        k.off = shapeMat * k.off
    }
    
    # get steady state values, from Supp eq (16) from Beaudouin (2006)
    F.ss = total / sum(1 + (k.on/k.off)[shapeMat])
    F.ss = shapeMat * F.ss
    
    C.ss = F.ss * k.on / k.off
    if(any(is.nan(C.ss))){
        C.ss[is.nan(C.ss)] = 0
    }
    FL.ss = F.ss + C.ss
    
    list("FL"=FL.ss, "F"=F.ss, "C"=C.ss)
}

#' Calculate the steady state with a possibly uneven distribution of k_on and k_off
#' C2 implementation: using two species of bound molecules
#' 
#' @param total Total fluorescence (number)
#' @param shapeMat a TRUE/FALSE matrix showing which pixels are inside the nucleus
#' @param k.on.c1 the K_on rate, either a single value or a matrix (first species)
#' @param k.on.c2 the K_on rate, either a single value or a matrix (second species)
#' @param k.off.c1 the K_off rate, either a signel value of a matrix (first species)
#' @param k.off.c2 the K_off rate, either a signel value of a matrix (second species)
#' @return a list of:
#'         - FL - predicted steady state fluorescence
#'         - F - free molecules
#'         - C1 - first species of bound molecules
#'         - C2 - second species of bound molecules
steadyStateC2 = function(total, shapeMat, k.on.c1, k.on.c2, k.off.c1, k.off.c2){
    if(length(k.on.c1)==1 && k.on.c1 == 0 && length(k.off.c1)==1 && k.off.c1 == 0){
        # Diffusion only model
        F.ss = shapeMat * total/sum(shapeMat)
        return(list("FL"=F.ss, "F"=F.ss, "C1"=shapeMat*0, "C2"=shapeMat*0))
    }
    
    if(!is.matrix(k.on.c1)){
        k.on.c1 = shapeMat * k.on.c1
    }
    if(!is.matrix(k.off.c1)){
        k.off.c1 = shapeMat * k.off.c1
    }
    if(!is.matrix(k.on.c2)){
        k.on.c2 = shapeMat * k.on.c2
    }
    if(!is.matrix(k.off.c2)){
        k.off.c2 = shapeMat * k.off.c2
    }
    
    # get steady state values, F + C1 + C2 = total
    F.ss = total / sum(1 + (k.on.c1/k.off.c1)[shapeMat] + (k.on.c2/k.off.c2)[shapeMat])
    F.ss = shapeMat * F.ss
    
    C1.ss = F.ss * k.on.c1 / k.off.c1
    C2.ss = F.ss * k.on.c2 / k.off.c2
    if(any(is.nan(C1.ss))){
        C1.ss[is.nan(C1.ss)] = 0
    }
    if(any(is.nan(C2.ss))){
        C2.ss[is.nan(C2.ss)] = 0
    }
    FL.ss = F.ss + C1.ss + C2.ss
    
    list("FL"=FL.ss, "F"=F.ss, "C1"=C1.ss, "C2"=C2.ss)
}

#' Extract the recovery curve from data 
#'
#' @param sim simulation results
#' @param roi region of interest. Defaults to the bleach mask. 
#' @return the recovery curve
recoveryCurve = function(sim, roi=NULL){
	FL = sim$FL
	p = sim$params
	
	if(is.null(roi))
		roi = p$bleachMask
	
	if("k.on.c1" %in% names(p) & "k.on.c2" %in% names(p)){
	    # this is the C2 implementation with two types of bound molecules
	    
	    # there can be only two values for total (depending if there are values before/after or only after)		
	    stopifnot(length(unique(round(apply(FL, 3, sum),3))) <= 2)
	    
	    # note: this code assume k.on is a matrix
	    if(!is.matrix(p$k.on.c1))
	        p$k.on.c1 = p$k.on.c1 * (p$shapeMat+0)
	    if(!is.matrix(p$k.on.c2))
	        p$k.on.c2 = p$k.on.c2 * (p$shapeMat+0)
	    
	    # since it might change before/after bleach, use the total for each data point
	    curve = rep(0, dim(FL)[3])
	    for(i in 1:dim(FL)[3]){
	        total = sum(FL[,,i])
	        ss = steadyStateC2(total, p$shapeMat, p$k.on.c1, p$k.on.c2, p$k.off.c1, p$k.off.c2)
	        
	        # calculate recovery relative to ss value
	        curve[i] = mean(FL[,,i][roi]) / mean(ss$FL[roi])
	    }
	} else if(all(p$k.off == 0)){	
	    curve = apply(FL, 3, function(x) mean(x[roi]) / mean(x[p$shapeMat]))    
	} else {
		# there can be only two values for total (depending if there are values before/after or only after)		
		stopifnot(length(unique(round(apply(FL, 3, sum),3))) <= 2)
		
		# note: this code assume k.on is a matrix
		if(!is.matrix(p$k.on))
			p$k.on = p$k.on * (p$shapeMat+0)

		
		# since it might change before/after bleach, use the total for each data point
		curve = rep(0, dim(FL)[3])
		for(i in 1:dim(FL)[3]){
			total = sum(FL[,,i])
			ss = steadyState(total, p$shapeMat, p$k.on, p$k.off)
				
			# calculate recovery relative to ss value
			curve[i] = mean(FL[,,i][roi]) / mean(ss$FL[roi])
			#curve = sapply(FL, 3, function(x) mean(FL[roi]) / mean(FL.ss[roi]))
		}
	}
	
	curve
}

#' Plot the recovery curve for an experiment
#'
#' @param sim Simulation results
#' @param ... parameters to pass to plot
plotRecoveryCurve = function(sim, ...){
	if(!("curve" %in% names(sim))){
		stop("No recovery curve found in simulation results. Did you run 'recoveryCurve'?")
	}
	
	plot(sim$t, sim$curve, ylim=c(0,max(1, sim$curve, na.rm=TRUE)), xlab="Time (s)", ylab="Recovery", ...)
	lines(sim$t, sim$curve)
	
	abline(h=1, lty=2)
}
