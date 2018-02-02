# Reporting functions

#' Force the bleach depth to the exact measured value
#' 
#' @param img set of images and status
#' @param p the set the simulation parameters to modify
#' @return modified parameters
roiBleachDepthCorrection = function(img, p){
    cat("Performing bleach depth correction\n")
    target.depth.inx = unique(apply(img$stats$rec, 2, function(x) which(is.na(x))+1))
    stopifnot(length(target.depth.inx)==1)
    target.depth = img$stats$rec[target.depth.inx,]
    
    c2 = "C1.init" %in% names(p)
    
    if(c2){
        FL = p$F.init + p$C1.init + p$C2.init
        total.sum = sum(FL[p$shapeMat])
        ss = steadyStateC2(total.sum, p$shapeMat, p$k.on.c1, p$k.on.c2, p$k.off.c1, p$k.off.c2)
    } else {
        FL = p$F.init + p$C.init
        total.sum = sum(FL[p$shapeMat])
        ss = steadyState(total.sum, p$shapeMat, p$k.on, p$k.off)
    }
    
    # we need to preserve the total, so we need to distribute any offsets to the rest of the nucleus
    
    total.diff = 0
    mask = p$shapeMat
    # loop over ROIs
    for(i in 1:length(target.depth)){
        roi = p$rois[[i]]
        mean.roi = mean(FL[roi])/mean(ss$FL[roi])
        s = mean.roi / target.depth[i]
        # record where made changes so we can offset so the total is the same
        mask = mask & !roi
        if(c2){
            p$F.init[roi] = p$F.init[roi] / s
            p$C1.init[roi] = p$C1.init[roi] / s
            p$C2.init[roi] = p$C2.init[roi] / s
        } else {
            p$F.init[roi] = p$F.init[roi] / s
            p$C.init[roi] = p$C.init[roi] / s
        }
    }
    
    if(c2){
        new.FL = p$F.init + p$C1.init + p$C2.init
    } else {
        new.FL = p$F.init + p$C.init
    }
    
    total.s = (total.sum - sum(new.FL[p$shapeMat & !mask])) / sum(new.FL[mask])
    
    # offset so the total is the same
    if(c2){
        p$F.init[mask] = p$F.init[mask] * total.s
        p$C1.init[mask] = p$C1.init[mask] * total.s
        p$C2.init[mask] = p$C2.init[mask] * total.s
        final.FL = p$F.init + p$C1.init + p$C2.init
    } else {
        p$F.init[mask] = p$F.init[mask] * total.s
        p$C.init[mask] = p$C.init[mask] * total.s
        final.FL = p$F.init + p$C.init
    }
    
    stopifnot(abs( sum(final.FL[p$shapeMat]) - total.sum ) < 1e-8)
    
    p
}

#' Evaluate parameter values
#'
#' @param img a set of images
#' @param p a set of simulation parameters
#' @param res.time.range the values to test for residence time
#' @param Free.range the values to test for the Free proportion
#' @param D.range the value of diffusion constants to check
#' @return the list of average error rates and best parameters
#' @param first.ROI.only if to only fit the first ROI
#' @param use.old.init.rates if to use the old post-bleach image initialisation
#' @param use.bleaching.init if to use the simulation of bleaching to initialise the image
#' @param bleaching.init.version the version of the bleaching initialisation algorithm
#' @param roi.bleach.depth.correction force the bleach depth to the exact measured value
evalParams = function(img, p, res.time.range, Free.range, D.range, first.ROI.only=TRUE, use.old.init.rates=FALSE,
                      use.bleaching.init=FALSE, bleaching.init.version=1, roi.bleach.depth.correction=TRUE,
                      use.uniform.init=FALSE) {
	err = array(0, dim=c(length(res.time.range), length(Free.range), length(D.range)))
	dimnames(err) = list(res.time.range, Free.range, D.range)

	sim = list()

	for(i in 1:dim(err)[1]){
		cat("Processing with residence time", i, "/", dim(err)[1], "\n")
		sim[[i]] = list()
		for(j in 1:dim(err)[2]){
			cat("=> Processing with Free range", j, "/", dim(err)[2], "\n")
			sim[[i]][[j]] = list()
			for(k in 1:dim(err)[3]){
			    cat("=> Processing with D range", k, "/", dim(err)[3], "\n")
			    if(use.bleaching.init){
			        cat("Using bleaching initialisation\n")
			        
			        # simulate the bleaching
			        p.sim = initRatesOld(p, residence.time=res.time.range[i], Free=Free.range[j], D=D.range[k])
			        p.sim$hx = p$L.x / p$grid.x
			        p.sim$hy = p$L.y / p$grid.y
			        bleach.sim = simulateCircularBleaching2D(p.sim)
			        
			        # re-calibrate the bleaching depths of F and C based on the bleaching simulations
			        if(bleaching.init.version == 1) {
			            # current implementation
			            final.fl = bleach.sim$final.F + bleach.sim$final.C
			            prop.c = bleach.sim$final.C / final.fl
			            prop.c[is.nan(prop.c)] = 0
			            
			            # re-calibrate the bleaching depths
			            p.cur = p.sim
			            total.FL = p.cur$F.init + p.cur$C.init
			            p.cur$C.init = total.FL * prop.c
			            p.cur$F.init = total.FL * (1-prop.c)
			        } else {
			            # alternative implementation: just use the F and C values
			            p.cur = p.sim
			            p.cur$F.init = bleach.sim$final.F
			            p.cur$C.init = bleach.sim$final.C
			        }
			        p.cur$bleach.sim = bleach.sim
			        
			        p.cur$bleachRate = bleach.sim$bleachRate
			    } else if(use.uniform.init){
			        p.cur = uniformInit(p, residence.time=res.time.range[i], Free=Free.range[j], D=D.range[k])
			    } else if(use.old.init.rates){
			        p.cur = initRatesOld(p, residence.time=res.time.range[i], Free=Free.range[j], D=D.range[k])
			    } else {
				    p.cur = initRates(p, residence.time=res.time.range[i], Free=Free.range[j], D=D.range[k])
			    }
			    
			    if(roi.bleach.depth.correction)
			        p.cur = roiBleachDepthCorrection(img, p.cur)
			    
				cat("Starting simulation\n")
				total.time = system.time({sim.mean = simulateMeans2D(p.cur)})
				cat("Simulation finished in ", total.time, "\n")

				# fit only the recovery of the first ROI
				rec = img$stats$rec[img$post.frame:(nrow(img$stats$rec)-1),]
				curve = sim.mean$curve
				if(first.ROI.only){
				    err[i,j,k] = mean((rec[,1]-curve[,1])^2)
				} else {
				    err[i,j,k] = mean((rec-curve)^2)
				}
				
				# record a stripped down object
				sim.mean$FL = NULL
				sim.mean$F = NULL
				sim.mean$C = NULL
				sim[[i]][[j]][[k]] = sim.mean
				
			}
		}
	}

	# get the best value and plot
	inx = which(err==min(err), arr.ind=TRUE)
	res.time = res.time.range[inx[1,1]]
	Free = Free.range[inx[1,2]]
	D = D.range[inx[1,3]]
	
	list(sim=sim, err=err, res.time=res.time, Free=Free, D=D)

}

#' Evaluate parameters using two different residence times
#'
#' @param img a set of images
#' @param p a set of simulation parameters
#' @param res.time1.range the values to test for first (global) residence time
#' @param res.time2.range the values to test for second residence time
#' @param res.time.roi roi for the second residence time
#' @param Free.range the values to test for the Free proportion
#' @param D.range the value of diffusion constants to check
#' @return the list of average error rates and best parameters
evalParams2Res = function(img, p, res.time1.range, res.time2.range, res.time.roi, Free.range, D.range, 
                          use.old.init.rates=FALSE, use.bleaching.init=FALSE, bleaching.init.version=1) {
	# init output variable
	err = array(0, dim=c(length(res.time1.range), length(res.time2.range), 
		length(Free.range), length(D.range)))
	dimnames(err) = list(res.time1.range, res.time2.range, Free.range, D.range)

	sim = list()

	for(i in 1:dim(err)[1]){
		cat("Processing with residence time1", i, "/", dim(err)[1], "\n")
		sim[[i]] = list()
		for(j in 1:dim(err)[2]){
			cat("=> Processing with residence time2", j, "/", dim(err)[2], "\n")
			sim[[i]][[j]] = list()
			for(k in 1:dim(err)[3]){
				cat("=> => Processing with Free range", k, "/", dim(err)[3], "\n")
				sim[[i]][[j]][[k]] = list()
				for(l in 1:dim(err)[4]){
					res.time = c(res.time1.range[i], res.time2.range[j])
					names(res.time) = c("", res.time.roi)
					if(use.bleaching.init){
					    cat("Using bleaching initialisation\n")
					    
					    # simulate the bleaching
					    p.sim = initRatesOld(p, residence.time=res.time, Free=Free.range[k], D=D.range[l])
					    p.sim$hx = p$L.x / p$grid.x
					    p.sim$hy = p$L.y / p$grid.y
					    bleach.sim = simulateCircularBleaching2D(p.sim)
					    
					    # re-calibrate the bleaching depths of F and C based on the bleaching simulations
					    if(bleaching.init.version == 1) {
					        # current implementation
					        final.fl = bleach.sim$final.F + bleach.sim$final.C
					        prop.c = bleach.sim$final.C / final.fl
					        prop.c[is.nan(prop.c)] = 0
					        
					        # re-calibrate the bleaching depths
					        p.cur = p.sim
					        total.FL = p.cur$F.init + p.cur$C.init
					        p.cur$C.init = total.FL * prop.c
					        p.cur$F.init = total.FL * (1-prop.c)
					    } else {
					        # alternative implementation: just use the F and C values
					        p.cur = p.sim
					        p.cur$F.init = bleach.sim$final.F
					        p.cur$C.init = bleach.sim$final.C
					    }
					    p.cur$bleach.sim = bleach.sim
					    
					    p.cur$bleachRate = bleach.sim$bleachRate
					} else if(use.old.init.rates){
					    p.cur = initRatesOld(p, residence.time=res.time, Free=Free.range[k], D=D.range[l])
					} else {
					    p.cur = initRates(p, residence.time=res.time, Free=Free.range[k], D=D.range[l])
					}
					
					cat("Starting simulation\n")
					total.time = system.time({sim.mean = simulateMeans2D(p.cur)})
					cat("Simulation finished in ", total.time, "\n")

					rec = img$stats$rec[img$post.frame:(nrow(img$stats$rec)-1),]
					curve = sim.mean$curve
					err[i,j,k,l] = mean((rec-curve)^2)
				
					# record a stripped down object
					sim.mean$FL = NULL
					sim.mean$F = NULL
					sim.mean$C = NULL
					sim[[i]][[j]][[k]][[l]] = sim.mean
				
				}
			}
		}
	}

	# get the best value and plot
	inx = which(err==min(err), arr.ind=TRUE)
	res.time1 = res.time1.range[inx[1,1]]
	res.time2 = res.time2.range[inx[1,2]]
	Free = Free.range[inx[1,3]]
	D = D.range[inx[1,4]]
	
	list(sim=sim, err=err, res.time1=res.time1, res.time2=res.time2, Free=Free, D=D)

}

#' Evaluate parameter values, two species of bound molecules
#'
#' @param img a set of images
#' @param p a set of simulation parameters
#' @param res.time.c1.range the values to test for residence time (first bound species)
#' @param res.time.c2.range the values to test for residence time (second bound species)
#' @param prop.c2.range range for propotion of second bound species (as propotional of total)
#' @param Free.range the values to test for the Free proportion
#' @param D.range the value of diffusion constants to check
#' @return the list of average error rates and best parameters
#' @param first.ROI.only if to only fit the first ROI
#' @param use.old.init.rates if to use the old post-bleach image initialisation
#' @param use.bleaching.init if to use the simulation of bleaching to initialise the image
#' @param bleaching.init.version the version of the bleaching initialisation algorithm
evalParamsC2 = function(img, p, res.time.c1.range, res.time.c2.range, prop.c2.range, 
                        Free.range, D.range, first.ROI.only=TRUE,
                        use.bleaching.init=FALSE, bleaching.init.version=1,
                        roi.bleach.depth.correction=TRUE, use.uniform.init=FALSE) {
    err = array(0, dim=c(length(res.time.c2.range), 
                         length(prop.c2.range), 
                         length(res.time.c1.range),
                         length(Free.range), 
                         length(D.range)))
    dimnames(err) = list(res.time.c2.range, 
                         prop.c2.range,
                         res.time.c1.range,
                         Free.range, 
                         D.range)
    
    sim = list()
    
    for(i1 in 1:dim(err)[1]){
        sim[[i1]] = list()
        for(i2 in 1:dim(err)[2]){
            sim[[i1]][[i2]] = list()
            for(i3 in 1:dim(err)[3]){
                sim[[i1]][[i2]][[i3]] = list()
                for(i4 in 1:dim(err)[4]){
                    sim[[i1]][[i2]][[i3]][[i4]] = list()
                    for(i5 in 1:dim(err)[5]){
                        cat("res.time.c2", i1, "/", dim(err)[1], 
                            ", prop.c2", i2, "/", dim(err)[2],
                            ", res.time.c1", i3, "/", dim(err)[3],
                            ", free", i4, "/", dim(err)[4],
                            ", D", i5, "/", dim(err)[5], "\n")
                        
                        if(prop.c2.range[i2] >= (1-Free.range[i4])){
                            # This combination of parameters doesn't make sense
                            sim[[i1]][[i2]][[i3]][[i4]][[i5]] = NULL 
                            err[i1,i2,i3,i4,i5] = NA
                        } else if(use.bleaching.init){
                            stop("Bleaching init not implemented yet")
                        } else {
                            if(use.uniform.init){
                                p.cur = uniformInitC2(p, residence.time=res.time.c1.range[i3], 
                                                      Free=Free.range[i4], 
                                                      D=D.range[i5],
                                                      residence.time.c2=res.time.c2.range[i1],
                                                      prop.c2=prop.c2.range[i2])
                            } else {
                                p.cur = initRatesC2(p, 
                                                    residence.time=res.time.c1.range[i3], 
                                                    Free=Free.range[i4], 
                                                    D=D.range[i5],
                                                    residence.time.c2=res.time.c2.range[i1],
                                                    prop.c2=prop.c2.range[i2])
                            }
                
                            if(roi.bleach.depth.correction)
                                p.cur = roiBleachDepthCorrection(img, p.cur)
                                        
                            cat("Starting simulation\n")
                            total.time = system.time({sim.mean = simulateMeans2D_C2(p.cur)})
                            cat("Simulation finished in ", total.time, "\n")
                            
                            # fit only the recovery of the first ROI
                            rec = img$stats$rec[img$post.frame:(nrow(img$stats$rec)-1),]
                            curve = sim.mean$curve
                            if(first.ROI.only){
                                err[i1,i2,i3,i4,i5] = mean((rec[,1]-curve[,1])^2)
                            } else {
                                err[i1,i2,i3,i4,i5] = mean((rec-curve)^2)
                            }
                            
                            # record a stripped down object
                            sim.mean$FL = NULL
                            sim.mean$F = NULL
                            sim.mean$C1 = NULL
                            sim.mean$C2 = NULL
                            sim[[i1]][[i2]][[i3]][[i4]][[i5]] = sim.mean
                        }
                    }
                }
            }
        }
    }
    
    # get the best value and plot
    inx = which(err==min(err, na.rm=TRUE), arr.ind=TRUE)
    res.time.c2 = res.time.c2.range[inx[1,1]]
    prop.c2 = prop.c2.range[inx[1,2]]
    res.time.c1 = res.time.c1.range[inx[1,3]]
    Free = Free.range[inx[1,4]]
    D = D.range[inx[1,5]]
    
    list(sim=sim, err=err, res.time.c2=res.time.c2, prop.c2=prop.c2, res.time.c1=res.time.c1, Free=Free, D=D)
    
}