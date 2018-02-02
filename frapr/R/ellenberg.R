# Ellenberg model implementation

#' Ellenberg 2D Integration function
#' 
#' @param t current time
#' @param state the state vector
#' @param p parameters
.ellenberg2D_full_int = function(t, state, p){
	# the last 0 is for the case when the value should be omitted
	state.mat = matrix(state, ncol=2)
	F.clean = state.mat[,1]
	F = c(state.mat[,1], 0)

	C = state.mat[,2]

	dF = p$fromL*F[p$inxL] + p$fromR*F[p$inxR] - p$leaveX * F.clean +
		 p$fromU*F[p$inxU] + p$fromD*F[p$inxD] - p$leaveY * F.clean - 
		 p$k.on * F.clean + p$k.off * C

	dC = p$k.on * F.clean - p$k.off * C
		  
	list(c(as.vector(dF), as.vector(dC)))
}

#' Ellenberg 2D Integration function, C2 version with two species of bound molecules
#' 
#' @param t current time
#' @param state the state vector
#' @param p parameters
.ellenberg2D_C2_full_int = function(t, state, p){
    # the last 0 is for the case when the value should be omitted
    state.mat = matrix(state, ncol=3)
    F.clean = state.mat[,1]
    F = c(state.mat[,1], 0)
    
    C1 = state.mat[,2]
    C2 = state.mat[,3]
    
    dF = p$fromL*F[p$inxL] + p$fromR*F[p$inxR] - p$leaveX * F.clean +
        p$fromU*F[p$inxU] + p$fromD*F[p$inxD] - p$leaveY * F.clean - 
        p$k.on.c1 * F.clean + p$k.off.c1 * C1 -
        p$k.on.c2 * F.clean + p$k.off.c2 * C2
    
    dC1 = p$k.on.c1 * F.clean - p$k.off.c1 * C1
    dC2 = p$k.on.c2 * F.clean - p$k.off.c2 * C2
    
    list(c(as.vector(dF), as.vector(dC1), as.vector(dC2)))
}


#' Ellenberg 2D Integration function for bleaching
#' 
#' Both C and F molecules are bleached with the rate of p$bleachRate (matrix of rates)
#' 
#' @param t current time
#' @param state the state vector
#' @param p parameters
.ellenberg2D_bleaching_full_int = function(t, state, p){
    # the last 0 is for the case when the value should be omitted
    state.mat = matrix(state, ncol=2)
    F.clean = state.mat[,1]
    F = c(state.mat[,1], 0)
    
    C = state.mat[,2]
    
    dF = p$fromL*F[p$inxL] + p$fromR*F[p$inxR] - p$leaveX * F.clean +
        p$fromU*F[p$inxU] + p$fromD*F[p$inxD] - p$leaveY * F.clean - 
        p$k.on * F.clean + p$k.off * C - p$bleachRate * F.clean
    
    dC = p$k.on * F.clean - p$k.off * C - p$bleachRate * C
    
    list(c(as.vector(dF), as.vector(dC)))
}


#' Simulate the Ellenberg model using a set of parameters
#'
#' @param params a set of parameters
#' @return sim object
simulateMeans2D = function(params){
	p = params
	.checkParams(p)
	
	shapeMat = p$shapeMat
	if("D.ncl" %in% names(p)){
		D.ncl = p$D.ncl
	} else {
		D.ncl = p$D
	}
	if("nclMask" %in% names(p)){
		nclMask = p$nclMask
	} else {
		nclMask = shapeMat
		nclMask[,] = FALSE
	}
	
	if(!("L.x" %in% names(p))){
		p$L.x = p$L		
	}
	if(!("L.y" %in% names(p))){
		p$L.y = p$L		
	}
	
	hx = p$L.x / p$grid.x
	hy = p$L.y / p$grid.y

	# diffusion: 10um^2/s
	dx = p$D / hx^2
	dy = p$D / hy^2
	
	dx.ncl = p$D.ncl / hx^2
	dy.ncl = p$D.ncl / hy^2

	# Matrices of movement INTO tile (i,j) - this is the oppositve of what we have
	# in our Gillepsie simulations (which is out of)
	fromL = fromR = matrix(dx, nrow=p$grid.y, ncol=p$grid.x)
	fromU = fromD = matrix(dy, nrow=p$grid.y, ncol=p$grid.x)

	# indicies of pixels that are left and right
	# the default is one out of range, which will be pre-filled with 0
	inxL = inxR = inxU = inxD = matrix(p$grid.x * p$grid.y + 1, nrow=p$grid.y, ncol=p$grid.x)

	#' Convert 2D to 1D coordinates
	#'
	#' @param i row coordinate
	#' @param j col coordinate
	#' @param nr number of rows
	to1D = function(i, j, nr=nrow(shapeMat)){
		i + (j-1)*nr
	}

	# create the diffusion matrices
	for(i in 1:nrow(shapeMat)){
		for(j in 1:ncol(shapeMat)){
			if(!shapeMat[i,j]){
				fromL[i,j] = fromR[i,j] = fromD[i,j] = fromU[i,j] = 0
				next
			}
			
			###########################################
			# cannot come from L
			if(j == 1 || !shapeMat[i,j-1]){
				fromL[i,j] = 0
			} else {
				inxL[i,j] = to1D(i,j-1)
			}
			
			# nucleolus diffusion from L
			if(j != 1 && nclMask[i,j] && nclMask[i,j-1])
				fromL[i,j] = dx.ncl
			
			###########################################
			# cannot come from R
			if(j == ncol(shapeMat) || !shapeMat[i,j+1]){
				fromR[i,j] = 0	
			} else {
				inxR[i,j] = to1D(i,j+1)
			}
			
			# nucleolus diffusion from R
			if(j != ncol(shapeMat) && nclMask[i,j] && nclMask[i,j+1])
				fromR[i,j] = dx.ncl
			
			###########################################
			# cannot come from U
			if(i == 1 || !shapeMat[i-1,j]){
				fromU[i,j] = 0
			} else {
				inxU[i,j] = to1D(i-1,j)
			}
			
			# nucleolus diffusion from U
			if(i != 1 && nclMask[i,j] && nclMask[i-1,j])
				fromU[i,j] = dy.ncl
			
			############################################
			# cannot come from D
			if(i == nrow(shapeMat) || !shapeMat[i+1,j]){
				fromD[i,j] = 0	
			} else {
				inxD[i,j] = to1D(i+1,j)
			}			
			
			# nucleolus diffusion from D
			if(i != nrow(shapeMat) && nclMask[i,j] && nclMask[i+1,j])
				fromD[i,j] = dy.ncl

		}
	}
	
	# zero flux parameters
	p$leaveX = fromL + fromR
	p$leaveY = fromU + fromD

	# put into parameters
	p$fromL = fromL
	p$fromR = fromR
	p$fromU = fromU
	p$fromD = fromD
	p$inxL = inxL
	p$inxR = inxR
	p$inxU = inxU
	p$inxD = inxD
	
	if(all(c("F.init", "C.init") %in% names(p))){
		# if specified use them
		F.init = p$F.init
		C.init = p$C.init
	} else {
	    cat("Generating F.init and C.init values\n")
		# initial state
		F.eq = p$k.off / (p$k.on + p$k.off)
		C.eq = p$k.on / (p$k.on + p$k.off)

		# this happens when k.on or k.off are 0... 	
		if(any(!is.finite(F.eq))){
			if(any(!is.finite(F.eq) & p$shapeMat))
				F.eq[!is.finite(F.eq) & p$shapeMat] = 1
			if(any(!is.finite(F.eq) & !p$shapeMat))
				F.eq[!is.finite(F.eq) & !p$shapeMat] = 0
		}
		if(any(!is.finite(C.eq))){
			C.eq[!is.finite(C.eq)] = 0
		}
	
		F.init = (shapeMat + 0) * F.eq
		F.init[p$bleachMask] = F.init[p$bleachMask] * p$bleach.depth

		C.init = (shapeMat + 0) * C.eq
		C.init[p$bleachMask] = C.init[p$bleachMask] * p$bleach.depth
		
		p$F.init = F.init
		p$C.init = C.init
	}

	F.init.vect = as.vector(F.init)
	C.init.vect = as.vector(C.init)

	init.vect = c(F.init.vect, C.init.vect)

	# time
	if("t" %in% names(p)){
		t = p$t[-length(p$t)]
	} else {
		t = seq(0, p$max.time-1/p$fps-p$bleach.time, 1/p$fps)
	}
	#out = ode(init.vect, t, cppEllenberg2D_full_int, p)
	out = ode(init.vect, t, .ellenberg2D_full_int, p)

	# need to convert 'out' into a sim object
	sim = list()
	sim$F.init = F.init
	sim$C.init = C.init
	sim$FL = array(0, dim=c(p$grid.y, p$grid.x, length(t)))
	sim$F = sim$C = array(0, dim=c(p$grid.y, p$grid.x, length(t)))
	for(i in 1:nrow(out)){
		vect = out[i,-1]
		m = matrix(vect, ncol=2)
		sim$F[,,i] = matrix(m[,1], ncol=p$grid.x)
		sim$C[,,i] = matrix(m[,2], ncol=p$grid.x)
		sim$FL[,,i] = sim$F[,,i] + sim$C[,,i]
	}
	sim$t = t
	sim$params = p
	if("rois" %in% names(p)){
		curve = matrix(0, nrow=dim(sim$F)[3], ncol=length(p$rois))
		colnames(curve) = names(p$rois)
		for(i in 1:length(p$rois)){
			curve[,i] = recoveryCurve(sim, p$rois[[i]])
		}
		sim$curve = curve
	} else {
		sim$curve = recoveryCurve(sim)
	}
	
	sim
}

#' Create a circular bleach rate profile for the bleaching simulation
#' 
#' The base rate will be set to 1, so it can be easily modified later on
#' 
#' @param params list of parameters, including the bleaching position
#' @return the bleach rate profile
createCircularBleachingRatelessProfile = function(params){
    bleachRate = matrix(0, nrow=nrow(p$shapeMat), ncol=ncol(p$shapeMat))
    for(i in 1:nrow(bleachRate)){
        for(j in 1:ncol(bleachRate)){
            # distance in real units, note that rows are "y" and cols are "x"
            d = sqrt( ( (i-p$bleach.center[1])*hy) ^2 + ( (j-p$bleach.center[2])*hx)^2) 
            if(d < p$bleach.radius){
                bleachRate[i,j] = 1
            }
        }
    }
    bleachRate
}

#' Create a realistic Gaussian bleaching profile
#' 
#' NOTE: in order to avoid the loss of precision due to scaling, the calculation is
#' done on the full-scale image which is then scaled down. 
#' 
#' @param p the parameters list
createGaussianBleachingRatelessProfile = function(p){
    norm.sd = p$bleach.norm.sd
    
    # use the original full-res ROI to create the bleaching mask
    bleachMask = p$originalBleachMask
    bleach.center = colMeans(which(bleachMask, arr.ind=TRUE))
    
    # this is because in our images the first dimension is y, second is x
    bleach.x = bleach.center[2]
    bleach.y = bleach.center[1]
    
    # because of scaling our resolution is bigger.. 
    hx = p$hx / p$originalScaleBy
    hy = p$hy / p$originalScaleBy
    
    # make the bleaching rate matrix
    bleachRate = matrix(0, nrow=nrow(bleachMask), ncol=ncol(bleachMask))
    
    for(y in 1:nrow(bleachRate)){
        for(x in 1:ncol(bleachRate)){
            d = sqrt( ( (y-bleach.y)*hy)^2 + ( (x-bleach.x)*hx)^2)
            bleachRate[y,x] = dnorm(d, mean=0, sd=norm.sd)
        }
    }
    
    # normalize to 1
    bleachRate = bleachRate / max(bleachRate)
    bleachRate = .scaleDownImage(bleachRate, p$originalScaleBy)
    if("shapeMat" %in% names(p))
        bleachRate[!p$shapeMat] = 0
    
    bleachRate
}

#' Convert the output object from ode simulation to a "sim" object
#' @param out the output of ode integration
#' @param p the parameters used in the simulation
.ode_output_to_sim = function(out, p){
    # need to convert 'out' into a sim object
    sim = list()
    sim$F.init = p$F.init
    sim$C.init = p$C.init
    sim$FL = array(0, dim=c(p$grid.y, p$grid.x, length(p$t)))
    sim$F = sim$C = array(0, dim=c(p$grid.y, p$grid.x, length(p$t)))
    for(i in 1:nrow(out)){
        vect = out[i,-1]
        m = matrix(vect, ncol=2)
        sim$F[,,i] = matrix(m[,1], ncol=p$grid.x)
        sim$C[,,i] = matrix(m[,2], ncol=p$grid.x)
        sim$FL[,,i] = sim$F[,,i] + sim$C[,,i]
    }
    sim$t = p$t
    sim$params = p
    
    sim
    
}

#' This is the actual bit of the code doing the simulation of bleaching
#' It is separated out here because it used to optimise the bleaching depth using optim()
#' 
#' @param p the parameters list
#' @param bleachRate the strength of bleaching (this is the one parameter we are optimising)
#' @param returns the list of bleaching, delay simulations and achieved bleaching depth
.simulateCircularBleaching2D_inner = function(p, bleachRate){
    sims = list()
    
    p.orig = p
    
    # for now use the bleach mask to initialise the bleach mask
    p$bleachRate = p$bleachRateless * bleachRate
    
    # do the actual integration 
    out = ode(p$FC.init.vect, p$t, .ellenberg2D_bleaching_full_int, p)
    sim = .ode_output_to_sim(out, p)
    
    sims = list(sim.bleach = sim)
    last.frame.inx = dim(sim$F)[3]
    final.F = sim$F[,,last.frame.inx]
    final.C = sim$C[,,last.frame.inx]
    final.FL = final.F+final.C
    initial.FL = sim$F[,,1] + sim$C[,,1]
    
    # second simulation of the delay
    if(p$bleach.delay > 0){
        p.delay = p.orig
        # take the values from the last frame to initialise the second simulation
        p.delay$F.init = final.F
        p.delay$C.init = final.C
        
        p.delay$t = c(0, p$bleach.delay)
        
        init.vect.delay = c(as.vector(p.delay$F.init), as.vector(p.delay$C.init))
        out.delay = ode(init.vect.delay, p.delay$t, .ellenberg2D_full_int, p.delay)
        sim.delay = .ode_output_to_sim(out.delay, p.delay)
        
        sims$sim.delay = sim.delay
        final.F = sim.delay$F[,,2]
        final.C = sim.delay$C[,,2]
        final.FL = final.F+final.C
    }
    
    # estimate the bleaching depth
    ss = steadyState(sum(final.FL[p$shapeMat]), p$shapeMat, p$k.on, p$k.off)
    
    sims$bleach.depth = sum(final.FL[p$bleachMask]) / sum(ss$FL[p$bleachMask])
    sims$bleach.fl.prop = sum(final.FL[p$shapeMat]) / sum(initial.FL[p$shapeMat]) 
    sims$bleach.profile = final.FL / ss$FL
    sims$bleach.profile[is.nan(sims$bleach.profile)] = 0
    sims$final.F = final.F
    sims$final.C = final.C
    sims$final.FL = final.FL
    sims$bleachRate = bleachRate
    sims$bleach.norm.sd = p$bleach.norm.sd
    
    sims
}

#' This function simulates a model with constant bleaching
#' 
#' It estimates the best bleaching rate based on the real bleaching depth and simulation of bleaching
#' 
#' @param params list with the simulation parmaeters
#' @param bleachRates min and max values of bleaching rates to try
simulateCircularBleaching2D = function(params, bleachRates=c(0.01, 10)){
    p = params
    .checkParams(p)
    
    shapeMat = p$shapeMat
    if("D.ncl" %in% names(p)){
        D.ncl = p$D.ncl
    } else {
        D.ncl = p$D
    }
    if("nclMask" %in% names(p)){
        nclMask = p$nclMask
    } else {
        nclMask = shapeMat
        nclMask[,] = FALSE
    }
    
    if(!("L.x" %in% names(p))){
        p$L.x = p$L		
    }
    if(!("L.y" %in% names(p))){
        p$L.y = p$L		
    }
    
    hx = p$L.x / p$grid.x
    hy = p$L.y / p$grid.y
    
    # diffusion: 10um^2/s
    dx = p$D / hx^2
    dy = p$D / hy^2
    
    dx.ncl = p$D.ncl / hx^2
    dy.ncl = p$D.ncl / hy^2
    
    # Matrices of movement INTO tile (i,j) - this is the oppositve of what we have
    # in our Gillepsie simulations (which is out of)
    fromL = fromR = matrix(dx, nrow=p$grid.y, ncol=p$grid.x)
    fromU = fromD = matrix(dy, nrow=p$grid.y, ncol=p$grid.x)
    
    # indicies of pixels that are left and right
    # the default is one out of range, which will be pre-filled with 0
    inxL = inxR = inxU = inxD = matrix(p$grid.x * p$grid.y + 1, nrow=p$grid.y, ncol=p$grid.x)
    
    #' Convert 2D to 1D coordinates
    #'
    #' @param i row coordinate
    #' @param j col coordinate
    #' @param nr number of rows
    to1D = function(i, j, nr=nrow(shapeMat)){
        i + (j-1)*nr
    }
    
    # create the diffusion matrices
    for(i in 1:nrow(shapeMat)){
        for(j in 1:ncol(shapeMat)){
            if(!shapeMat[i,j]){
                fromL[i,j] = fromR[i,j] = fromD[i,j] = fromU[i,j] = 0
                next
            }
            
            ###########################################
            # cannot come from L
            if(j == 1 || !shapeMat[i,j-1]){
                fromL[i,j] = 0
            } else {
                inxL[i,j] = to1D(i,j-1)
            }
            
            # nucleolus diffusion from L
            if(j != 1 && nclMask[i,j] && nclMask[i,j-1])
                fromL[i,j] = dx.ncl
            
            ###########################################
            # cannot come from R
            if(j == ncol(shapeMat) || !shapeMat[i,j+1]){
                fromR[i,j] = 0	
            } else {
                inxR[i,j] = to1D(i,j+1)
            }
            
            # nucleolus diffusion from R
            if(j != ncol(shapeMat) && nclMask[i,j] && nclMask[i,j+1])
                fromR[i,j] = dx.ncl
            
            ###########################################
            # cannot come from U
            if(i == 1 || !shapeMat[i-1,j]){
                fromU[i,j] = 0
            } else {
                inxU[i,j] = to1D(i-1,j)
            }
            
            # nucleolus diffusion from U
            if(i != 1 && nclMask[i,j] && nclMask[i-1,j])
                fromU[i,j] = dy.ncl
            
            ############################################
            # cannot come from D
            if(i == nrow(shapeMat) || !shapeMat[i+1,j]){
                fromD[i,j] = 0	
            } else {
                inxD[i,j] = to1D(i+1,j)
            }			
            
            # nucleolus diffusion from D
            if(i != nrow(shapeMat) && nclMask[i,j] && nclMask[i+1,j])
                fromD[i,j] = dy.ncl
            
        }
    }
    
    # zero flux parameters
    p$leaveX = fromL + fromR
    p$leaveY = fromU + fromD
    
    # put into parameters
    p$fromL = fromL
    p$fromR = fromR
    p$fromU = fromU
    p$fromD = fromD
    p$inxL = inxL
    p$inxR = inxR
    p$inxU = inxU
    p$inxD = inxD
 
    # initialisation is to steady states calculated from the pre-bleachig images
    ss = steadyState(sum(p$pre.image.sim), p$shapeMat, p$k.on, p$k.off)
    F.init = ss$F
    C.init = ss$C
    
    p$F.init = F.init
    p$C.init = C.init
       
    F.init.vect = as.vector(F.init)
    C.init.vect = as.vector(C.init)
    
    p$FC.init.vect = c(F.init.vect, C.init.vect)
    
    # time
    t = c(0, p$bleach.duration)
    p$t = t
    
    p$bleach.profile = p$post.image.sim / p$pre.image.sim
    p$bleach.profile[is.nan(p$bleach.profile)] = 0
    
    # try all the values for bleach.norm.sd.range and record the best results
    cat("Estimating bleach.norm.sd range with ", length(p$bleach.norm.sd.range), "values\n")
    opt.all = list()
    for(norm.sd.inx in 1:length(p$bleach.norm.sd.range)){
        p.norm.sd = p
        p.norm.sd$bleach.norm.sd = p$bleach.norm.sd.range[norm.sd.inx]
        p.norm.sd$bleachRateless = createGaussianBleachingRatelessProfile(p.norm.sd)
    
        p.orig = p.norm.sd
        sims = list()
        
        opt = optimize(function(bleachRate){
            s = .simulateCircularBleaching2D_inner(p.orig, bleachRate)
            #(p.orig$bleach.depth - s$bleach.depth)^2 + (p.orig$bleach.fl.prop - s$bleach.fl.prop)^2
            sum((p$bleach.profile - s$bleach.profile)^2)
        }, bleachRates)
            
        opt.all[[norm.sd.inx]] = opt
    }
    
    # see which parameter combination has the smallest value of the error
    opt.obj = sapply(opt.all, function(x) x$objective)
    opt.inx = which.min(opt.obj)
    opt = opt.all[[opt.inx]]
    # initialise bleach.norm.sd to this value
    p.orig = p
    p.orig$bleach.norm.sd = p$bleach.norm.sd.range[opt.inx]
    p.orig$bleachRateless = createGaussianBleachingRatelessProfile(p.orig)
    
    # use the best parameter and return the values
    .simulateCircularBleaching2D_inner(p.orig, opt$minimum)
}


#' Find a set of parameters by fitting the simulated recovery curve
#'
#' @param sim simulation data
#' @param D.init initial value for diffusion constant D in the fitting search
#' @param k.on.init initial value for k.on
#' @param k.off.init initial value for k.off
fitMeans2D = function(sim, D.init=5*1e-6^2, k.on.init=0.01, k.off.init=0.01){
	p = sim$params
	
	shapeMat = p$shapeMat

	hx = p$L / p$grid.x
	hy = p$L / p$grid.y

	# Matrices of movement INTO tile (i,j) - this is the oppositve of what we have
	# in our Gillepsie simulations (which is out of)
	fromL = fromR = matrix(1, nrow=p$grid.y, ncol=p$grid.x)
	fromU = fromD = matrix(1, nrow=p$grid.y, ncol=p$grid.x)

	# leaveX and leaveY are how many ways a state can be left in the x/y directions
	leaveX = matrix(2, nrow=p$grid.y, ncol=p$grid.x)
	leaveY = matrix(2, nrow=p$grid.y, ncol=p$grid.x)

	# indicies of pixels that are left and right
	# the default is one out of range, which will be pre-filled with 0
	inxL = inxR = inxU = inxD = matrix(p$grid.x * p$grid.y + 1, nrow=p$grid.y, ncol=p$grid.x)

	#' Convert 2D to 1D coordinates
	#'
	#' @param i row coordinate
	#' @param j col coordinate
	#' @param nr number of rows
	to1D = function(i, j, nr=nrow(shapeMat)){
		i + (j-1)*nr
	}

	# create the diffusion matrices
	for(i in 1:nrow(shapeMat)){
		for(j in 1:ncol(shapeMat)){
			if(!shapeMat[i,j]){
				fromL[i,j] = fromR[i,j] = fromD[i,j] = fromU[i,j] = 0
				leaveX[i,j] = leaveY[i,j] = 0
				next
			}
			
			if(j == 1 || !shapeMat[i,j-1]){
				fromL[i,j] = 0
				# cannot leave to L
				leaveX[i,j] = leaveX[i,j] - 1
			} else {
				inxL[i,j] = to1D(i,j-1)
			}
			
			if(j == ncol(shapeMat) || !shapeMat[i,j+1]){
				fromR[i,j] = 0	
				# cannot leave to R
				leaveX[i,j] = leaveX[i,j] - 1
			} else {
				inxR[i,j] = to1D(i,j+1)
			}
			
			if(i == 1 || !shapeMat[i-1,j]){
				fromU[i,j] = 0
				# cannot leave U
				leaveY[i,j] = leaveY[i,j] - 1
			} else {
				inxU[i,j] = to1D(i-1,j)
			}
			
			if(i == nrow(shapeMat) || !shapeMat[i+1,j]){
				fromD[i,j] = 0	
				# cannot leave D
				leaveY[i,j] = leaveY[i,j] - 1
			} else {
				inxD[i,j] = to1D(i+1,j)
			}
		}
	}

	# initial state
	if(p$k.off == 0 | p$k.on == 0){
		F.eq = 1
		C.eq = 0
	} else {
		F.eq = p$k.off / (p$k.on + p$k.off)
		C.eq = p$k.on / (p$k.on + p$k.off)
	}
	
	F.init = (shapeMat + 0) * F.eq
	F.init[p$bleachMask] = 0

	C.init = (shapeMat + 0) * C.eq
	C.init[p$bleachMask] = 0

	F.init.vect = as.vector(F.init)
	C.init.vect = as.vector(C.init)

	init.vect = c(F.init.vect, C.init.vect)

	# time
	t = sim$t[sim$t > p$bleach.time] - p$bleach.time
	
	# bleach depth
	theta = sim$curve[min(which(sim$t > p$bleach.time))]	

	#' Function to minimise by optim()
	#' @param par parameters being optimised 
	frapLS = function(par){
		#cat("Called with", par, "\n")
	
		D = par[1]
		k.on = par[2]
		k.off = par[3]
		
		# fill in the rest of the parameters
		dx = D / hx^2
		dy = D / hy^2
		
		## full parameters
		pp = p
		pp$fromL = fromL * dx
		pp$fromR = fromR * dx
		pp$fromU = fromU * dy
		pp$fromD = fromD * dy
		pp$inxL = inxL
		pp$inxR = inxR
		pp$inxU = inxU
		pp$inxD = inxD
		pp$leaveX = leaveX * dx
		pp$leaveY = leaveY * dy
		pp$D = D
		pp$k.on = k.on
		pp$k.off = k.off
		
		# generate solution
		out = ode(init.vect, t, cppEllenberg2D_full_int, pp)
		
		# extract fluorescence
		s = list()
		s$FL = array(0, dim=c(p$grid.y, p$grid.x, length(t)))
		for(i in 1:nrow(out)){
			vect = out[i,-1]
			m = matrix(vect, ncol=2)
			s$FL[,,i] = matrix(m[,1], ncol=p$grid.x) + matrix(m[,2], ncol=p$grid.x)
		}
		s$params = pp
		
		# get the recovery curve
		s$curve = recoveryCurve(s)
		
		# scale and compare
		sum((sim$curve[sim$t > p$bleach.time] - (s$curve*(1-theta)+theta))^2)
	}
	
	best.param = optim(c(D.init, k.on.init, k.off.init), frapLS, method="L-BFGS-B",
		lower=c(1*1e-6^2, 0, 0), upper=c(100*1e-6^2, 1, 1), control=list(trace=1, REPORT=1))
}

#' Simulate the Ellenberg 3D model using a set of parameters
#'
#' @param params a set of parameters
#' @return sim object
simulateMeans3D = function(params){
	p = params
	.checkParams(p)
	
	shapeMat = p$shapeMat
	if("D.ncl" %in% names(p)){
		D.ncl = p$D.ncl
	} else {
		D.ncl = p$D
	}
	if("nclMask" %in% names(p)){
		nclMask = p$nclMask
	} else {
		nclMask = shapeMat
		nclMask[,] = FALSE
	}
	
	if(!("L.x" %in% names(p))){
		p$L.x = p$L		
	}
	if(!("L.y" %in% names(p))){
		p$L.y = p$L		
	}
	
	hx = p$L.x / p$grid.x
	hy = p$L.y / p$grid.y

	# diffusion: 10um^2/s
	dx = p$D / hx^2
	dy = p$D / hy^2
	
	dx.ncl = p$D.ncl / hx^2
	dy.ncl = p$D.ncl / hy^2

	# Matrices of movement INTO tile (i,j) - this is the oppositve of what we have
	# in our Gillepsie simulations (which is out of)
	fromL = fromR = matrix(dx, nrow=p$grid.y, ncol=p$grid.x)
	fromU = fromD = matrix(dy, nrow=p$grid.y, ncol=p$grid.x)

	# indicies of pixels that are left and right
	# the default is one out of range, which will be pre-filled with 0
	inxL = inxR = inxU = inxD = matrix(p$grid.x * p$grid.y + 1, nrow=p$grid.y, ncol=p$grid.x)

	#' Convert 2D to 1D coordinates
	#'
	#' @param i row coordinate
	#' @param j col coordinate
	#' @param nr number of rows
	to1D = function(i, j, nr=nrow(shapeMat)){
		i + (j-1)*nr
	}

	# create the diffusion matrices
	for(i in 1:nrow(shapeMat)){
		for(j in 1:ncol(shapeMat)){
			if(!shapeMat[i,j]){
				fromL[i,j] = fromR[i,j] = fromD[i,j] = fromU[i,j] = 0
				next
			}
			
			###########################################
			# cannot come from L
			if(j == 1 || !shapeMat[i,j-1]){
				fromL[i,j] = 0
			} else {
				inxL[i,j] = to1D(i,j-1)
			}
			
			# nucleolus diffusion from L
			if(j != 1 && nclMask[i,j] && nclMask[i,j-1])
				fromL[i,j] = dx.ncl
			
			###########################################
			# cannot come from R
			if(j == ncol(shapeMat) || !shapeMat[i,j+1]){
				fromR[i,j] = 0	
			} else {
				inxR[i,j] = to1D(i,j+1)
			}
			
			# nucleolus diffusion from R
			if(j != ncol(shapeMat) && nclMask[i,j] && nclMask[i,j+1])
				fromR[i,j] = dx.ncl
			
			###########################################
			# cannot come from U
			if(i == 1 || !shapeMat[i-1,j]){
				fromU[i,j] = 0
			} else {
				inxU[i,j] = to1D(i-1,j)
			}
			
			# nucleolus diffusion from U
			if(i != 1 && nclMask[i,j] && nclMask[i-1,j])
				fromU[i,j] = dy.ncl
			
			############################################
			# cannot come from D
			if(i == nrow(shapeMat) || !shapeMat[i+1,j]){
				fromD[i,j] = 0	
			} else {
				inxD[i,j] = to1D(i+1,j)
			}			
			
			# nucleolus diffusion from D
			if(i != nrow(shapeMat) && nclMask[i,j] && nclMask[i+1,j])
				fromD[i,j] = dy.ncl

		}
	}
	
	# zero flux parameters
	p$leaveX = fromL + fromR
	p$leaveY = fromU + fromD

	# put into parameters
	p$fromL = fromL
	p$fromR = fromR
	p$fromU = fromU
	p$fromD = fromD
	p$inxL = inxL
	p$inxR = inxR
	p$inxU = inxU
	p$inxD = inxD
	
	if(all(c("F.init", "C.init") %in% names(p))){
		# if specified use them
		F.init = p$F.init
		C.init = p$C.init
	} else {
		# initial state
		F.eq = p$k.off / (p$k.on + p$k.off)
		C.eq = p$k.on / (p$k.on + p$k.off)

		# this happens when k.on or k.off are 0... 	
		if(any(!is.finite(F.eq))){
			if(any(!is.finite(F.eq) & p$shapeMat))
				F.eq[!is.finite(F.eq) & p$shapeMat] = 1
			if(any(!is.finite(F.eq) & !p$shapeMat))
				F.eq[!is.finite(F.eq) & !p$shapeMat] = 0
		}
		if(any(!is.finite(C.eq))){
			C.eq[!is.finite(C.eq)] = 0
		}
	
		F.init = (shapeMat + 0) * F.eq
		F.init[p$bleachMask] = F.init[p$bleachMask] * p$bleach.depth

		C.init = (shapeMat + 0) * C.eq
		C.init[p$bleachMask] = C.init[p$bleachMask] * p$bleach.depth
	}

	F.init.vect = as.vector(F.init)
	C.init.vect = as.vector(C.init)

	init.vect = c(F.init.vect, C.init.vect)

	# time
	if("t" %in% names(p)){
		t = p$t[-length(p$t)]
	} else {
		t = seq(0, p$max.time-1/p$fps-p$bleach.time, 1/p$fps)
	}
	#out = ode(init.vect, t, cppEllenberg2D_full_int, p)
	out = ode(init.vect, t, .ellenberg2D_full_int, p)

	# need to convert 'out' into a sim object
	sim = list()
	sim$FL = array(0, dim=c(p$grid.y, p$grid.x, length(t)))
	sim$F = sim$C = array(0, dim=c(p$grid.y, p$grid.x, length(t)))
	for(i in 1:nrow(out)){
		vect = out[i,-1]
		m = matrix(vect, ncol=2)
		sim$F[,,i] = matrix(m[,1], ncol=p$grid.x)
		sim$C[,,i] = matrix(m[,2], ncol=p$grid.x)
		sim$FL[,,i] = sim$F[,,i] + sim$C[,,i]
	}
	sim$t = t
	sim$params = p
	if("rois" %in% names(p)){
		curve = matrix(0, nrow=dim(sim$F)[3], ncol=length(p$rois))
		colnames(curve) = names(p$rois)
		for(i in 1:length(p$rois)){
			curve[,i] = recoveryCurve(sim, p$rois[[i]])
		}
		sim$curve = curve
	} else {
		sim$curve = recoveryCurve(sim)
	}
	
	sim
}

#' Simulate the Ellenberg model using a set of parameters
#'
#' This version uses two population of bound molecules (C2 implementation)
#'
#' @param params a set of parameters
#' @return sim object
simulateMeans2D_C2 = function(params){
    p = params
    .checkParams(p)
    
    shapeMat = p$shapeMat
    if("D.ncl" %in% names(p)){
        D.ncl = p$D.ncl
    } else {
        D.ncl = p$D
    }
    if("nclMask" %in% names(p)){
        nclMask = p$nclMask
    } else {
        nclMask = shapeMat
        nclMask[,] = FALSE
    }
    
    if(!("L.x" %in% names(p))){
        p$L.x = p$L		
    }
    if(!("L.y" %in% names(p))){
        p$L.y = p$L		
    }
    
    hx = p$L.x / p$grid.x
    hy = p$L.y / p$grid.y
    
    # diffusion: 10um^2/s
    dx = p$D / hx^2
    dy = p$D / hy^2
    
    dx.ncl = p$D.ncl / hx^2
    dy.ncl = p$D.ncl / hy^2
    
    # Matrices of movement INTO tile (i,j) - this is the oppositve of what we have
    # in our Gillepsie simulations (which is out of)
    fromL = fromR = matrix(dx, nrow=p$grid.y, ncol=p$grid.x)
    fromU = fromD = matrix(dy, nrow=p$grid.y, ncol=p$grid.x)
    
    # indicies of pixels that are left and right
    # the default is one out of range, which will be pre-filled with 0
    inxL = inxR = inxU = inxD = matrix(p$grid.x * p$grid.y + 1, nrow=p$grid.y, ncol=p$grid.x)
    
    #' Convert 2D to 1D coordinates
    #'
    #' @param i row coordinate
    #' @param j col coordinate
    #' @param nr number of rows
    to1D = function(i, j, nr=nrow(shapeMat)){
        i + (j-1)*nr
    }
    
    # create the diffusion matrices
    for(i in 1:nrow(shapeMat)){
        for(j in 1:ncol(shapeMat)){
            if(!shapeMat[i,j]){
                fromL[i,j] = fromR[i,j] = fromD[i,j] = fromU[i,j] = 0
                next
            }
            
            ###########################################
            # cannot come from L
            if(j == 1 || !shapeMat[i,j-1]){
                fromL[i,j] = 0
            } else {
                inxL[i,j] = to1D(i,j-1)
            }
            
            # nucleolus diffusion from L
            if(j != 1 && nclMask[i,j] && nclMask[i,j-1])
                fromL[i,j] = dx.ncl
            
            ###########################################
            # cannot come from R
            if(j == ncol(shapeMat) || !shapeMat[i,j+1]){
                fromR[i,j] = 0	
            } else {
                inxR[i,j] = to1D(i,j+1)
            }
            
            # nucleolus diffusion from R
            if(j != ncol(shapeMat) && nclMask[i,j] && nclMask[i,j+1])
                fromR[i,j] = dx.ncl
            
            ###########################################
            # cannot come from U
            if(i == 1 || !shapeMat[i-1,j]){
                fromU[i,j] = 0
            } else {
                inxU[i,j] = to1D(i-1,j)
            }
            
            # nucleolus diffusion from U
            if(i != 1 && nclMask[i,j] && nclMask[i-1,j])
                fromU[i,j] = dy.ncl
            
            ############################################
            # cannot come from D
            if(i == nrow(shapeMat) || !shapeMat[i+1,j]){
                fromD[i,j] = 0	
            } else {
                inxD[i,j] = to1D(i+1,j)
            }			
            
            # nucleolus diffusion from D
            if(i != nrow(shapeMat) && nclMask[i,j] && nclMask[i+1,j])
                fromD[i,j] = dy.ncl
            
        }
    }
    
    # zero flux parameters
    p$leaveX = fromL + fromR
    p$leaveY = fromU + fromD
    
    # put into parameters
    p$fromL = fromL
    p$fromR = fromR
    p$fromU = fromU
    p$fromD = fromD
    p$inxL = inxL
    p$inxR = inxR
    p$inxU = inxU
    p$inxD = inxD
    
    if(all(c("F.init", "C1.init", "C2.init") %in% names(p))){
        # if specified use them
        F.init = p$F.init
        C1.init = p$C1.init
        C2.init = p$C2.init
    } else {
        cat("Generating F.init, C1.init, C2.init values assuming homogeneity adding up to 1\n")
        # initial state, from F + C1 + C2 = 1
        F.eq = 1 / (1 + p$k.on.c1/p$k.off.c1 + p$k.on.c2/p$k.off.c2)
        C1.eq = 1 - F.eq - p$k.on.c2/p$k.off.c2*F.eq
        C2.eq = 1 - F.eq - p$k.on.c1/p$k.off.c1*F.eq
        
        # this happens when k.on or k.off are 0... 	
        if(any(!is.finite(F.eq))){
            if(any(!is.finite(F.eq) & p$shapeMat))
                F.eq[!is.finite(F.eq) & p$shapeMat] = 1
            if(any(!is.finite(F.eq) & !p$shapeMat))
                F.eq[!is.finite(F.eq) & !p$shapeMat] = 0
        }
        if(any(!is.finite(C1.eq))){
            C1.eq[!is.finite(C1.eq)] = 0
        }
        if(any(!is.finite(C2.eq))){
            C2.eq[!is.finite(C2.eq)] = 0
        }
        
        F.init = (shapeMat + 0) * F.eq
        F.init[p$bleachMask] = F.init[p$bleachMask] * p$bleach.depth
        
        C1.init = (shapeMat + 0) * C1.eq
        C1.init[p$bleachMask] = C1.init[p$bleachMask] * p$bleach.depth
        
        C2.init = (shapeMat + 0) * C2.eq
        C2.init[p$bleachMask] = C2.init[p$bleachMask] * p$bleach.depth
    }
    
    F.init.vect = as.vector(F.init)
    C1.init.vect = as.vector(C1.init)
    C2.init.vect = as.vector(C2.init)
    
    init.vect = c(F.init.vect, C1.init.vect, C2.init.vect)
    
    # time
    if("t" %in% names(p)){
        t = p$t[-length(p$t)]
    } else {
        t = seq(0, p$max.time-1/p$fps-p$bleach.time, 1/p$fps)
    }
    # Numerical integration
    out = ode(init.vect, t, .ellenberg2D_C2_full_int, p)
    
    # need to convert 'out' into a sim object
    sim = list()
    sim$F.init = F.init
    sim$C1.init = C1.init
    sim$C2.init = C2.init
    sim$FL = array(0, dim=c(p$grid.y, p$grid.x, length(t)))
    sim$F = sim$C1 = sim$C2 = array(0, dim=c(p$grid.y, p$grid.x, length(t)))
    for(i in 1:nrow(out)){
        vect = out[i,-1]
        m = matrix(vect, ncol=3)
        sim$F[,,i] = matrix(m[,1], ncol=p$grid.x)
        sim$C1[,,i] = matrix(m[,2], ncol=p$grid.x)
        sim$C2[,,i] = matrix(m[,3], ncol=p$grid.x)
        sim$FL[,,i] = sim$F[,,i] + sim$C1[,,i] + sim$C2[,,i]
    }
    sim$t = t
    sim$params = p
    if("rois" %in% names(p)){
        curve = matrix(0, nrow=dim(sim$F)[3], ncol=length(p$rois))
        colnames(curve) = names(p$rois)
        for(i in 1:length(p$rois)){
            curve[,i] = recoveryCurve(sim, p$rois[[i]])
        }
        sim$curve = curve
    } else {
        sim$curve = recoveryCurve(sim)
    }
    
    sim
}
