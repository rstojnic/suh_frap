# Curve fitting methods from McNally 2008 paper

#' Predict FRAP curve based on diffusion (S.22)
#' 
#' Note this function assumes that during the time given
#' the curve reaches its steady state!!!
#'
#' @param t time vector
#' @param Df diffusion constant
#' @param alpha alpha parameters (chi / Rn)
#' @param U U parameters (S.21)
#' @param J0 J0 parameters (S.17)
.mcNally2008_frapDiffusion = function(t, Df, alpha, U, J0, num.terms){
	frap = rep(0, length(t))
	alpha = c(0, alpha)
	for(i in 0:(num.terms-1)){
		frap = frap + U[i+1] * exp(-Df*alpha[i+1]^2*t) * J0[i+1]
	}      
	frap = frap / max(frap)

	frap
}

#' Fit a diffusion only model from the McNally 2008 paper
#'
#' @param sim Simulation object with calculated recovery curve
#' @param Df.all Diffusion constant ranges to consider in fitting, default from 1 to 100 um^2/s
#' @param num.terms Number of terms to use in the expansion
#' @param list of:
#'          t - time
#'          fitted - predicted recovery
#'          D - inferred diffusion constant
fitMcNally2008_diffusionOnly = function(sim, Df.all=seq(1 * (1e-6)^2, 100 * (1e-6)^2, 0.25 * (1e-6)^2), num.terms=500){
	p = sim$params
	bleach.time = p$bleach.time
	
	t = sim$t[sim$t>bleach.time]
	recovery = sim$curve[sim$t>bleach.time]
	
	# normalise to start from 0
	t = t - min(t)
	
	# calculate radius parameters
	if(!("L.x" %in% names(p)))
		p$L.x = p$L
	if(!("L.y" %in% names(p)))
		p$L.y = p$L		
	pixel.x = p$L.x / p$grid.x
	pixel.y = p$L.y / p$grid.y
	
	# parameters
	Rn = sqrt(sum(p$shapeMat)*pixel.x*pixel.y / pi)
	Rm = sqrt(sum(p$bleachMask)*pixel.x*pixel.y / pi)
	Rc = Rm

	# do the calculation
	chi = bessel_zero_J1(1:(num.terms-1))
	alpha = chi / Rn

	# S.17
	J0 = c(1, 2 / (alpha * Rm) * bessel_J1(alpha*Rm))

	# bleach depth
	theta = p$bleach.depth

	# S.21
	U = c(1 + (theta - 1) * Rc^2 / Rn^2,
		  (theta - 1) * 2 * Rc / (alpha * Rn^2) * bessel_J1(alpha * Rc) / (bessel_J0(chi)^2))
	
	#' Least square different to the data
	#'
	#' @param Df diffusion rate
	frapLS = function(Df){
		f = .mcNally2008_frapDiffusion(t, Df, alpha, U, J0, num.terms)
		sum((f - recovery)^2)
	}

	# Try diffusion rates between 1 and 100, and select the best one
	ls.all = sapply(Df.all, frapLS)

	Df.out = Df.all[which.min(ls.all)]
	
	# Calculate the predicted recovery curve for the predicted Df and plot
	frap = .mcNally2008_frapDiffusion(t, Df.out, alpha, U, J0, num.terms)
	
	list(t=sim$t[sim$t>bleach.time], fitted=frap, D=Df.out)
	
}

#' Predict the FRAP curve based on the full model (S.16)
#'
#' Note this function assumes that during the time given
#' the curve reaches its steady state!!!
#'
#' @param t time vector (starting at 0)
#' @param Df diffusion constant
#' @param k_on pseudo-on rate
#' @param k_off off rate
#' @param chi the zeros of J1
#' @param alpha alpha parameters (chi/Rn)
#' @param J0 J0 parameters (S.17)
#' @param theta bleach deapth
#' @param Rc radius of the bleach spot
#' @param Rn radius of the nucleus
#' @param num.terms number of terms to consider in the expansion
.mcNally2008_frapFull = function(t, Df, k_on, k_off, chi, alpha, J0, theta, Rc, Rn, num.terms){
	# equilibirum concentrations
	F_eq = k_off / (k_on + k_off)
	C_eq = k_on / (k_on + k_off)

	# S.4
	w = 0.5 * (Df*alpha^2 + k_on + k_off)
	v = sqrt(0.25 * (Df*alpha^2 + k_on + k_off)^2 - k_off * Df * alpha^2)

	# S.5 
	beta = w + v
	gamma = w - v

	# S.15
	Z = c(F_eq * (1 + (theta-1) * Rc^2 / Rn^2),  # k = 0
		  (theta-1) * F_eq * 2*Rc/(alpha[-1]*Rn^2) * bessel_J1(alpha[-1]*Rc) / bessel_J0(chi[-1])^2) # k != 0

	# S.11
	U = 1 / (-2*k_off*v) * (-w - v + k_off) * (w - v) * Z
	V = 1 / (2*k_off*v) * (-w + v + k_off) * (w + v) * Z

	# S.10
	W = U * k_on / (-beta + k_off)
	X = V * k_on / (-gamma + k_off)

	# now FRAP curve!
	frap = rep(0, length(t))
	for(k in 1:num.terms){ # note this is really 0..(num.terms-1) because of R indexing from 1
		frap = frap + ((U[k] + W[k]) * exp(-beta[k]*t) + (V[k] + X[k]) * exp(-gamma[k]*t)) * J0[k]
	}      
	frap = frap / max(frap)

	frap
}


#' Fit a full model from the McNally 2008 paper
#'
#' @param sim Simulation object with calculated recovery curve
#' @param Df.all Diffusion constant ranges to consider in fitting, default from 1 to 100 um^2/s
#' @param k.off.all the k.off values
#' @param k.on.all either provide the k.on grid, or the Free grid
#' @param Free.all Free molecules grid
#' @param num.terms Number of terms to use in the expansion
#' @param list of:
#'          t - time
#'          fitted - predicted recovery
#'          D - inferred diffusion constant
fitMcNally2008_full = function(sim, Df.all=seq(1 * (1e-6)^2, 50 * (1e-6)^2, 1 * (1e-6)^2), 
	k.on.all=NULL, Free.all=NULL, k.off.all=seq(1/30, 0.5, 1/30), num.terms=500){
	
	p = sim$params
	bleach.time = p$bleach.time

	t = sim$t[sim$t>bleach.time]
	recovery = sim$curve[sim$t>bleach.time]

	# normalise to start from 0
	t = t - min(t)

	# calculate radius parameters
	if( !("L.x" %in% names(p)) ) p$L.x = p$L
	if( !("L.y" %in% names(p)) ) p$L.y = p$L
	pixel.x = p$L.x / p$grid.x
	pixel.y = p$L.y / p$grid.y

	# parameters
	Rn = sqrt(sum(p$shapeMat)*pixel.x*pixel.y / pi)
	Rm = sqrt(sum(p$bleachMask)*pixel.x*pixel.y / pi)
	Rc = Rm

	# do the calculation
	chi = bessel_zero_J1(0:(num.terms-1))
	alpha = chi / Rn

	# S.17
	J0 = c(1, 2 / (alpha[-1] * Rm) * bessel_J1(alpha[-1]*Rm))

	# bleach depth
	theta = p$bleach.depth
	
	#' Least square different to the data
	#'
	#' @param Df diffusion rate
	frapLS = function(Df, k.on, k.off){
		f = .mcNally2008_frapFull(t, Df, k.on, k.off, chi, alpha, J0, theta, Rc, Rn, num.terms)
		sum((f - recovery)^2)
	}
	
	# fit the values with LS
	if(!is.null(k.on.all)){
	    ls.all = array(NA, dim=c(length(Df.all), length(k.on.all), length(k.off.all)))
	} else {
	    ls.all = array(NA, dim=c(length(Df.all), length(Free.all), length(k.off.all)))   
	}
	for(i in 1:length(Df.all)){
		cat("Doing", i, "/", length(Df.all), "\n")
	    for(k in 1:length(k.off.all)){
	        if(!is.null(k.on.all)){
	            for(j in 1:length(k.on.all)){
	                ls.all[i,j,k] = frapLS(Df.all[i], k.on.all[j], k.off.all[k])
	            }
	        } else {
	            for(j in 1:length(Free.all)){
	                k.on.j = ((1-Free.all[j])/Free.all[j]) * k.off.all[k]
	                ls.all[i,j,k] = frapLS(Df.all[i], k.on.j, k.off.all[k])
	            }
	        }
		}
	}

	# best parameter combination
	inx = which(ls.all == min(ls.all), arr.ind=TRUE)
	Df.out = Df.all[inx[1,1]]
	k.off.out = k.off.all[inx[1,3]]
	if(!is.null(k.on.all)){
	    k.on.out = k.on.all[inx[1,2]]
	    Free.out = k.off.out / (k.on.out + k.off.out)
	} else {
	    Free.out = Free.all[inx[1,2]]
	    k.on.out = ((1-Free.out)/Free.out) * k.off.out
	}
	
	# Calculate the predicted recovery curve for the predicted Df and plot
	frap = .mcNally2008_frapFull(t, Df.out, k.on.out, k.off.out, chi, alpha, J0, theta, Rc, Rn, num.terms)
	
	list(t=sim$t[sim$t>bleach.time], fitted=frap, D=Df.out, k.on=k.on.out, Free=Free.out, k.off=k.off.out, res.time=1/k.off.out)
}
