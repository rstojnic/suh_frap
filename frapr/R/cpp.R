# Interface for Cpp functions

#' Convert a list of matrices into a 3D array
#'
#' @param x a list of matrices
.listToArray = function(x){
	# lists might be actually shorter due to how Rcpp handles lists
	real.len = max(which(!sapply(x, is.null)))
	# copy
	out = array(0, dim=c(dim(x[[1]]), real.len))
	for(i in 1:real.len){
		out[,,i] = x[[i]]
	}
	out
}

#' Convert a list of matrices into a 4D array
#'
#' @param x a list of 3D arrays
.listToArray3D = function(x){
	# lists might be actually shorter due to how Rcpp handles lists
	real.len = max(which(!sapply(x, is.null)))
	# copy
	out = array(0, dim=c(dim(x[[1]]), real.len))
	for(i in 1:real.len){
		out[,,,i] = x[[i]]
	}
	out
}

#' Gillespie simulation of a simple 1 binding site system
#'
#' @param params a list of parameters:
#'   grid.x number of horizontal tiles
#'   grid.y number of vertical tiles
#'   mol molecules per tile
#'   k.on the ON rate of rection F -> C
#'   k.off the OFF Rate of rection C -> F
#'   L the length of the grid (both vertical and horizontal)
#'   D diffusion constant
#'   shapeMat a logical matrix specifying where the molecules can be
#'   max.time end time of the simulation
#'   bleach.time time at which photobleaching should be done
#'   bleachMask a logical matrix specifying which tiles should be bleached
#'   bleach.depth the value to which bleaching should be done. Default is 0.5. A value of 0 will produce 100% bleaching. 
#'   fps frames per second to record 
gillespieSim = function(params){
	p = params

	.checkParams(p)
	
	F.total = sum(p$shapeMat) * p$mol
			
	# call an implementation in Rcpp
	res = .Call("gillespieSimCpp", p$grid.x, p$grid.y, F.total, p$k.on, p$k.off, p$L, p$D, p$shapeMat + 0, p$max.time, 
		p$bleach.time, p$bleachMask, p$bleach.depth, p$fps)
		
	# convert results into an array
	out = list()	
	out$F = .listToArray(lapply(res, function(x) x$F))
	out$FL_F = .listToArray(lapply(res, function(x) x$FL))
	out$C = .listToArray(lapply(res, function(x) x$C))
	out$FL_C = .listToArray(lapply(res, function(x) x$FLC))
	out$FL = out$FL_F + out$FL_C
	out$t = unlist(sapply(res, function(x) x$t))
	out$params = p
	
	out$curve = recoveryCurve(out)
		
	out
}
#' 3D version of gillespieSim, with the same parameters
#'
#' The only additional parameter is params$grid.z which specified number of depth tiles
#'
gillespieSim3D = function(params){
	p = params

	.checkParams3D(p)
	
	F.total = sum(p$shapeMat) * p$mol
			
	# call an implementation in Rcpp
	res = .Call("gillespieSim3DCpp", p$grid.x, p$grid.y, p$grid.z, F.total, p$k.on, p$k.off, p$L, p$D, p$shapeMat + 0, p$max.time, 
		p$bleach.time, p$bleachMask, p$bleach.depth, p$fps)
		
	# convert results into an array
	out = list()	
	out$F = .listToArray3D(lapply(res, function(x) x$F))
	out$FL_F = .listToArray3D(lapply(res, function(x) x$FL))
	out$C = .listToArray3D(lapply(res, function(x) x$C))
	out$FL_C = .listToArray3D(lapply(res, function(x) x$FLC))
	out$FL = out$FL_F + out$FL_C
	out$t = unlist(sapply(res, function(x) x$t))
	out$params = p
	#out$raw = res
	
	#out$curve = recoveryCurve(out)
		
	out
}

cppEllenberg2D_full_int = function(t, state, p){
	.Call("ellenberg2D_full_intCpp", t, state, p)
}
