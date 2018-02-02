# functions for FRAP

#' Solve mass action equilibrium for A+B <-> C
#'
#' @param A total A molecules
#' @param B total B molecules
#' @param k dissociation constant (koff/kon)
solve.ma = function(A, B, k){
	c1 = ((A+B+k) + sqrt((A+B+k)^2 - 4*A*B))/2
	c2 = ((A+B+k) - sqrt((A+B+k)^2 - 4*A*B))/2

	if(c1 <= 0 | (c1 > A | c1 > B)){
		c = c2
	} else {
		c = c1
	}	

	a = A - c
	b = B - c

	#cat("a=", a, "b=", b, "c=", c, "\n")
	list(A=A, B=B, k=k, a = a, b = b, ab = c)
}

#' Draw a filled circle inside a logical matrix
#'
#' @param ci the row of the center
#' @param cj the column of the center
#' @param radius the radius of the circle
drawCircle = function(mat, ci, cj, radius){
	for(i in 1:nrow(mat)){
		for(j in 1:ncol(mat)){
			mi = i - 0.5
			mj = j - 0.5
			
			# see if the midpoint is within the radius, of so, colour the pixel
			if((sqrt((mi-ci)^2 + (mj-cj)^2))<=radius){
				mat[i,j] = TRUE
			} else {
				mat[i,j] = FALSE
			}
		}
	}
	mat
}

#' Draw a filled ball inside a logical matrix
#'
#' @param ci the row of the center
#' @param cj the column of the center
#' @param ck the z of the center
#' @param radius the radius of the circle
drawCircle3D = function(mat, ci, cj, ck, radius){
	for(i in 1:dim(mat)[[1]]){
		for(j in 1:dim(mat)[[2]]){
			for(k in 1:dim(mat)[[3]]){
				mi = i - 0.5
				mj = j - 0.5
				mk = k - 0.5
			
				# see if the midpoint is within the radius, of so, colour the pixel
				if((sqrt((mi-ci)^2 + (mj-cj)^2 + (mk-ck)^2))<=radius){
					mat[i,j,k] = TRUE
				} else {
					mat[i,j,k] = FALSE
				}
			}
		}
	}
	mat
}

#' Check if a set of parameters is valid
#'
#' @param p a list of parameters
.checkParams = function(p){
	# check all the parameters !!!
	required.params = c("grid.x", "grid.y", "mol", "k.on", "k.off", "L", "D", "shapeMat", "max.time", "bleach.time",
		"bleachMask", "bleach.depth", "fps")
		
	params.diff = setdiff(required.params, names(p))
	params.extra = setdiff(names(p), required.params)
	
	if(length(params.diff)>0){
		stop(paste("Parameters:", paste(params.diff, collapse=","), "are missing"))
	}
	
	# check types and ranges
	if(!is.numeric(p$grid.x) || floor(p$grid.x)!=p$grid.x){
		stop("grid.x parameter needs to be an integer value")
	}
	if(!is.numeric(p$grid.y) || floor(p$grid.y)!=p$grid.y){
		stop("grid.y parameter needs to be an integer value")
	}
	# all of these should be numeric
	for(name in c("mol", "k.on", "k.off", "L", "D", "max.time", "bleach.time", "bleach.depth", "fps")){
		if(!is.numeric(p[[name]])){
			stop(paste(name, "needs to be a numeric value"))
		}
	}
	for(name in c("shapeMat", "bleachMask")){
		if(!is.matrix(p[[name]]) || !is.logical(p[[name]])){
			stop(paste(name, "needs to be a matrix of logical values"))
		}
		
		if(!all(dim(p[[name]]) == c(p$grid.y, p$grid.x))){
			stop(paste(name, "needs to be a matrix of size (grid.y, grid.x)"))
		}
	}
}

#' Check if a set of parameters is valid (3D simulations)
#'
#' @param p a list of parameters
.checkParams3D = function(p){
	# check all the parameters !!!
	required.params = c("grid.x", "grid.y", "grid.z", "mol", "k.on", "k.off", "L", "D", "shapeMat", "max.time", "bleach.time",
		"bleachMask", "bleach.depth", "fps")
		
	params.diff = setdiff(required.params, names(p))
	params.extra = setdiff(names(p), required.params)
	
	if(length(params.diff)>0){
		stop(paste("Parameters:", paste(params.diff, collapse=","), "are missing"))
	}
	
	# check types and ranges
	if(!is.numeric(p$grid.x) || floor(p$grid.x)!=p$grid.x){
		stop("grid.x parameter needs to be an integer value")
	}
	if(!is.numeric(p$grid.y) || floor(p$grid.y)!=p$grid.y){
		stop("grid.y parameter needs to be an integer value")
	}
	if(!is.numeric(p$grid.z) || floor(p$grid.z)!=p$grid.z){
		stop("grid.y parameter needs to be an integer value")
	}
	# all of these should be numeric
	for(name in c("mol", "k.on", "k.off", "L", "D", "max.time", "bleach.time", "bleach.depth", "fps")){
		if(!is.numeric(p[[name]])){
			stop(paste(name, "needs to be a numeric value"))
		}
	}
	for(name in c("shapeMat", "bleachMask")){
		if(!is.array(p[[name]]) || !is.logical(p[[name]])){
			stop(paste(name, "needs to be an array of logical values"))
		}
		
		if(!all(dim(p[[name]]) == c(p$grid.y, p$grid.x, p$grid.z))){
			stop(paste(name, "needs to be an array of size (grid.y, grid.x, grid.z)"))
		}
	}
}

#' Rescale the values of x so that the minimam value is the same, and the new value
#' is mapped to 1
#'
#' @param x vector or matrix of values
#' @param new.one value that needs to be mapped to 1
rescaleOne = function(x, new.one){
	min.val = min(na.omit(x))
	
	x.norm = x - min.val
	x.one = 1 - min.val
	one.val = new.one - min.val
	
	min.val + (x.norm * x.one/one.val)
}

#' Normalise into the 0-1 range preserving the 1
#'
#' @param x a vector
norm01 = function(x) {
	mx = min(na.omit(x))
	norm = x - mx
	norm.1 = 1 - mx
	
	norm / norm.1
}
