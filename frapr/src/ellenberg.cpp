#include <Rcpp.h>

#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>

using namespace Rcpp;

RcppExport SEXP ellenberg2D_full_intCpp(SEXP inT, SEXP inState, SEXP inP) {
	// initialise the input parameters	
	Rcpp::NumericVector t(inT);
	Rcpp::NumericVector state(inState);
	Rcpp::List p(inP);
	
	// this will be the output
	Rcpp::NumericVector out(state.size());
	
	// extract grid parameters
	Rcpp::NumericVector gridX = p["grid.x"];
	Rcpp::NumericVector gridY = p["grid.y"];

	// size of both F and C	
	int offset = gridX[0] * gridY[0];
	
	// copy F and C otherwise the code is too messy
	Rcpp::NumericVector F(offset+1);
	Rcpp::NumericVector C(offset);
	
	for(int i=0;i<offset;i++){
		F[i] = state[i];
		C[i] = state[offset+i];
	}
	// missing value
	F[offset] = 0;
	
	// get out the variables we need from p
	Rcpp::NumericMatrix fromL = p["fromL"];
	Rcpp::NumericMatrix fromR = p["fromR"];
	Rcpp::NumericMatrix fromU = p["fromU"];
	Rcpp::NumericMatrix fromD = p["fromD"];
	Rcpp::NumericMatrix inxL = p["inxL"];
	Rcpp::NumericMatrix inxR = p["inxR"];
	Rcpp::NumericMatrix inxU = p["inxU"];
	Rcpp::NumericMatrix inxD = p["inxD"];
	Rcpp::NumericMatrix leaveX = p["leaveX"];
	Rcpp::NumericMatrix leaveY = p["leaveY"];
	Rcpp::NumericVector kOn = p["k.on"];
	Rcpp::NumericVector kOff = p["k.off"];

	double k_on = kOn[0];
	double k_off = kOff[0];
	
	// the main loop
	for(int i=0;i<offset;i++){
		out[i] = fromL[i] * F[inxL[i]-1] + fromR[i] * F[inxR[i]-1] +
		         fromD[i] * F[inxD[i]-1] + fromU[i] * F[inxU[i]-1] -
		         leaveX[i] * F[i] - leaveY[i] * F[i] -
		         k_on * F[i] + k_off * C[i];
		         
		out[offset+i] = k_on * F[i] - k_off * C[i];
	}
	
	
	return(List::create(out));		
}

/** R code:
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
**/



