#include <RcppArmadillo.h>

#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>

#define DEBUG_3D false
#define CONSISTENCY_CHECK_3D false

/*
	This macro is used as a general template for diffusion operations
	
	@param PROP propensity matrix, e.g. LR
	@param i current i value
	@param j current j value
	@param k current k value
	@param ni i coordinate that is modified
	@param nj j coordinate that is modified
	@param nk k coordinate that is modified
*/
#define DIFFUSION3D(PROP, I, J, K, NI, NJ, NK) ({ \
	if(PROP(I,J,K) != 0 && action_sum < action_val && action_val <= action_sum_next){ \
		/* how many are fluorescent */ \
		fl_prop = FL(I,J,K) / F(I,J,K); \
		\
		/* calculate the difference in total propensities */ \
		alpha_diff = - LR(I,J,K) - RL(I,J,K) - UD(I,J,K) - DU(I,J,K) - AB(I,J,K) - BA(I,J,K) - r_on(I,J,K) - \
		               LR(NI,NJ,NK) - RL(NI,NJ,NK) - UD(NI,NJ,NK) - DU(NI,NJ,NK) - AB(NI,NJ,NK) - BA(NI,NJ,NK) - r_on(NI,NJ,NK); \
		\
		/* update */ \
		F(I,J,K) -= 1; \
		F(NI,NJ,NK) += 1; \
		/* make sure to update all propensities!!! */ \
		LR(I,J,K) = D_LR(I,J,K) * F(I,J,K); \
		RL(I,J,K) = D_RL(I,J,K) * F(I,J,K); \
		UD(I,J,K) = D_UD(I,J,K) * F(I,J,K); \
		DU(I,J,K) = D_DU(I,J,K) * F(I,J,K); \
		AB(I,J,K) = D_AB(I,J,K) * F(I,J,K); \
		BA(I,J,K) = D_BA(I,J,K) * F(I,J,K); \
		r_on(I,J,K) = k_on * F(I,J,K); \
		/* and for ni,nj */\
		LR(NI,NJ,NK) = D_LR(NI,NJ,NK) * F(NI,NJ,NK); \
		RL(NI,NJ,NK) = D_RL(NI,NJ,NK) * F(NI,NJ,NK); \
		UD(NI,NJ,NK) = D_UD(NI,NJ,NK) * F(NI,NJ,NK); \
		DU(NI,NJ,NK) = D_DU(NI,NJ,NK) * F(NI,NJ,NK); \
		AB(NI,NJ,NK) = D_AB(NI,NJ,NK) * F(NI,NJ,NK); \
		BA(NI,NJ,NK) = D_BA(NI,NJ,NK) * F(NI,NJ,NK); \
		r_on(NI,NJ,NK) = k_on * F(NI,NJ,NK);\
		\
		/* add new values to propensitiy */ \
		alpha_diff +=  LR(I,J,K) + RL(I,J,K) + UD(I,J,K) + DU(I,J,K) + AB(I,J,K) + BA(I,J,K) + r_on(I,J,K) + \
		               LR(NI,NJ,NK) + RL(NI,NJ,NK) + UD(NI,NJ,NK) + DU(NI,NJ,NK) + AB(NI,NJ,NK) + BA(NI,NJ,NK) + r_on(NI,NJ,NK); \
		\
		/* update the fluorescent ones as well */ \
		if(runif(1)[0] < fl_prop){ \
			FL(I,J,K) -= 1; \
			FL(NI,NJ,NK) += 1; \
		} \
		\
		/* update alpha0 */ \
		alpha0 += alpha_diff; \
		\
		finished = true; \
		break; \
	} \
})

/* 
	This macro is used to model a general reaction operation
	
	@param PROP the propensity matrix, e.g. r_on
	@param FLOR fluorescence species, e.g. FL
	@param FLOR_TOTAL total molecules corresponding to FLOR, e.g. F
	@param F_ADD value to add to F, e.g. -1
	@param C_ADD value to add to C, e.g. 1
	
*/
#define REACTION(PROP, FLOR, FLOR_TOTAL, F_ADD, C_ADD) ({ \
	if(PROP(i,j,k) != 0 && action_sum < action_val && action_val <= action_sum_next){ \
		/* how many are fluorescent */ \
		fl_prop = FLOR(i,j,k) / FLOR_TOTAL(i,j,k); \
		\
		/* calculate the difference in total propensities */ \
		alpha_diff = - r_on(i,j,k) - r_off(i,j,k) - LR(i,j,k) - RL(i,j,k) - UD(i,j,k) - DU(i,j,k) - AB(i,j,k) - BA(i,j,k); \
		\
		/* update */ \
		F(i,j,k) += F_ADD; \
		C(i,j,k) += C_ADD; \
		LR(i,j,k) = D_LR(i,j,k) * F(i,j,k); \
		RL(i,j,k) = D_RL(i,j,k) * F(i,j,k); \
		UD(i,j,k) = D_UD(i,j,k) * F(i,j,k); \
		DU(i,j,k) = D_DU(i,j,k) * F(i,j,k); \
		AB(i,j,k) = D_AB(i,j,k) * F(i,j,k); \
		BA(i,j,k) = D_BA(i,j,k) * F(i,j,k); \
		r_on(i,j,k) = F(i,j,k) * k_on; \
		r_off(i,j,k) = C(i,j,k) * k_off; \
		\
		/* add new values to propensitiy */ \
		alpha_diff += r_on(i,j,k) + r_off(i,j,k) + LR(i,j,k) + RL(i,j,k) + UD(i,j,k) + DU(i,j,k) + AB(i,j,k) + BA(i,j,k); \
		\
		/* update the fluorescent ones as well */ \
		if(runif(1)[0] < fl_prop){ \
			FL(i,j,k) += F_ADD; \
			FLC(i,j,k) += C_ADD; \
		} \
		\
		/* update alpha0 */ \
		alpha0 += alpha_diff; \
		\
		finished = true; \
		break; \
	} \
})

using namespace Rcpp;

arma::cube cloneCube(arma::cube cubeArray){
	arma::cube cloned(cubeArray);
	
	return cloned;
}

RcppExport SEXP gillespieSim3DCpp(SEXP inGridX, SEXP inGridY, SEXP inGridZ, SEXP inFtotal, SEXP inKon, SEXP inKoff, SEXP inL, SEXP inD, SEXP inShapeMat, SEXP inMaxTime, 
	SEXP inFrapTime, SEXP inFrapMask, SEXP inBleachDepth, SEXP inFps) {
	// initialise the input parameters
	Rcpp::IntegerVector rGridX(inGridX);
	Rcpp::IntegerVector rGridY(inGridY);
	Rcpp::IntegerVector rGridZ(inGridZ);
	Rcpp::IntegerVector rFtotal(inFtotal);
	Rcpp::NumericVector rkOn(inKon);
	Rcpp::NumericVector rkOff(inKoff);
	Rcpp::NumericVector rL(inL);
	Rcpp::NumericVector rD(inD);
	Rcpp::NumericVector rMaxTime(inMaxTime);
	Rcpp::NumericVector rFrapTime(inFrapTime);
	Rcpp::NumericVector rFps(inFps);
	Rcpp::NumericVector rBleachDepth(inBleachDepth);
	
	// convert frap mask to 3D matrix
	Rcpp::NumericVector frapMaskVec(inFrapMask);
	Rcpp::IntegerVector frapMaskArrayDims = frapMaskVec.attr("dim"); 
  	arma::cube frapMask(frapMaskVec.begin(), frapMaskArrayDims[0], frapMaskArrayDims[1], frapMaskArrayDims[2], false);
	
	// convert shapeMat to 3D matrix
	Rcpp::NumericVector shapeMatVec(inShapeMat);
	Rcpp::IntegerVector shapeMatArrayDims = shapeMatVec.attr("dim"); 
  	arma::cube shapeMat(shapeMatVec.begin(), shapeMatArrayDims[0], shapeMatArrayDims[1], shapeMatArrayDims[2], false);

	// convert into normal C++ data structures
	int grid_x = rGridX[0];
	int grid_y = rGridY[0];
	int grid_z = rGridZ[0];
	int F_total = rFtotal[0];
	double k_on = rkOn[0];
	double k_off = rkOff[0];
	double L = rL[0];
	double D = rD[0];
	double max_time = rMaxTime[0];
	double frap_time = rFrapTime[0];
	double fps = rFps[0];
	double bleach_depth = rBleachDepth[0];
	
	// size of each compartment	
	double hx = L / grid_x;
	double hy = L / grid_y;
	double hz = L / grid_z;

	// calculate diffusion coefficients
	double dx = D / (hx*hx);
	double dy = D / (hy*hy);
	double dz = D / (hz*hz);
	
	// create the diffusion matrices for the 6 directions
	arma::cube D_LR(grid_y, grid_x, grid_z);
	arma::cube D_RL(grid_y, grid_x, grid_z);
	arma::cube D_UD(grid_y, grid_x, grid_z);
	arma::cube D_DU(grid_y, grid_x, grid_z);
	arma::cube D_AB(grid_y, grid_x, grid_z);
	arma::cube D_BA(grid_y, grid_x, grid_z);

	
	for(int i=0;i<grid_y;i++){
		for(int j=0;j<grid_x;j++){
			for(int k=0;k<grid_z;k++){
				// LR
				if(j == (grid_x - 1) || shapeMat(i,j,k) == 0 || shapeMat(i,j+1,k) == 0){
					D_LR(i,j,k) = 0;
				} else {
					D_LR(i,j,k) = dx;
				}
			
				// RL
				if(j == 0 || shapeMat(i,j,k) == 0 || shapeMat(i,j-1,k) == 0){
					D_RL(i,j,k) = 0;
				} else {
					D_RL(i,j,k) = dx;
				}
			
				// UD
				if(i == (grid_y - 1) || shapeMat(i,j,k)==0 || shapeMat(i+1,j,k) == 0){
					D_UD(i,j,k) = 0;
				} else {
					D_UD(i,j,k) = dy;
				}
			
				// DU
				if(i == 0 || shapeMat(i,j,k)==0 || shapeMat(i-1,j,k) == 0){
					D_DU(i,j,k) = 0;
				} else {
					D_DU(i,j,k) = dy;
				}
				
				// AB - smaller indicies are more UP, i.e the shapeMate(i,j,0) is the top Z-stack
				if(k == (grid_z - 1) || shapeMat(i,j,k) == 0 || shapeMat(i,j,k+1) == 0){
					D_AB(i,j,k) = 0;
				} else {
					D_AB(i,j,k) = dz;
				}
			
				// BA
				if(k == 0 || shapeMat(i,j,k) == 0 || shapeMat(i,j,k-1) == 0){
					D_BA(i,j,k) = 0;
				} else {
					D_BA(i,j,k) = dz;
				}
				
			}
		}
	}
	
	//////////////////////////////////////////////////
	
	// init
	bool frapped = false;
	int event_inx = 0;
	double t = 0;
	
	double r1,r2;
	
	// initial distribution of molecules
	arma::cube F(grid_y, grid_x, grid_z);
	arma::cube FL(grid_y, grid_x, grid_z);
	arma::cube C(grid_y, grid_x, grid_z);
	arma::cube FLC(grid_y, grid_x, grid_z);
	
	// propensity matrices
	arma::cube LR(grid_y, grid_x, grid_z);
	arma::cube RL(grid_y, grid_x, grid_z);
	arma::cube UD(grid_y, grid_x, grid_z);
	arma::cube DU(grid_y, grid_x, grid_z);
	arma::cube AB(grid_y, grid_x, grid_z);
	arma::cube BA(grid_y, grid_x, grid_z);
	arma::cube r_on(grid_y, grid_x, grid_z);
	arma::cube r_off(grid_y, grid_x, grid_z);
	double alpha0 = 0; // total propensity
	
	// total number of tiles (according to the shape matrix)
	int tiles_total = sum(shapeMatVec);
	
	// equilibrum values
	double F_eq = 0;
	double C_eq = 0;
	// diassociation rate
	double kd = 0;
	// equilibirum under F <-> C
	if(k_on != 0 && k_off != 0){		
		kd = k_off / k_on; 
		F_eq = round(F_total/tiles_total/(kd+1)*kd);
		C_eq = round(F_total/tiles_total) - F_eq;
	} else {
		F_eq = round(F_total/tiles_total);
		C_eq = 0;
	}
	
	int last_frame = 0;
	// calculate the initial distributions and propensities
	for(int i=0;i<grid_y;i++){
		for(int j=0;j<grid_x;j++){
			for(int k=0;k<grid_z;k++){
				if(shapeMat(i,j,k) != 0){
					// numbers of molecules
					F(i,j,k) = F_eq;
					C(i,j,k) = C_eq;
					FL(i,j,k) = F_eq;
					FLC(i,j,k) = C_eq;				
				} else {
					F(i,j,k) = 0;
					C(i,j,k) = 0;
					FL(i,j,k) = 0;
					FLC(i,j,k) = 0;
				}
				// propensities			
				LR(i,j,k) = D_LR(i,j,k) * F(i,j,k);
				RL(i,j,k) = D_RL(i,j,k) * F(i,j,k);
				UD(i,j,k) = D_UD(i,j,k) * F(i,j,k);
				DU(i,j,k) = D_DU(i,j,k) * F(i,j,k);
				AB(i,j,k) = D_AB(i,j,k) * F(i,j,k);
				BA(i,j,k) = D_BA(i,j,k) * F(i,j,k);

				// reaction propensities
				r_on(i,j,k) = k_on * F(i,j,k);
				r_off(i,j,k) = k_off * C(i,j,k);
			
				// update total
				alpha0 += LR(i,j,k) + RL(i,j,k) + UD(i,j,k) + DU(i,j,k) + AB(i,j,k) + BA(i,j,k) + r_on(i,j,k) + r_off(i,j,k);			
			}
		}
	}
	
	double dt = 0;
	double action_val=0;
	double action_sum, action_sum_next, alpha_diff, fl_prop;
	Rcpp::List FLOR((int)floor(max_time*fps)+1);
	// iterate over time and simulate!
	while(t <= max_time){	
		Rcpp::NumericVector r = runif(2);
		
		r1 = r[0];
		r2 = r[1];
		
		// (1) calculate new time
		dt = 1.0 / alpha0 * log( 1.0 / r1);
		
		// (2) get the new event
		action_val = r2 * alpha0;
		
		// re-calculate the sum to see which action it is... 
		// TOTRY: in principle it could be also done by a binary search... 
		action_sum = 0;
		action_sum_next = 0;
		int action = 0; // identifies the current action (debug only)
		bool finished = false;
		for(int i=0;i<grid_y;i++){
			for(int j=0;j<grid_x;j++){
				for(int k=0;k<grid_z;k++){
					action = 0;
					action_sum_next = action_sum;
					action_sum_next += LR(i,j,k);				
					DIFFUSION3D(LR, i, j, k, i, j+1, k);
					action++;
				
					action_sum = action_sum_next;
					action_sum_next += RL(i,j,k);
					DIFFUSION3D(RL, i, j, k, i, j-1, k);
					action++;
				
					action_sum = action_sum_next;
					action_sum_next += UD(i,j,k);
					DIFFUSION3D(UD, i, j, k, i+1, j, k);
					action++;
				
					action_sum = action_sum_next;
					action_sum_next += DU(i,j,k);
					DIFFUSION3D(DU, i, j, k, i-1, j, k);
					action++;
					
					action_sum = action_sum_next;
					action_sum_next += AB(i,j,k);
					DIFFUSION3D(AB, i, j, k, i, j, k+1);
					action++;
				
					action_sum = action_sum_next;
					action_sum_next += BA(i,j,k);
					DIFFUSION3D(BA, i, j, k, i, j, k-1);
					action++;
				
					action_sum = action_sum_next;
					action_sum_next += r_on(i,j,k);
					REACTION(r_on, FL, F, -1, 1);
					action++;
				
					action_sum = action_sum_next;
					action_sum_next += r_off(i,j,k);
					REACTION(r_off, FLC, C, 1, -1);
					action++;
				
					action_sum = action_sum_next;
				}
				
				// if action was performed break out of everything... 
				if(finished)
					break;
				
			}
			
			// if action was performed break out of everything... 
			if(finished){
				break;
			}
		}
		
		/*if(abs(alpha0-action_sum) > 1e-8)
			return(wrap(false)); */
		
		if(CONSISTENCY_CHECK_3D){
			bool broken = false;
			double alpha0rep = 0;
			for(int i=0;i<grid_y;i++){
				for(int j=0;j<grid_x;j++){
					for(int k=0;k<grid_z;k++){
						// propensities			
						if(LR(i,j,k) != D_LR(i,j,k) * F(i,j,k))
							broken = true;
						if(RL(i,j,k) != D_RL(i,j,k) * F(i,j,k))
							broken = true;
						if(UD(i,j,k) != D_UD(i,j,k) * F(i,j,k))
							broken = true;
						if(DU(i,j,k) != D_DU(i,j,k) * F(i,j,k))
							broken = true;
						if(AB(i,j,k) != D_AB(i,j,k) * F(i,j,k))
							broken = true;
						if(BA(i,j,k) != D_BA(i,j,k) * F(i,j,k))
							broken = true;

						// reaction propensities
						if(r_on(i,j,k) != k_on * F(i,j,k))
							broken = true;
						if(r_off(i,j,k) != k_off * C(i,j,k))
							broken = true;
						
						// absolute numbers
						if(F(i,j,k) < 0 || C(i,j,k) < 0 || FL(i,j,k) < 0 || FLC(i,j,k) < 0)
							broken = true;
						
						alpha0rep += LR(i,j,k) + RL(i,j,k) + UD(i,j,k) + DU(i,j,k) + AB(i,j,k) + BA(i,j,k) + r_on(i,j,k) + r_off(i,j,k);						
					}
				}
			}
			
			if(abs(alpha0-alpha0rep) > 1e-8)
				broken = true;
			
			if(broken){
				Rcpp::List out;
			
				out["F"] = cloneCube(F);
				out["FL"] = cloneCube(FL);
				out["C"] = cloneCube(C);
				out["FLC"] = cloneCube(FLC);
				out["t"] = t;
			
				// debug stuffs
				out["LR"] = cloneCube(LR);
				out["RL"] = cloneCube(RL);
				out["UD"] = cloneCube(UD);				
				out["DU"] = cloneCube(DU);
				out["AB"] = cloneCube(AB);
				out["BA"] = cloneCube(BA);
			
				out["D_LR"] = cloneCube(D_LR);
				out["D_RL"] = cloneCube(D_RL);
				out["D_UD"] = cloneCube(D_UD);
				out["D_DU"] = cloneCube(D_DU);
				out["D_AB"] = cloneCube(D_AB);
				out["D_BA"] = cloneCube(D_BA);
				out["r_on"] = cloneCube(r_on);
				out["r_off"] = cloneCube(r_off);
			
				out["alpha0"] = alpha0;
				out["action"] = action;
				out["error"] = true;
				
				return(out);
			}
		}
		
		
		// update time
		t = t + dt;	
		
		// bleach if needed!!!
		if(!frapped && t > frap_time){			
			for(int i=0;i<grid_y;i++){
				for(int j=0;j<grid_x;j++){
					for(int k=0;k<grid_z;k++){
						if(frapMask(i,j,k)){
							// bleach!
							FL(i,j,k) = round(F(i,j,k) * bleach_depth);
							FLC(i,j,k) = round(C(i,j,k) * bleach_depth);
						}
					}
				}
			}
			frapped = true;
		}
				
		//std::cout << "t: " << t << " fps: " << fps << " last_frame: " << last_frame << " num: " << (int)floor(t / (1.0/fps)) << "\n";
		
		// plot and output current simulation time		
		if((int)floor(t / (1.0/fps)) != last_frame){
			last_frame = floor(t / (1.0/fps));

			// save a snapshot of time timepoint			
			Rcpp::List out;
			
			out["F"] = cloneCube(F);
			out["FL"] = cloneCube(FL);
			out["C"] = cloneCube(C);
			out["FLC"] = cloneCube(FLC);
			out["t"] = t;
			
			// debug stuffs
			if(DEBUG_3D){
				out["LR"] = cloneCube(LR);
				out["RL"] = cloneCube(RL);
				out["UD"] = cloneCube(UD);
				out["DU"] = cloneCube(DU);
				out["AB"] = cloneCube(AB);
				out["BA"] = cloneCube(BA);
				
				out["D_LR"] = cloneCube(D_LR);
				out["D_RL"] = cloneCube(D_RL);
				out["D_UD"] = cloneCube(D_UD);
				out["D_DU"] = cloneCube(D_DU);
				out["D_AB"] = cloneCube(D_AB);
				out["D_BA"] = cloneCube(D_BA);
				out["r_on"] = cloneCube(r_on);
				out["r_off"] = cloneCube(r_off);
				
				out["alpha0"] = alpha0;
			}
			
			FLOR[last_frame-1] = out;
		
			std::cout << "Time: " << t << "\n";
		}
		
		// update event indicator			
		event_inx++;
	}
	
	return(FLOR);		
}



