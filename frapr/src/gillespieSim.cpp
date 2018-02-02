#include <Rcpp.h>

#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>

#define DEBUG false
#define CONSISTENCY_CHECK false

/*
	This macro is used as a general template for diffusion operations
	
	@param PROP propensity matrix, e.g. LR
	@param i current i value
	@param j current j value
	@param ni i coordinate that is modified
	@param nj j coordinate that is modified
*/
#define DIFFUSION(PROP, I, J, NI, NJ) ({ \
	if(PROP(I,J) != 0 && action_sum < action_val && action_val <= action_sum_next){ \
		/* how many are fluorescent */ \
		fl_prop = FL(I,J) / F(I,J); \
		\
		/* calculate the difference in total propensities */ \
		alpha_diff = - LR(I,J) - RL(I,J) - UD(I,J) - DU(I,J) - r_on(I,J) - \
		               LR(NI,NJ) - RL(NI,NJ) - UD(NI,NJ) - DU(NI,NJ) - r_on(NI,NJ); \
		\
		/* update */ \
		F(I,J) -= 1; \
		F(NI,NJ) += 1; \
		/* make sure to update all propensities!!! */ \
		LR(I,J) = D_LR(I,J) * F(I,J); \
		RL(I,J) = D_RL(I,J) * F(I,J); \
		UD(I,J) = D_UD(I,J) * F(I,J); \
		DU(I,J) = D_DU(I,J) * F(I,J); \
		r_on(I,J) = k_on * F(I,J); \
		/* and for ni,nj */\
		LR(NI,NJ) = D_LR(NI,NJ) * F(NI,NJ); \
		RL(NI,NJ) = D_RL(NI,NJ) * F(NI,NJ); \
		UD(NI,NJ) = D_UD(NI,NJ) * F(NI,NJ); \
		DU(NI,NJ) = D_DU(NI,NJ) * F(NI,NJ); \
		r_on(NI,NJ) = k_on * F(NI,NJ);\
		\
		/* add new values to propensitiy */ \
		alpha_diff +=  LR(I,J) + RL(I,J) + UD(I,J) + DU(I,J) + r_on(I,J) + \
		               LR(NI,NJ) + RL(NI,NJ) + UD(NI,NJ) + DU(NI,NJ) + r_on(NI,NJ); \
		\
		/* update the fluorescent ones as well */ \
		if(runif(1)[0] < fl_prop){ \
			FL(I,J) -= 1; \
			FL(NI,NJ) += 1; \
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
	if(PROP(i,j) != 0 && action_sum < action_val && action_val <= action_sum_next){ \
		/* how many are fluorescent */ \
		fl_prop = FLOR(i,j) / FLOR_TOTAL(i,j); \
		\
		/* calculate the difference in total propensities */ \
		alpha_diff = - r_on(i,j) - r_off(i,j) - LR(i,j) - RL(i,j) - UD(i,j) - DU(i,j); \
		\
		/* update */ \
		F(i,j) += F_ADD; \
		C(i,j) += C_ADD; \
		LR(i,j) = D_LR(i,j) * F(i,j); \
		RL(i,j) = D_RL(i,j) * F(i,j); \
		UD(i,j) = D_UD(i,j) * F(i,j); \
		DU(i,j) = D_DU(i,j) * F(i,j); \
		r_on(i,j) = F(i,j) * k_on; \
		r_off(i,j) = C(i,j) * k_off; \
		\
		/* add new values to propensitiy */ \
		alpha_diff += r_on(i,j) + r_off(i,j) + LR(i,j) + RL(i,j) + UD(i,j) + DU(i,j); \
		\
		/* update the fluorescent ones as well */ \
		if(runif(1)[0] < fl_prop){ \
			FL(i,j) += F_ADD; \
			FLC(i,j) += C_ADD; \
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

RcppExport SEXP gillespieSimCpp(SEXP inGridX, SEXP inGridY, SEXP inFtotal, SEXP inKon, SEXP inKoff, SEXP inL, SEXP inD, SEXP inShapeMat, SEXP inMaxTime, 
	SEXP inFrapTime, SEXP inFrapMask, SEXP inBleachDepth, SEXP inFps) {
	// initialise the input parameters
	Rcpp::IntegerVector rGridX(inGridX);
	Rcpp::IntegerVector rGridY(inGridY);
	Rcpp::IntegerVector rFtotal(inFtotal);
	Rcpp::NumericVector rkOn(inKon);
	Rcpp::NumericVector rkOff(inKoff);
	Rcpp::NumericVector rL(inL);
	Rcpp::NumericVector rD(inD);
	Rcpp::NumericMatrix shapeMat(inShapeMat);
	Rcpp::NumericVector rMaxTime(inMaxTime);
	Rcpp::NumericVector rFrapTime(inFrapTime);
	Rcpp::NumericMatrix frapMask(inFrapMask);
	Rcpp::NumericVector rFps(inFps);
	Rcpp::NumericVector rBleachDepth(inBleachDepth);
	
	// convert into normal C++ data structures
	int grid_x = rGridX[0];
	int grid_y = rGridY[0];
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

	// calculate diffusion coefficients
	double dx = D / (hx*hx);
	double dy = D / (hy*hx);
	
	// create the diffusion matrices for the 4 directions
	Rcpp::NumericMatrix D_LR(grid_y, grid_x);
	Rcpp::NumericMatrix D_RL(grid_y, grid_x);
	Rcpp::NumericMatrix D_UD(grid_y, grid_x);
	Rcpp::NumericMatrix D_DU(grid_y, grid_x);
	
	for(int i=0;i<grid_y;i++){
		for(int j=0;j<grid_x;j++){
			// LR
			if(j == (grid_x - 1) || shapeMat(i,j) == 0 || shapeMat(i,j+1) == 0){
				D_LR(i,j) = 0;
			} else {
				D_LR(i,j) = dx;
			}
			
			// RL
			if(j == 0 || shapeMat(i,j) == 0 || shapeMat(i,j-1) == 0){
				D_RL(i,j) = 0;
			} else {
				D_RL(i,j) = dx;
			}
			
			// UD
			if(i == (grid_y - 1) || shapeMat(i,j)==0 || shapeMat(i+1,j) == 0){
				D_UD(i,j) = 0;
			} else {
				D_UD(i,j) = dy;
			}
			
			// DU
			if(i == 0 || shapeMat(i,j)==0 || shapeMat(i-1,j) == 0){
				D_DU(i,j) = 0;
			} else {
				D_DU(i,j) = dy;
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
	Rcpp::NumericMatrix F(grid_y, grid_x);
	Rcpp::NumericMatrix FL(grid_y, grid_x);
	Rcpp::NumericMatrix C(grid_y, grid_x);
	Rcpp::NumericMatrix FLC(grid_y, grid_x);
	
	// propensity matrices
	Rcpp::NumericMatrix LR(grid_y, grid_x);
	Rcpp::NumericMatrix RL(grid_y, grid_x);
	Rcpp::NumericMatrix UD(grid_y, grid_x);
	Rcpp::NumericMatrix DU(grid_y, grid_x);
	Rcpp::NumericMatrix r_on(grid_y, grid_x);
	Rcpp::NumericMatrix r_off(grid_y, grid_x);
	double alpha0 = 0; // total propensity
	
	// total number of tiles (according to the shape matrix)
	int tiles_total = sum(shapeMat);
	
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
			if(shapeMat(i,j) != 0){
				// numbers of molecules
				F(i,j) = F_eq;
				C(i,j) = C_eq;
				FL(i,j) = F_eq;
				FLC(i,j) = C_eq;				
			} else {
				F(i,j) = 0;
				C(i,j) = 0;
				FL(i,j) = 0;
				FLC(i,j) = 0;
			}
			// propensities			
			LR(i,j) = D_LR(i,j) * F(i,j);
			RL(i,j) = D_RL(i,j) * F(i,j);
			UD(i,j) = D_UD(i,j) * F(i,j);
			DU(i,j) = D_DU(i,j) * F(i,j);

			// reaction propensities
			r_on(i,j) = k_on * F(i,j);
			r_off(i,j) = k_off * C(i,j);
			
			// update total
			alpha0 += LR(i,j) + RL(i,j) + UD(i,j) + DU(i,j) + r_on(i,j) + r_off(i,j);			
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
		int action = 0; // identifies the current action
		bool finished = false;
		for(int i=0;i<grid_y;i++){
			for(int j=0;j<grid_x;j++){
				action = 0;
				action_sum_next = action_sum;
				action_sum_next += LR(i,j);				
				DIFFUSION(LR, i, j, i, j+1);
				action++;
				
				action_sum = action_sum_next;
				action_sum_next += RL(i,j);
				DIFFUSION(RL, i, j, i, j-1);
				action++;
				
				action_sum = action_sum_next;
				action_sum_next += UD(i,j);
				DIFFUSION(UD, i, j, i+1, j);
				action++;
				
				action_sum = action_sum_next;
				action_sum_next += DU(i,j);
				DIFFUSION(DU, i, j, i-1, j);
				action++;
				
				action_sum = action_sum_next;
				action_sum_next += r_on(i,j);
				REACTION(r_on, FL, F, -1, 1);
				action++;
				
				action_sum = action_sum_next;
				action_sum_next += r_off(i,j);
				REACTION(r_off, FLC, C, 1, -1);
				action++;
				
				action_sum = action_sum_next;
				
			}
			
			// if action was performed break out of everything... 
			if(finished)
				break;
		}
		
		/*if(abs(alpha0-action_sum) > 1e-8)
			return(wrap(false)); */
		
		if(CONSISTENCY_CHECK){
			bool broken = false;
			double alpha0rep = 0;
			for(int i=0;i<grid_y;i++){
				for(int j=0;j<grid_x;j++){
					// propensities			
					if(LR(i,j) != D_LR(i,j) * F(i,j))
						broken = true;
					if(RL(i,j) != D_RL(i,j) * F(i,j))
						broken = true;
					if(UD(i,j) != D_UD(i,j) * F(i,j))
						broken = true;
					if(DU(i,j) != D_DU(i,j) * F(i,j))
						broken = true;

					// reaction propensities
					if(r_on(i,j) != k_on * F(i,j))
						broken = true;
					if(r_off(i,j) != k_off * C(i,j))
						broken = true;
						
					// absolute numbers
					if(F(i,j) < 0 || C(i,j) < 0 || FL(i,j) < 0 || FLC(i,j) < 0)
						broken = true;
						
					alpha0rep += LR(i,j) + RL(i,j) + UD(i,j) + DU(i,j) + r_on(i,j) + r_off(i,j);						
				}
			}
			
			if(abs(alpha0-alpha0rep) > 1e-8)
				broken = true;
			
			if(broken){
				Rcpp::List out;
			
				out["F"] = clone(F);
				out["FL"] = clone(FL);
				out["C"] = clone(C);
				out["FLC"] = clone(FLC);
				out["t"] = t;
			
				// debug stuffs
				out["LR"] = clone(LR);
				out["RL"] = clone(RL);
				out["UD"] = clone(UD);
				out["DU"] = clone(DU);
			
				out["D_LR"] = clone(D_LR);
				out["D_RL"] = clone(D_RL);
				out["D_UD"] = clone(D_UD);
				out["D_DU"] = clone(D_DU);
				out["r_on"] = clone(r_on);
				out["r_off"] = clone(r_off);
			
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
					if(frapMask(i,j)){
						// bleach!
						FL(i,j) = round(F(i,j) * bleach_depth);
						FLC(i,j) = round(C(i,j) * bleach_depth);
					}
				}
			}
			frapped = true;
		}
				
		// std::cout << "t: " << t << " fps: " << fps << " last_frame: " << last_frame << " num: " << (int)floor(t / (1.0/fps)) << "\n";
		
		// plot and output current simulation time		
		if((int)floor(t / (1.0/fps)) != last_frame){
			last_frame = floor(t / (1.0/fps));

			// save a snapshot of time timepoint			
			Rcpp::List out;
			
			out["F"] = clone(F);
			out["FL"] = clone(FL);
			out["C"] = clone(C);
			out["FLC"] = clone(FLC);
			out["t"] = t;
			
			// debug stuffs
			if(DEBUG){
				out["LR"] = clone(LR);
				out["RL"] = clone(RL);
				out["UD"] = clone(UD);
				out["DU"] = clone(DU);
				
				out["D_LR"] = clone(D_LR);
				out["D_RL"] = clone(D_RL);
				out["D_UD"] = clone(D_UD);
				out["D_DU"] = clone(D_DU);
				out["r_on"] = clone(r_on);
				out["r_off"] = clone(r_off);
				
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



