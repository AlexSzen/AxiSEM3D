// class to handle various processing operations 
// e.g. filter, taper, frequency derivative, convolution, ... 

#pragma once 


#include "eigenc.h"
#include "eigenp.h"
#include "XMPI.h"
#include <iostream>
class Processor {

public:	

	static void initialize(int totSteps, const RColX &bufTime, const RMatX2 filtParams, Real begWin, Real endWin);
	static void finalize();
	
	
	static void zeroPad(RColX &trace, int npad);
	static void taper(vec_vec_ar3_CMatPP &u);
	
	static void transformT2F(const vec_ar3_RMatPP& ut, vec_ar3_CMatPP& uf);
	static void transformT2F(const vec_ar6_RMatPP& ut, vec_ar6_CMatPP& uf);
	static void transformT2F(const vec_ar9_RMatPP& ut, vec_ar9_CMatPP& uf);
	
	static void transformF2T(const vec_ar3_CMatPP& uf, vec_ar3_RMatPP& ut);
	static void transformF2T(const vec_ar6_CMatPP& uf, vec_ar6_RMatPP& ut);
	static void transformF2T(const vec_ar9_CMatPP& uf, vec_ar9_RMatPP& ut);

	template<class vec_arY_CMatPP>
	static void filter(vec_arY_CMatPP &uf, int ifilt) {
		
		if (sFilters.cols() != uf.size()) throw std::runtime_error("Processor::filter error : filter and trace of different lengths.");
		
		for (int i = 0; i < uf.size(); i++)
			for (int ic = 0; ic < uf[0].size(); ic++)
				uf[i][ic] *= sFilters(ifilt,i);
	};

	template<class vec_arY_CMatPP> 
	static void sumAndConvolve(const vec_arY_CMatPP &uf1, const vec_arY_CMatPP &uf2, vec_ar6_CMatPP &conv, const IColX &inds1, const IColX &inds2, Real prefactor, int indc) {
	// conv specifically for kernels. ufs are fwd and bwd freq domain for 1 slice
	// inds1 are indices from uf1 to be summed, then convolved with summed inds from uf2
	// conv is size 6 because 6 base kernels in elastic medium with radial anisotropy, Fichtner p.169 
	
		if (uf1.size() != uf2.size()) throw std::runtime_error("Processor::convolveSum : input sizes differ");
		
		vec_CMatPP temp1(uf1.size(), CMatPP::Zero());
		vec_CMatPP temp2(uf1.size(), CMatPP::Zero());
		
		for (int i = 0; i < uf1.size(); i++) {
			
			for (int i1 = 0; i1 < inds1.size(); i1++) {
				temp1[i] += uf1[i][inds1(i1)]; 
			}
			
			for (int i2 = 0; i2 < inds2.size(); i2++) {
				temp2[i] += uf2[i][inds2(i2)]; 
			}
			
		}
		
		for (int i = 0; i < uf1.size(); i++) {
			conv[i][indc] += prefactor * temp1[i].schur(temp2[i]); //we use += because for some kernels convolution is done in several calls to sumAndConvolve. 
		}
		
		
		
	};
	
	template<class vec_vec_arY_CMatPP, class vec_arY_CMatPP>
	static void timeWindow(const vec_vec_arY_CMatPP &utf, vec_arY_CMatPP &uf) {
		
		for (int it = 0; it < sTime.size(); it++) 
			if (sTime(it) > sWindowBeg && sTime(it) < sWindowEnd) 
				for (int inu = 0; inu < uf.size(); inu ++) 
					for (int ic = 0; ic < uf[0].size(); ic++)
						uf[inu][ic] += utf[it][inu][ic];
				
			
		
		
	};
	
	
	template<class vec_arY_TMatPP, class TMatXNY>
	static void makeFlat(const vec_arY_TMatPP &ucStruct, TMatXNY &ucFlat) {
		for (int alpha = 0; alpha < ucStruct.size(); alpha++) {
			for (int i = 0; i < ucStruct[0].size(); i++) {
				for (int j = 0; j < nPntEdge; j++) {
					ucFlat.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge) 
						= ucStruct[alpha][i].block(j, 0, 1, nPntEdge);
				}
			} 
		} 
	};
	
	template<class vec_arY_TMatPP, class TMatXNY>
	static void makeStruct(vec_arY_TMatPP &ucStruct, const TMatXNY &ucFlat) {
		for (int alpha = 0; alpha < ucStruct.size(); alpha++) {
			for (int i = 0; i < ucStruct[0].size(); i++) {
				for (int j = 0; j < nPntEdge; j++) {
					ucStruct[alpha][i].block(j, 0, 1, nPntEdge)
						= ucFlat.block(alpha, nPE * i + nPntEdge * j, 1, nPntEdge);
				}
			} 
		} 
	};
	
	static int sNumFilters;
private:
	static RColX sTime;
	static RColX sFreq;
	static RMatXX sFilters;
	static RColX sTaper;
	
	static Real sDt, sDf, sT;
	
	static Real sWindowBeg, sWindowEnd;
	
	
};