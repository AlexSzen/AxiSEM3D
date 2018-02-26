// class to handle various processing operations 
// e.g. filter, taper, frequency derivative, convolution, ... 

#pragma once 


#include "eigenc.h"
#include "eigenp.h"
#include "XMPI.h"
#include <iostream>
class Processor {

public:	
	// we can initialize once in preloop, and re initialize during kernels 
	static void initialize(int totSteps, const RColX &bufTime);
	static void finalize();
	
	
	static void zeroPad(RColX &trace, int npad);
	static void taper(vec_vec_ar3_CMatPP &u);
	static void taper(RMatX3 &trace, Real begWin, Real endWin);
	static void createFilters(const RMatXX &filter_params);
	
	// processing for kernels 
	static void transformT2F(const vec_ar3_RMatPP& ut, vec_ar3_CMatPP& uf);
	static void transformT2F(const vec_ar6_RMatPP& ut, vec_ar6_CMatPP& uf);
	static void transformT2F(const vec_ar9_RMatPP& ut, vec_ar9_CMatPP& uf);
	
	static void transformF2T(const vec_ar3_CMatPP& uf, vec_ar3_RMatPP& ut);
	static void transformF2T(const vec_ar6_CMatPP& uf, vec_ar6_RMatPP& ut);
	static void transformF2T(const vec_ar9_CMatPP& uf, vec_ar9_RMatPP& ut);
	
	// processing for seismograms 
	static void transformT2F(const RMatX3 &ut, CMatX3 &uf);
	static void transformF2T(const CMatX3 &uf, RMatX3 &ut);



	template<class vec_arY_CMatPP>
	static void filter(vec_arY_CMatPP &uf, int ifilt) {
		
		if (sFilters.cols() != uf.size()) throw std::runtime_error("Processor::filter error : filter and trace of different lengths.");
		
		for (int i = 0; i < uf.size(); i++)
			for (int ic = 0; ic < uf[0].size(); ic++)
				uf[i][ic] *= sFilters(ifilt,i);
	};
	
	static void filter(CMatX3 &trace, int ifilter) {
		
		if (sFilters.cols() != trace.rows()) throw std::runtime_error("Processor::filter error : filter and trace of different lengths.");
		for (int i = 0; i < trace.rows(); i++) {
			for (int ic = 0; ic < trace.cols(); ic++) {
				trace(i,ic) *= sFilters(ifilter, i);
			}			
		}
		
	}
	// take (time) derivative in frequency domain, i.e multiply by 2*pi*i*f
	template<class vec_arY_CMatPP>
	static void derivate(vec_arY_CMatPP &uf) {
		
		if (sFreq.size() != uf.size()) throw std::runtime_error("Processor::derivate error : frequency and trace of different lengths.");
		
		for (int i = 0; i < uf.size(); i++)
			for (int ic = 0; ic < uf[0].size(); ic++)
				uf[i][ic] *= two * (Real) pi * ii * sFreq(i);
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