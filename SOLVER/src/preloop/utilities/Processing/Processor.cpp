// class to handle various processing operations 
// e.g. filter, taper, frequency derivative, convolution, ... 


#include "Processor.h"
#include "KernerFFTW_N3.h"
#include "KernerFFTW_N6.h"
#include "KernerFFTW_N9.h"
#include "PreloopFFTW_time.h"
#include "eigenc.h"
#include "Tapers.h"
#include "Filters.h"

RColX Processor::sTime;
RColX Processor::sFreq;
RMatXX Processor::sFilters;
RColX Processor::sTaper;
Real Processor::sT = 0.;
Real Processor::sDt = 0.;
Real Processor::sDf = 0.;
Real Processor::sWindowBeg = 0.;
Real Processor::sWindowEnd = 0.;
int Processor::sNumFilters = 0;

void Processor::initialize(int totSteps, const RColX &bufTime, const RMatX2 filtParams) {

	// get dt 
	sTime = bufTime;
	sDt = bufTime(1) - bufTime(0); //sampling 
	int temp_size = bufTime.size();
	Real temp_T = bufTime(temp_size - 1); //current max Time 

	// create taper with current time 
	Tapers::cosineTaper(sTaper, temp_size);

	// prolong taper 
	zeroPad(sTaper, totSteps);

	// prolong time 
	zeroPad(sTime, totSteps);
	
	for (int i = 0; i < totSteps - temp_size; i++) {
		sTime(temp_size + i) = temp_T + (i+1) * sDt;
	}
	
	// new total time
	sT = sTime(totSteps - 1) - sTime(0); 
	
	//create freq 
	sFreq = RColX(totSteps/2 + 1);
	sDf = one/sT;
	
	for (int i = 0; i < totSteps/2 + 1; i++) {
		sFreq(i) = i * sDf;
	}
	
	// create filters 
	sNumFilters = filtParams.rows();
	sFilters = RMatXX(sNumFilters, sFreq.size());
	
	for (int i = 0; i< sNumFilters; i++) {
		Filters::logGabor(sFilters, sFreq, one/filtParams(i,0), filtParams(i,1), i);
	}
	
}

void Processor::finalize() {

	
}

void Processor::taper(vec_vec_ar3_CMatPP &u) {
	
	if (u.size() != sTaper.size()) throw std::runtime_error("Processor::taper error : trace and taper have different sizes.");

	for (int it = 0; it < u.size(); it++)
		for (int inu = 0; inu < u[0].size(); inu++)
			for (int ic = 0; ic < 3; ic++)
				u[it][inu][ic] *= sTaper(it);
	
}

void Processor::zeroPad(RColX &trace, int npad) {
	
	int temp_size = trace.size();
	trace.conservativeResize( npad );
	
	for (int i = temp_size; i < trace.size(); i++) {
		trace(i) = zero;
	}
}



void Processor::transformT2F(const vec_ar3_RMatPP& ut, vec_ar3_CMatPP& uf) {
	

	RMatXN3 &tempT = KernerFFTW_N3::getR2C_RMat();
	makeFlat<vec_ar3_RMatPP, RMatXN3>(ut, tempT);
	KernerFFTW_N3::computeR2C();
	CMatXN3 tempF = KernerFFTW_N3::getR2C_CMat();
	makeStruct<vec_ar3_CMatPP, CMatXN3>(uf, tempF);

}

void Processor::transformT2F(const vec_ar6_RMatPP& ut, vec_ar6_CMatPP& uf) {
	

	RMatXN6 &tempT = KernerFFTW_N6::getR2C_RMat();
	makeFlat<vec_ar6_RMatPP, RMatXN6>(ut, tempT);
	KernerFFTW_N6::computeR2C();
	CMatXN6 tempF = KernerFFTW_N6::getR2C_CMat();
	makeStruct<vec_ar6_CMatPP, CMatXN6>(uf, tempF);
	

}

void Processor::transformT2F(const vec_ar9_RMatPP& ut, vec_ar9_CMatPP& uf) {
	

	RMatXN9 &tempT = KernerFFTW_N9::getR2C_RMat();

	makeFlat<vec_ar9_RMatPP, RMatXN9>(ut, tempT);
	KernerFFTW_N9::computeR2C();
	CMatXN9 tempF = KernerFFTW_N9::getR2C_CMat();
	makeStruct<vec_ar9_CMatPP, CMatXN9>(uf, tempF);

}

void Processor::transformF2T(const vec_ar3_CMatPP& uf, vec_ar3_RMatPP& ut) {
	
	CMatXN3 &tempF = KernerFFTW_N3::getC2R_CMat();
	makeFlat<vec_ar3_CMatPP, CMatXN3>(uf, tempF);
	KernerFFTW_N3::computeC2R();
	RMatXN3 tempT = KernerFFTW_N3::getC2R_RMat();
	makeStruct<vec_ar3_RMatPP, RMatXN3>(ut, tempT);

}

void Processor::transformF2T(const vec_ar6_CMatPP& uf, vec_ar6_RMatPP& ut) {

	CMatXN6 &tempF = KernerFFTW_N6::getC2R_CMat();
	makeFlat<vec_ar6_CMatPP, CMatXN6>(uf, tempF);
	KernerFFTW_N6::computeC2R();
	RMatXN6 tempT = KernerFFTW_N6::getC2R_RMat();
	makeStruct<vec_ar6_RMatPP, RMatXN6>(ut, tempT);

}

void Processor::transformF2T(const vec_ar9_CMatPP& uf, vec_ar9_RMatPP& ut) {

	CMatXN9 &tempF = KernerFFTW_N9::getC2R_CMat(); 
	makeFlat<vec_ar9_CMatPP, CMatXN9>(uf, tempF);
	KernerFFTW_N9::computeC2R();
	RMatXN9 tempT = KernerFFTW_N9::getC2R_RMat();
	makeStruct<vec_ar9_RMatPP, RMatXN9>(ut, tempT);

}

void Processor::transformT2F_preloop(const RMatX3& ut, CMatX3& uf) {
	
	RMatX3 &tempT = PreloopFFTW_time::getR2C_RMat();
	tempT = ut;
	PreloopFFTW_time::computeR2C();
	uf = PreloopFFTW_time::getR2C_CMat();
		
}

void Processor::transformF2T_preloop(const CMatX3& uf, RMatX3& ut) {
	
	CMatX3 &tempF = PreloopFFTW_time::getC2R_CMat();
	tempF = uf;
	PreloopFFTW_time::computeC2R();
	ut = PreloopFFTW_time::getC2R_RMat();
		
}