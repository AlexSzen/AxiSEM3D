// class to handle various processing operations 
// e.g. filter, taper, frequency derivative, convolution, ... 


#include "Processor.h"
#include "KernerFFTW_N3.h"
#include "KernerFFTW_N6.h"
#include "KernerFFTW_N9.h"
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

void Processor::initialize(int totSteps, const RColX &bufTime, const RMatX2 filtParams, Real begWin, Real endWin) {
	
	// init FFTW
	KernerFFTW_N3::initialize(totSteps);
	KernerFFTW_N6::initialize(totSteps);
	KernerFFTW_N9::initialize(totSteps);
	
	// get dt 
	sTime = bufTime;
	sDt = bufTime(1) - bufTime(0); //sampling 
	int temp_size = bufTime.size();
	Real temp_T = bufTime(temp_size - 1); //current max Time 

	// create taper with current time 
	Tapers::cosineTaper(sTaper, temp_size);
	
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
	
	// window 
	sWindowBeg = begWin;
	sWindowEnd = endWin;

	
}

void Processor::finalize() {
	
	KernerFFTW_N3::finalize();
	KernerFFTW_N6::finalize();
	KernerFFTW_N9::finalize();
	
}

void Processor::zeroPad(RColX &trace, int npad) {
	
	trace.conservativeResize( npad );
		
}


void Processor::transformT2F(const vec_ar3_RMatPP& ut, vec_ar3_CMatPP& uf) {
	

	RMatXN3 &tempT = KernerFFTW_N3::getR2C_RMat();
	makeFlat<vec_ar3_RMatPP, RMatXN3>(ut, tempT);
	KernerFFTW_N3::computeR2C();
	CMatXN3 &tempF = KernerFFTW_N3::getR2C_CMat();
	makeStruct<vec_ar3_CMatPP, CMatXN3>(uf, tempF);
	
	
}

void Processor::transformT2F(const vec_ar6_RMatPP& ut, vec_ar6_CMatPP& uf) {
	

	RMatXN6 &tempT = KernerFFTW_N6::getR2C_RMat();
	makeFlat<vec_ar6_RMatPP, RMatXN6>(ut, tempT);
	KernerFFTW_N6::computeR2C();
	CMatXN6 &tempF = KernerFFTW_N6::getR2C_CMat();
	makeStruct<vec_ar6_CMatPP, CMatXN6>(uf, tempF);
	
	
}

void Processor::transformT2F(const vec_ar9_RMatPP& ut, vec_ar9_CMatPP& uf) {
	

	RMatXN9 &tempT = KernerFFTW_N9::getR2C_RMat();
	makeFlat<vec_ar9_RMatPP, RMatXN9>(ut, tempT);
	KernerFFTW_N9::computeR2C();
	CMatXN9 &tempF = KernerFFTW_N9::getR2C_CMat();
	makeStruct<vec_ar9_CMatPP, CMatXN9>(uf, tempF);
	
}

void Processor::transformF2T(const vec_ar3_CMatPP& uf, vec_ar3_RMatPP& ut) {
	
	CMatXN3 &tempF = KernerFFTW_N3::getC2R_CMat();
	makeFlat<vec_ar3_CMatPP, CMatXN3>(uf, tempF);
	KernerFFTW_N3::computeR2C();
	RMatXN3 &tempT = KernerFFTW_N3::getC2R_RMat();
	makeStruct<vec_ar3_RMatPP, RMatXN3>(ut, tempT);
	
}

void Processor::transformF2T(const vec_ar6_CMatPP& uf, vec_ar6_RMatPP& ut) {
	
	CMatXN6 &tempF = KernerFFTW_N6::getC2R_CMat();
	makeFlat<vec_ar6_CMatPP, CMatXN6>(uf, tempF);
	KernerFFTW_N6::computeR2C();
	RMatXN6 &tempT = KernerFFTW_N6::getC2R_RMat();
	makeStruct<vec_ar6_RMatPP, RMatXN6>(ut, tempT);
	
}

void Processor::transformF2T(const vec_ar9_CMatPP& uf, vec_ar9_RMatPP& ut) {

	CMatXN9 &tempF = KernerFFTW_N9::getC2R_CMat();
	makeFlat<vec_ar9_CMatPP, CMatXN9>(uf, tempF);
	KernerFFTW_N9::computeR2C();
	RMatXN9 &tempT = KernerFFTW_N9::getC2R_RMat();
	makeStruct<vec_ar9_RMatPP, RMatXN9>(ut, tempT);
	
}