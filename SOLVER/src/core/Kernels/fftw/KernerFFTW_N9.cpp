// fftw for kernels 
// we fft after going fourier->physical

#include "KernerFFTW_N9.h"

int KernerFFTW_N9::sNmax = 0; 
int KernerFFTW_N9::sTotStepsTime = 0;
int KernerFFTW_N9::sTotStepsFreq = 0;
PlanFFTW KernerFFTW_N9::sR2CPlan;
PlanFFTW KernerFFTW_N9::sC2RPlan;
RMatXN9 KernerFFTW_N9::sR2C_RMat;
CMatXN9 KernerFFTW_N9::sR2C_CMat;
RMatXN9 KernerFFTW_N9::sC2R_RMat;
CMatXN9 KernerFFTW_N9::sC2R_CMat;

void KernerFFTW_N9::initialize(int totSteps) {
	
	sTotStepsTime = totSteps;
	sTotStepsFreq = totSteps/2 +1;
	
	int NT = sTotStepsTime;
	int n[] = {NT};
		
	int xx = nPntElem * 9;
	
	sR2C_RMat = RMatXN9(NT, xx);
	sR2C_CMat = CMatXN9(NT, xx);
	sC2R_RMat = RMatXN9(NT, xx);
	sC2R_CMat = CMatXN9(NT, xx);
	Real *r2c_r = &(sR2C_RMat(0, 0));
    Complex *r2c_c = &(sR2C_CMat(0, 0));
	sR2CPlan = planR2CFFTW(1, n, xx, r2c_r, n, 1, sTotStepsTime, complexFFTW(r2c_c), n, 1, sTotStepsFreq, FFTW_PATIENT);   
	Real *c2r_r = &(sC2R_RMat(0, 0));
	Complex *c2r_c = &(sC2R_CMat(0, 0));
	sC2RPlan = planC2RFFTW(1, n, xx, complexFFTW(c2r_c), n, 1, sTotStepsFreq, c2r_r, n, 1, sTotStepsTime, FFTW_PATIENT); 
		

	
}

void KernerFFTW_N9::finalize() {
	
    destroyFFTW(sR2CPlan);
    destroyFFTW(sC2RPlan);
    
}

void KernerFFTW_N9::computeR2C() {

    execFFTW(sR2CPlan);
    Real inv_time = one / (Real)sTotStepsTime;
    sR2C_CMat *= inv_time;
}

void KernerFFTW_N9::computeC2R() {

    execFFTW(sC2RPlan);
}



