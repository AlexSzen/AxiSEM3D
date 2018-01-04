// fftw for kernels 
// we fft after going fourier->physical

#include "KernerFFTW_N6.h"

int KernerFFTW_N6::sNmax = 0; 
int KernerFFTW_N6::sTotStepsTime = 0;
int KernerFFTW_N6::sTotStepsFreq = 0;
PlanFFTW KernerFFTW_N6::sR2CPlan;
PlanFFTW KernerFFTW_N6::sC2RPlan;
RMatXN6 KernerFFTW_N6::sR2C_RMat;
CMatXN6 KernerFFTW_N6::sR2C_CMat;
RMatXN6 KernerFFTW_N6::sC2R_RMat;
CMatXN6 KernerFFTW_N6::sC2R_CMat;

void KernerFFTW_N6::initialize(int totSteps) {
	
	int ndim = 6;
	
	sTotStepsTime = totSteps;
	sTotStepsFreq = totSteps/2 + 1;
	
	int NT = sTotStepsTime;
	int NF = sTotStepsFreq;
	int n[] = {NT};
		
	int xx = nPntElem * ndim;
	
	sR2C_RMat = RMatXN6(NT, xx);
	sR2C_CMat = CMatXN6(NF, xx);
	sC2R_RMat = RMatXN6(NT, xx);
	sC2R_CMat = CMatXN6(NF, xx);
	Real *r2c_r = &(sR2C_RMat(0, 0));
    Complex *r2c_c = &(sR2C_CMat(0, 0));
	sR2CPlan = planR2CFFTW(1, n, xx, r2c_r, n, 1, NT, complexFFTW(r2c_c), n, 1, NF, FFTW_PATIENT);   
	Real *c2r_r = &(sC2R_RMat(0, 0));
	Complex *c2r_c = &(sC2R_CMat(0, 0));
	sC2RPlan = planC2RFFTW(1, n, xx, complexFFTW(c2r_c), n, 1, NF, c2r_r, n, 1, NT, FFTW_PATIENT); 
		

	
}

void KernerFFTW_N6::finalize() {
	
    destroyFFTW(sR2CPlan);
    destroyFFTW(sC2RPlan);
    
}

void KernerFFTW_N6::computeR2C() {

    execFFTW(sR2CPlan);
    Real inv_time = one / (Real)sTotStepsTime;
    sR2C_CMat *= inv_time;
}

void KernerFFTW_N6::computeC2R() {

    execFFTW(sC2RPlan);
}



