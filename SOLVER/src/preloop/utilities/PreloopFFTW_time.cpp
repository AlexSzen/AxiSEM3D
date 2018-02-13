// fftw for kernels 
// we fft after going fourier->physical

#include "PreloopFFTW_time.h"

int PreloopFFTW_time::sNmax = 0; 
int PreloopFFTW_time::sTotStepsTime = 0;
int PreloopFFTW_time::sTotStepsFreq = 0;
PlanFFTW PreloopFFTW_time::sR2CPlan;
PlanFFTW PreloopFFTW_time::sC2RPlan;
RMatX3 PreloopFFTW_time::sR2C_RMat;
CMatX3 PreloopFFTW_time::sR2C_CMat;
RMatX3 PreloopFFTW_time::sC2R_RMat;
CMatX3 PreloopFFTW_time::sC2R_CMat;

void PreloopFFTW_time::initialize(int totSteps) {
	
	
	sTotStepsTime = totSteps;
	sTotStepsFreq = totSteps/2 +1;
	
	int NT = sTotStepsTime;
	int NF = sTotStepsFreq;
	int n[] = {NT};
		
	int xx = 3;
	
	sR2C_RMat = RMatX3(NT, xx);
	sR2C_CMat = CMatX3(NF, xx);
	sC2R_RMat = RMatX3(NT, xx);
	sC2R_CMat = CMatX3(NF, xx);
	Real *r2c_r = &(sR2C_RMat(0, 0));
    Complex *r2c_c = &(sR2C_CMat(0, 0));
	sR2CPlan = planR2CFFTW(1, n, xx, r2c_r, n, 1, NT, complexFFTW(r2c_c), n, 1, NF, FFTW_PATIENT);   
	Real *c2r_r = &(sC2R_RMat(0, 0));
	Complex *c2r_c = &(sC2R_CMat(0, 0));
	sC2RPlan = planC2RFFTW(1, n, xx, complexFFTW(c2r_c), n, 1, NF, c2r_r, n, 1, NT, FFTW_PATIENT); 
		

	
}

void PreloopFFTW_time::finalize() {
	
    destroyFFTW(sR2CPlan);
    destroyFFTW(sC2RPlan);
    
}

void PreloopFFTW_time::computeR2C() {

    execFFTW(sR2CPlan);
    Real inv_time = one / (Real)sTotStepsTime;
    sR2C_CMat *= inv_time;
}

void PreloopFFTW_time::computeC2R() {

    execFFTW(sC2RPlan);
}



