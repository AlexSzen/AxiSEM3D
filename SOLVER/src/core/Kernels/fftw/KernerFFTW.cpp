// fftw for kernels 
// we fft after going fourier->physical

#include "KernerFFTW.h"

static int KernerFFTW::sNmax = 0; 
static int KernerFFTW::sTotStepsTime = 0;
static std::vector<fftw_plan> KernerFFTW::sR2CPlans;
static std::vector<fftw_plan> KernerFFTW::sC2RPlans;
static std::vector<RMatXX> KernerFFTW::sR2C_RMats;
static std::vector<CMatXX> KernerFFTW::sR2C_CMats;
static std::vector<RMatXX> KernerFFTW::sC2R_RMats;
static std::vector<CMatXX> KernerFFTW::sC2R_CMats;

void KernerFFTW::initialize(int Nmax, int totSteps) {
	
	sTotStepsTime = totSteps;
	sTotStepsFreq = totSteps/2 +1;
	sNmax = Nmax;
	sR2C_RMats.reserve(sNmax);
	sR2C_CMats.reserve(sNmax);
	sC2R_RMats.reserve(sNmax);
	sC2R_CMats.reserve(sNmax);
	
	int NT = mTotSteps;
	int n[] = {NT};
	
	for (int NR = 1; NR <= Nmax; NR++) {
		
		xx = nPntElem * NR;
		
		sR2C_RMats.push_back(RMatXX(NT, xx));
		sR2C_CMats.push_back(CMatXX(NT, xx));
		sC2R_RMats.push_back(RMatXX(NT, xx));
		sC2R_CMats.push_back(CMatXX(NT, xx));
		Real *t2f_r = &(sR2C_RMats[NR - 1](0, 0));
        Complex *r2c_c = &(sR2C_CMats[NR - 1](0, 0));
		sR2CPlans.push_back(fftw_plan_many_dft_r2c(
			1, n, xx, r2c_r, n, 1, sTotStepsTime, reinterpret_cast<fftw_complex*>(r2c_c), n, 1, sTotStepsFreq, FFTW_PATIENT));   
		Real *c2r_r = &(sC2R_RMats[NR - 1](0, 0));
		Complex *c2r_c = &(sC2R_CMats[NR - 1](0, 0));
		sC2RPlans.push_back(fftw_plan_many_dft_c2r(
			1, n, xx, reinterpret_cast<fftw_complex*>(c2r_c), n, 1, sTotStepsFreq, c2r_r, n, 1, sTotStepsTime, FFTW_PATIENT)); 
			
	}
	
}

void KernerFFTW::finalize() {
    for (int i = 0; i < sNmax; i++) {
        fftw_destroy_plan(sR2CPlans[i]);
        fftw_destroy_plan(sC2RPlans[i]);
    }
    sNmax = 0;
}

void KernerFFTW::computeR2C(int nr) {

    fftw_execute(sR2CPlans[nr - 1]);
    Real inv_time = one / (Real)sTotStepsTime;
    sR2C_CMats[nr - 1] *= inv_time;
}

void KernerFFTW::computeC2R(int nr) {

    fftw_execute(sC2RPlans[nr - 1]);
}



