// PreloopFFTW_time.h
// perform FFT using fftw
// used for creating source time functions from seismograms. 
// although in preloop STF is double, we use Real since it is being cast to Real when it's realeased anyway 
 
#pragma once

#include <fftw3.h>
#include <vector>
#include "eigenc.h"

#ifdef _USE_DOUBLE
    typedef fftw_plan PlanFFTW;
    #define complexFFTW reinterpret_cast<fftw_complex*>
    #define planR2CFFTW fftw_plan_many_dft_r2c
    #define planC2RFFTW fftw_plan_many_dft_c2r
    #define destroyFFTW fftw_destroy_plan
    #define execFFTW fftw_execute
#else
    typedef fftwf_plan PlanFFTW;
    #define complexFFTW reinterpret_cast<fftwf_complex*>
    #define planR2CFFTW fftwf_plan_many_dft_r2c
    #define planC2RFFTW fftwf_plan_many_dft_c2r
    #define destroyFFTW fftwf_destroy_plan
    #define execFFTW fftwf_execute
#endif


class PreloopFFTW_time {
public:
    // initialize plans
    static void initialize(int totSteps);
    // finalize plans
    static void finalize();
    
	static RMatX3 &getR2C_RMat() { return sR2C_RMat;};
    static CMatX3 &getR2C_CMat() { return sR2C_CMat;};
    static RMatX3 &getC2R_RMat() { return sC2R_RMat;};    
    static CMatX3 &getC2R_CMat() { return sC2R_CMat;};
     
    // forward, real => complex
    static void computeR2C();
    // backward, complex => real
    static void computeC2R();
    

    
private:
    static int sNmax;
	static int sTotStepsTime;
	static int sTotStepsFreq;
    static PlanFFTW sR2CPlan;
    static PlanFFTW sC2RPlan;
    static RMatX3 sR2C_RMat;
    static CMatX3 sR2C_CMat;
    static RMatX3 sC2R_RMat;
    static CMatX3 sC2R_CMat;

};