

#pragma once 

#include "KernerFFTW.h"


class KernerFFTW_N9 {
public:
    // initialize plans
    static void initialize(int totSteps);
    // finalize plans
    static void finalize();
    
	static RMatXN9 &getR2C_RMat() { return sR2C_RMat;};
    static CMatXN9 &getR2C_CMat() { return sR2C_CMat;};
    static RMatXN9 &getC2R_RMat() { return sC2R_RMat;};    
    static CMatXN9 &getC2R_CMat() { return sC2R_CMat;};
     
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
    static RMatXN9 sR2C_RMat;
    static CMatXN9 sR2C_CMat;
    static RMatXN9 sC2R_RMat;
    static CMatXN9 sC2R_CMat;

};