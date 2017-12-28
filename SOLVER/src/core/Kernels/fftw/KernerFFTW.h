//time to freq fftw for kernels 
// it is performed once we are already in Nr space.
// DO I NEED TO SWITCH TO RMatXX_RM???

#include <fftw3.h>
#include <vector>
#include "eigenp.h"

class KernerFFTW {
public:
    // initialize plans
    static void initialize(int Nmax, int totSteps);
    // finalize plans
    static void finalize();
    
	static RMatXX &getR2C_RMat(int nr) { return sR2C_RMats[nr - 1];};
    static CMatXX &getR2C_CMat(int nr) { return sR2C_CMats[nr - 1];};
    static RMatXX &getC2R_RMat(int nr) { return sC2R_RMats[nr - 1];};    
    static CMatXX &getC2R_CMat(int nr) { return sC2R_CMats[nr - 1];};
     
    // forward, real => complex
    static void computeR2C(int nr);
    // backward, complex => real
    static void computeC2R(int nr);
    

    
private:
    static int sNmax;
	static int sTotStepsTime;
	static int sTotStepsFreq;
    static std::vector<fftw_plan> sR2CPlans;
    static std::vector<fftw_plan> sC2RPlans;
    static std::vector<RMatXX> sR2C_RMats;
    static std::vector<CMatXX> sR2C_CMats;
    static std::vector<RMatXX> sC2R_RMats;
    static std::vector<CMatXX> sC2R_CMats;

};
