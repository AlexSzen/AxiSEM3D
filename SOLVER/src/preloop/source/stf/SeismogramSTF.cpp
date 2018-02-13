// SeismogramSTF.cpp
// STF produced by seismogram 
#include "SeismogramSTF.h"
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
SeismogramSTF::SeismogramSTF(const RMatX3 trace, double dt_fwd, double duration_fwd, double hdur_fwd, double decay_fwd):
mHalfDuration(hdur_fwd), mDecay(decay_fwd) {
	
	mDeltaT = dt_fwd;
    int nStepBeforeZero = ceil(1.5 * mHalfDuration / mDeltaT);
    int nStepAfterZero = ceil(duration_fwd / mDeltaT);
    mShift = nStepBeforeZero * mDeltaT;
    int nStep = nStepBeforeZero + nStepAfterZero;
	
	if (nStep + 1 != trace.rows()) throw std::runtime_error("SeismogramSTF::SeismogramSTF : length of seismogram and STF inconsistent.");
	
	for (int i = 0; i <= nStep; i++) {
   	 mSTFs.push_back(trace(i,0));
   	 mSTFp.push_back(trace(i,1));
   	 mSTFz.push_back(trace(i,2));
    }
	
}

std::string SeismogramSTF::verbose() const {
	std::stringstream ss;
	ss << "\n=================== Source Time Function ===================" << std::endl;
	ss << "  Time Step               =   " << mDeltaT << std::endl;
	ss << "  Number of Steps         =   " << mSTFs.size() << std::endl;
	ss << "  Total Duration          =   " << mDeltaT * mSTFs.size() << std::endl;
	ss << "  Duration after Origin   =   " << mDeltaT * mSTFs.size() - mShift << std::endl;
	ss << "  Shift before Origin     =   " << mShift << std::endl;
	ss << "  Time Series Type        =   Seismogram" << std::endl;
	ss << "  Half Duration           =   N/A"  << std::endl;
	ss << "  Decay Factor            =   N/A" << std::endl;
	ss << "=================== Source Time Function ===================\n" << std::endl;
	return ss.str();
	
}