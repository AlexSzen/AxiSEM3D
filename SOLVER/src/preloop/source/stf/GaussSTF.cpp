// GaussSTF.cpp
// created by Alex on 8-May-2016
// gaussian stf

#include "GaussSTF.h"
#include <cmath>
#include <string>
#include <sstream>

GaussSTF::GaussSTF(double dt, double duration, double hdur, double decay):
mHalfDuration(hdur), mDecay(decay) {
    mDeltaT = dt;
    int nStepBeforeZero = ceil(1.5 * mHalfDuration / mDeltaT);
    int nStepAfterZero = ceil(duration / mDeltaT);
    mShift = nStepBeforeZero * mDeltaT;
    int nStep = nStepBeforeZero + nStepAfterZero;
    for (int i = 0; i <= nStep; i++) {
        double t = -mShift + i * mDeltaT;
        mSTFs.push_back(exp(-pow((mDecay / mHalfDuration * t), 2.)) * mDecay / (mHalfDuration * sqrt(pi)));
		mSTFp.push_back(exp(-pow((mDecay / mHalfDuration * t), 2.)) * mDecay / (mHalfDuration * sqrt(pi)));
		mSTFz.push_back(exp(-pow((mDecay / mHalfDuration * t), 2.)) * mDecay / (mHalfDuration * sqrt(pi)));
   
    }
}

std::string GaussSTF::verbose() const {
    std::stringstream ss;
    ss << "\n=================== Source Time Function ===================" << std::endl;
    ss << "  Time Step               =   " << mDeltaT << std::endl;
    ss << "  Number of Steps         =   " << mSTFs.size() << std::endl;
    ss << "  Total Duration          =   " << mDeltaT * mSTFs.size() << std::endl;
    ss << "  Duration after Origin   =   " << mDeltaT * mSTFs.size() - mShift << std::endl;
    ss << "  Shift before Origin     =   " << mShift << std::endl;
    ss << "  Time Series Type        =   Gaussian" << std::endl;
    ss << "  Half Duration           =   " << mHalfDuration << std::endl;
    ss << "  Decay Factor            =   " << mDecay << std::endl;
    ss << "=================== Source Time Function ===================\n" << std::endl;
    return ss.str();
}
