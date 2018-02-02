// SourceTimeFunction.h
// created by Kuangdai on 7-Apr-2016 
// source time function


#pragma once

#include <vector>
#include "global.h"

class SourceTimeFunction {
public:
    SourceTimeFunction(const std::vector<Real> &stfs, const std::vector<Real> &stfp, const std::vector<Real> &stfz, Real dt, Real shift);
    
    int getSize() const {return mSTFs.size();};
    Real getFactorS(int tstep) const {return mSTFs[tstep];};
	Real getFactorP(int tstep) const {return mSTFp[tstep];};
	Real getFactorZ(int tstep) const {return mSTFz[tstep];};

    Real getDeltaT() const {return mDeltaT;};
    Real getShift() const {return mShift;};
    
private:
    std::vector<Real> mSTFs;
	std::vector<Real> mSTFp;
	std::vector<Real> mSTFz;

    Real mDeltaT;
    Real mShift;
};