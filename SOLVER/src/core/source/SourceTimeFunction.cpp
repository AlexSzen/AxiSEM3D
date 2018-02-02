// SourceTimeFunction.cpp
// created by Kuangdai on 7-Apr-2016 
// source time function

#include "SourceTimeFunction.h"

SourceTimeFunction::SourceTimeFunction(const std::vector<Real> &stfs, const std::vector<Real> &stfp, const std::vector<Real> &stfz, Real dt, Real shift):
mSTFs(stfs), mSTFp(stfp), mSTFz(stfz), mDeltaT(dt), mShift(shift) {
    // nothing
}
