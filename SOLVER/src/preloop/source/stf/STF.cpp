// STF.cpp
// created by Kuangdai on 9-May-2016
// source time function

#include "STF.h"
#include "SourceTimeFunction.h"
#include "Domain.h"

void STF::release(Domain &domain) const {
    std::vector<Real> tss(mSTFs.begin(), mSTFs.end());
	std::vector<Real> tsp(mSTFp.begin(), mSTFp.end());
	std::vector<Real> tsz(mSTFz.begin(), mSTFz.end());

    domain.addSTF(new SourceTimeFunction(tss, tsp, tsz, mDeltaT, mShift));
}
