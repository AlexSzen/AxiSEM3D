// STFCollection.h
// source time function collection 
// we now have a collection because of off axis sources 
#pragma once

#include <vector>
#include <string>
#include "global.h"
#include "STF.h"

class Domain;
class Parameters;

class STFCollection {
public:
	STFCollection(double hdur, double duration, std::string mstf, double dt, int enforceMaxSteps, std::string offaxis_file, bool kernels);
    ~STFCollection();

    void release(Domain &domain) const;

    std::string verbose();

    static void buildInparam(STFCollection *&stf, const Parameters &par, double dt, int verbose);
    
    int getSize() const {return mSTFs[0]->getSize();};

protected:
    double mDeltaT;
    double mShift;
    std::vector<STF *> mSTFs;
};
