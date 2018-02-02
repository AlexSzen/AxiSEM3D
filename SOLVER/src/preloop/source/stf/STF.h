// STF.h
// created by Kuangdai on 9-May-2016
// source time function

#pragma once

#include <vector>
#include <string>
#include "global.h"

class Domain;
class Parameters;


class STF {
	
friend class STFCollection; 

public:
    virtual ~STF() {};

    void release(Domain &domain) const;

    virtual std::string verbose() const = 0;
    
    int getSize() const {return mSTFs.size();};

protected:
    double mDeltaT;
    double mShift;
    std::vector<double> mSTFs;
	std::vector<double> mSTFp;
	std::vector<double> mSTFz;

};
