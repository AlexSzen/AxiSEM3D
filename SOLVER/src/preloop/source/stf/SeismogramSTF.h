// SeismogramSTF.h
// Used for off-axis sources 

#pragma once 

#include "STF.h"

class SeismogramSTF : public STF {

public:
	SeismogramSTF(RMatX3 trace, int downsampling, double dt_fwd, double length_fwd, double hdur_fwd, double decay_fwd);
	std::string verbose();
	
private:
	
	int mDownsampling;
	
};