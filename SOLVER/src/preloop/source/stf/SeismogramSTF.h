// SeismogramSTF.h
// Used for off-axis sources 

#pragma once 

#include "STF.h"
#include "eigenc.h"
#include "eigenp.h"

class SeismogramSTF : public STF {
public:
	SeismogramSTF(const RMatX3 trace, double dt_fwd, double duration_fwd, double hdur_fwd, double decay_fwd, 
		const RDMatXX &adjoint_params, const IColX &filter_types);
	std::string verbose() const;
	
private:
    double mHalfDuration;
    double mDecay;
	
};