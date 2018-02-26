// SeismogramSTF.cpp
// STF produced by seismogram 
#include "SeismogramSTF.h"
#include "Processor.h"
#include <cmath>
#include <string>
#include <sstream>

SeismogramSTF::SeismogramSTF(const RMatX3 trace, double dt_fwd, double duration_fwd, double hdur_fwd, double decay_fwd, 
	const RDMatXX &adjoint_params, const IColX &filter_types):
mHalfDuration(hdur_fwd), mDecay(decay_fwd) {
	
	mDeltaT = dt_fwd;
    int nStepBeforeZero = ceil(1.5 * mHalfDuration / mDeltaT);
    int nStepAfterZero = ceil(duration_fwd / mDeltaT);
    mShift = nStepBeforeZero * mDeltaT;
    int nStep = nStepBeforeZero + nStepAfterZero;
	
	if (nStep + 1 != trace.rows()) throw std::runtime_error("SeismogramSTF::SeismogramSTF : length of seismogram and STF inconsistent.");
	
	mSTFs.assign(nStep + 1, (double) 0.);
	mSTFp.assign(nStep + 1, (double) 0.);
	mSTFz.assign(nStep + 1, (double) 0.);
	
	// apply each measurement and add to the STF 
	for (int i_measurement = 0; i_measurement < adjoint_params.rows(); i_measurement++) {
		
		RDCol2 window = adjoint_params.block(i_measurement, 0, 1, 2);		
		
		RMatX3 trace_measurement_T = trace;
		CMatX3 trace_measurement_F;
		
		Processor::transformT2F(trace_measurement_T, trace_measurement_F);
		Processor::filter(trace_measurement_F, filter_types[i_measurement]);
		Processor::transformF2T(trace_measurement_F, trace_measurement_T);
		Processor::taper(trace_measurement_T, (Real) window(0), (Real) window(1));

		for (int i = 0; i <= nStep; i++) { //time reversed seismogram 
	   	 mSTFs[i] += 1e15 * trace_measurement_T(nStep - i,0);
	   	 mSTFp[i] += 1e15 * trace_measurement_T(nStep - i,1);
	   	 mSTFz[i] += 1e15 * trace_measurement_T(nStep - i,2);
	    }
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