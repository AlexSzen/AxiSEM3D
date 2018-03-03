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
		double measurement = adjoint_params(i_measurement, 2);
		
		RMatX3 trace_measurement_T = trace;
		CMatX3 trace_measurement_F;
		
		Processor::transformT2F(trace_measurement_T, trace_measurement_F);
		Processor::derivate(trace_measurement_F);
		Processor::filter(trace_measurement_F, filter_types[i_measurement]);
		Processor::transformF2T(trace_measurement_F, trace_measurement_T);
		Processor::taper(trace_measurement_T, (Real) window(0), (Real) window(1));
		
		Real norm_s = 0.;
		Real norm_p = 0.;
		Real norm_z = 0.;
		
		for (int it = 0; it <= nStep; it++) { //normalization factor : for traveltime tomo it's time integrated squared velocity
			norm_s += mDeltaT * trace_measurement_T(it, 0) * trace_measurement_T(it, 0);
			norm_p += mDeltaT * trace_measurement_T(it, 1) * trace_measurement_T(it, 1);
			norm_z += mDeltaT * trace_measurement_T(it, 2) * trace_measurement_T(it, 2);
		}
		
		if (norm_s == 0.) norm_s = 1.;
		if (norm_p == 0.) norm_p = 1.;
		if (norm_z == 0.) norm_z = 1.;
		
		for (int i = 0; i <= nStep; i++) { //time reversed seismogram. 
		 
	   	 mSTFs[i] += measurement * trace_measurement_T(nStep - i, 0) / norm_s;
	   	 mSTFp[i] += measurement * trace_measurement_T(nStep - i, 1) / norm_p;
	   	 mSTFz[i] += measurement * trace_measurement_T(nStep - i, 2) / norm_z;
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