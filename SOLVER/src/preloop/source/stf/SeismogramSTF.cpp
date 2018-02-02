// SeismogramSTF.cpp
// STF produced by seismogram 
#include "SeismogramSTF.h"
#include <cmath>
#include <string>
#include <sstream>

SeismogramSTF::SeismogramSTF(RMatX3 trace, int downsampling, double dt_fwd, double length_fwd, double hdur_fwd, double decay_fwd):
mDownsampling(downsampling) {
	
	int length_downsampled = trace.rows() / downsampling + 1; // +1 from time 0
	RMatX3 trace_downsampled(length_downsampled, 3);
	
	int it_downsampled = 0;
	for (int it = 0; it < trace.rows(); it++) {
		
		if ( it % downsampling == 0) {
			mSTFs.push_back(trace(it, 0));
			mSTFp.push_back(trace(it, 1));
			mSTFz.push_back(trace(it, 2));
			it_downsampled ++;
		}	
	}
	
	
	
}