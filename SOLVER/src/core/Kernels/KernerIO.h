// KernerIO.h

#pragma once 

#include "eigenc.h"

class NetCDF_Reader;
class NetCDF_Writer;

class KernerIO {
	
public:
		
	void initialize(int dumpTimeKernels, int numFilters, int temp_startElem, int temp_countElem, int totSteps);
	void finalize();
	void dumpToFile(vec_vec_ar12_RMatPP &kernels, int numFilters);
	
	void loadNus(std::vector<int> &Nus);
	void loadWavefield(vec_vec_ar6_RMatPP &disp, std::vector<int> &Nus);
	void loadMaterial(vec_ar12_RMatPP &materials,  std::vector<int> &Nus);
private:
	
	NetCDF_Writer *mNetCDF_w; 
	NetCDF_Reader *mNetCDF_r; //used to load fwd wvf 
	
	size_t mStartElem, mCountElem, mStartElemNu, mCountElemNu; //storing these to avoid redundant computations
	int mTotSteps; 
};