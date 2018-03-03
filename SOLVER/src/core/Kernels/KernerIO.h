// KernerIO.h

#pragma once 

#include "eigenc.h"

class NetCDF_Reader;
class NetCDF_Writer;

class KernerIO {
	
public:
	
	KernerIO(bool dumpTimeKernels, int temp_startElem, int temp_countElem, int totSteps, int bufferSize);	
	void initialize(int totNu, int startElemNu, int countElemNu, int totElem, const std::vector<int> &nusKer,const std::vector<int> &nrsKer);
	void finalize();
	void dumpToFile(const vec_vec_ar6_RMatPP &kernels, int bufferSize);
	
	void loadNus(std::vector<int> &Nus);
	void loadNrs(std::vector<int> &Nrs); 
	void loadWavefield(vec_vec_ar6_RMatPP &disp, int bufferSize);
	void loadMaterial(vec_ar12_RMatPP &materials);
private:
	
	NetCDF_Writer *mNetCDF_w; 
	NetCDF_Reader *mNetCDF_r; //used to load fwd wvf 
	
	size_t mStartElem, mCountElem, mStartElemNuFwd, mCountElemNuFwd, mStartElemNuKernels, mCountElemNuKernels; //storing these to avoid redundant computations
	int mTotSteps; 
	int mBufferSize;
	int mNumLoads; //load forward wavefield in chunks
	bool mDumpTimeKernels;
	
	
};