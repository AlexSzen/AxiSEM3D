// KernerIO.h

#pragma once 

#include "eigenc.h"

class NetCDF_Reader;
class NetCDF_Writer;

class KernerIO {
	
public:
		
	void initialize(int dumpTimeKernels, int numKernels);
	void finalize();
	void dumpToFile();
	
	void loadNus(std::vector<int> &Nus, int temp_startElem, int temp_countElem);
	void loadWavefield(vec_vec_ar6_RMatPP &disp, std::vector<int> &Nus, int &totSteps, int temp_startElem, int temp_countElem);
	void loadMaterial(vec_ar12_RMatPP &materials,  std::vector<int> &Nus,int temp_startElem, int temp_countElem);
private:
	
	NetCDF_Writer *mNetCDF_w; 
	NetCDF_Reader *mNetCDF_r; //used to load fwd wvf 
};