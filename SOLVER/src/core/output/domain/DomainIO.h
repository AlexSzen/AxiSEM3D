//Alex
// IO for domain wide wavefields 

#pragma once 

#include "eigenc.h"

class NetCDF_Writer;

class DomainIO {
	
public:
	
	void initialize(int totalRecordSteps, int recordInterval, int bufferSize, int totNuProc,
		double srcLat, double srcLon, double srcDep);
	void finalize();
	
	void dumpToFile(const vec_vec_ar6_RMatPP& bufferDisp, const RColX& bufferTime, int bufferLineTime);
	
private:
	
	// writer 
	NetCDF_Writer *mNetCDF = 0;
	
	// location in time  
	int mCurrentRow = 0;
	
	
	std::vector<size_t> mStartWvf, mCountWvf, mStartTime, mCountTime;
	
	// source location
	double mSrcLat, mSrcLon, mSrcDep;
	
};