// DomainIO.cpp
// IO for domain wise wavefields 

#include "DomainIO.h"
#include "Parameters.h"
#include "NetCDF_Writer.h"
#include "XMPI.h"
#include "iostream"

void DomainIO::initialize(int totalRecordSteps, int recordInterval, int bufferSize, int totNuProc,
	double srcLat, double srcLon, double srcDep) {
	
	// source location
	mSrcLat = srcLat;
	mSrcLon = srcLon;
	mSrcDep = srcDep;
	
	mNetCDF = new NetCDF_Writer();
	std::string fname = Parameters::sOutputDirectory + "/wavefields/wavefield_db.nc4";
	
	
	// fill start and counts
	int temp_startElemNu; 
	std::vector<int> temp_countElemNu(XMPI::nproc(),0);

	XMPI::gather(totNuProc, temp_countElemNu, true); //fourier fields 
	temp_startElemNu = std::accumulate(temp_countElemNu.begin(), temp_countElemNu.begin()+XMPI::rank(),0);
	
	int totNu = XMPI::sum(totNuProc);
	
	mStartTime.push_back(0);
	mCountTime.push_back(bufferSize); 

	mStartWvf.push_back(0);
	mStartWvf.push_back(temp_startElemNu);
	mStartWvf.push_back(0);
	mStartWvf.push_back(0);
	mStartWvf.push_back(0);
	
	mCountWvf.push_back(1); // AT THE MOMENT CANNOT DUMP BUFFER SIZE BECAUSE OF NETCDF PROBLEM, SEE ARRAY TEST.
	mCountWvf.push_back(temp_countElemNu[XMPI::rank()]);
	mCountWvf.push_back(6);
	mCountWvf.push_back(nPntEdge);
	mCountWvf.push_back(nPntEdge);

	if (XMPI::root()) {
		
		//dimensions 
		std::vector<size_t> dimsTime;
		std::vector<size_t> dimsWvf;
		dimsTime.push_back(totalRecordSteps);
		dimsWvf.push_back(totalRecordSteps);
		dimsWvf.push_back(totNu);
		dimsWvf.push_back(6);
		dimsWvf.push_back(nPntEdge);
		dimsWvf.push_back(nPntEdge);
		
		mNetCDF->open(fname, false);
		mNetCDF->defModeOn();
		
		mNetCDF->defineVariable<Real>("time", dimsTime);
		mNetCDF->defineVariable<Real>("displacement_wavefield", dimsWvf);
		
		//number of cores 
		mNetCDF->addAttribute("", "number_procs", XMPI::nproc());
		// source location
		mNetCDF->addAttribute("", "source_latitude", mSrcLat);
		mNetCDF->addAttribute("", "source_longitude", mSrcLon);
		mNetCDF->addAttribute("", "source_depth", mSrcDep);
		mNetCDF->flush();
		
		mNetCDF->close();
	}
	mNetCDF->openParallel(fname);
}

void DomainIO::finalize() {
	
	mNetCDF->close();
	delete mNetCDF;
}

void DomainIO::dumpToFile(const vec_vec_ar6_RMatPP& bufferDisp, const RColX& bufferTime, int bufferLineTime) {
	
	mCountTime[0] = bufferLineTime;
	
	mNetCDF->writeVariableChunk("time", bufferTime, mStartTime, mCountTime);
	
	for (int it = 0; it < bufferLineTime; it++) { //AGAIN, NEED TO INVESTIGATE THIS THING
		mNetCDF->writeVariableChunk("displacement_wavefield", bufferDisp[it], mStartWvf, mCountWvf);
        mStartWvf[0]++;
    }
	
	mNetCDF->flush();

	mStartTime[0] += bufferLineTime;
	
}

