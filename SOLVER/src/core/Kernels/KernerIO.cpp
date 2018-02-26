//KernerIO.cpp

#include "KernerIO.h"
#include "NetCDF_Writer.h"
#include "NetCDF_Reader.h"
#include "Parameters.h"
#include "XMPI.h"
#include <iostream>

KernerIO::KernerIO(bool dumpTimeKernels, int temp_startElem, int temp_countElem, int totSteps):
mDumpTimeKernels(dumpTimeKernels), mStartElem(temp_startElem), mCountElem(temp_countElem), mTotSteps(totSteps) 
 {
	
	mNetCDF_w = new NetCDF_Writer();
	mNetCDF_r = new NetCDF_Reader();
	
}
void KernerIO::initialize(int totNu, int startElemNu, int countElemNu, int totElem) {


	std::string fname_ker = Parameters::sOutputDirectory + "/kernels/kernels_db.nc4";
	mStartElemNuKernels = startElemNu;
	mCountElemNuKernels = countElemNu;
	
	if (XMPI::root()) {
	
		//create dims and vars in kernel file 
		mNetCDF_w->open(fname_ker, true);
		
		std::vector<size_t> dims_kernels; 
		std::vector<size_t> dims_nus;
		
		if (mDumpTimeKernels) {
			dims_kernels.push_back( mTotSteps );
		} else {
			dims_kernels.push_back( 1 );
		}		
		dims_kernels.push_back( totNu );
		dims_kernels.push_back( 6 ); //real and imag parts of elastic radial ani kernels. 
		dims_kernels.push_back( nPntEdge );
		dims_kernels.push_back( nPntEdge );
		
		dims_nus.push_back(totElem);
		
		mNetCDF_w->defModeOn();
		mNetCDF_w->defineVariable<Real>("Kernels", dims_kernels);
		mNetCDF_w->defineVariable<int>("Nus", dims_nus);
		mNetCDF_w->defineVariable<int>("Nrs", dims_nus);

		
		mNetCDF_w->defModeOff();
		mNetCDF_w->close();
		
	}

	
	
}

void KernerIO::finalize() {
	

	delete mNetCDF_w;
	delete mNetCDF_r;
	
}

void KernerIO::dumpToFile(const vec_vec_ar6_RMatPP &kernels, const std::vector<int> &nusKer,const std::vector<int> &nrsKer) {
	
	std::vector<size_t> startKernels, countKernels, startNus, countNus;
	startKernels.push_back(0);
	startKernels.push_back(mStartElemNuKernels);
	startKernels.push_back(0);
	startKernels.push_back(0);
	startKernels.push_back(0);
	
	countKernels.push_back(1); //have to write one by one because of netcdf issue
	countKernels.push_back(mCountElemNuKernels);
	countKernels.push_back(6);
	countKernels.push_back(nPntEdge);
	countKernels.push_back(nPntEdge); 
	
	startNus.push_back(mStartElem);
	countNus.push_back(mCountElem);
	
	std::string fname = Parameters::sOutputDirectory + "/kernels/kernels_db.nc4";
	mNetCDF_w->openParallel(fname);
	
	int dumpSize = 1;
	if (mDumpTimeKernels) {
		dumpSize = kernels.size();
	}
	for (int it = 0; it < dumpSize ; it ++) {
		
		mNetCDF_w->writeVariableChunk("Kernels", kernels[it], startKernels, countKernels);
		startKernels[0]++;
		
	}
	
	mNetCDF_w->writeVariableChunk("Nus", nusKer, startNus, countNus);
	mNetCDF_w->writeVariableChunk("Nrs", nrsKer, startNus, countNus);

	mNetCDF_w->close();
	
	
	
}

void KernerIO::loadNus(std::vector<int> &Nus) {
	
	std::string fname_wvf = Parameters::sOutputDirectory + "/wavefields/wavefield_db_fwd.nc4";
	mNetCDF_r->openParallel(fname_wvf);

	//read nus 
	std::vector<size_t> startElem, countElem;
	startElem.push_back(mStartElem);
	countElem.push_back(mCountElem);
	Nus.assign(mCountElem, 0);
	
	mNetCDF_r->readVariableChunk("Nus", Nus, startElem, countElem);
	
	mNetCDF_r->close();
	
}

void KernerIO::loadNrs(std::vector<int> &Nrs) {
	
	std::string fname_wvf = Parameters::sOutputDirectory + "/wavefields/wavefield_db_fwd.nc4";
	mNetCDF_r->openParallel(fname_wvf);

	//read nus 
	std::vector<size_t> startElem, countElem;
	startElem.push_back(mStartElem);
	countElem.push_back(mCountElem);
	Nrs.assign(mCountElem, 0);
	
	mNetCDF_r->readVariableChunk("Nrs", Nrs, startElem, countElem);
	
	mNetCDF_r->close();
	
}

void KernerIO::loadWavefield(vec_vec_ar6_RMatPP &disp, std::vector<int> &Nus, std::vector<int> &Nrs) {
	
	std::string fname_wvf = Parameters::sOutputDirectory + "/wavefields/wavefield_db_fwd.nc4";

	mNetCDF_r->openParallel(fname_wvf);

	//read nus and nrs 
	std::vector<size_t> startElem, countElem;
	startElem.push_back(mStartElem);
	countElem.push_back(mCountElem);
	Nus.assign(mCountElem, 0);
	Nrs.assign(mCountElem, 0);
	
	mNetCDF_r->readVariableChunk("Nus", Nus, startElem, countElem);
	mNetCDF_r->readVariableChunk("Nrs", Nrs, startElem, countElem);

	// create start and count for elemNu
	int totNuProc = 0;
	for (int i = 0; i<Nus.size(); i++) totNuProc+=Nus[i];
	int temp_startElemNu;
	std::vector<int> temp_countElemNu(XMPI::nproc(),0);
	XMPI::gather(totNuProc, temp_countElemNu, true);
	temp_startElemNu = std::accumulate(temp_countElemNu.begin(), temp_countElemNu.begin() + XMPI::rank(), 0);
	
	mCountElemNuFwd = temp_countElemNu[XMPI::rank()];
	mStartElemNuFwd = temp_startElemNu;
	
	std::vector<size_t> startElemNu, countElemNu;
	
	startElemNu.push_back(0);
	startElemNu.push_back(temp_startElemNu);
	startElemNu.push_back(0);
	startElemNu.push_back(0);
	startElemNu.push_back(0);
	
	countElemNu.push_back(1); // have to load one by one because of netcdf unresolved issue 
	countElemNu.push_back(temp_countElemNu[XMPI::rank()]);
	countElemNu.push_back(6);
	countElemNu.push_back(nPntEdge);
	countElemNu.push_back(nPntEdge);
	// fill disp with 0 
	vec_ar6_RMatPP initBuf(mCountElemNuFwd, zero_ar6_RMatPP);
	for (int it = 0; it < mTotSteps; it++) {//have to read one by one 
		mNetCDF_r->readVariableChunk("displacement_wavefield", initBuf, startElemNu, countElemNu);
		disp.push_back(initBuf);
		startElemNu[0]++;
	}
	
	mNetCDF_r->close();
	
	
	
	
			
}

void KernerIO::loadMaterial(vec_ar12_RMatPP &materials, std::vector<int> &Nus) {
	
	std::string fname_wvf = Parameters::sOutputDirectory + "/wavefields/wavefield_db_fwd.nc4";
	
	mNetCDF_r->openParallel(fname_wvf);

	//read nus 
	std::vector<size_t> startElem, countElem;
	startElem.push_back(mStartElem);
	countElem.push_back(mCountElem);
	Nus.assign(mCountElem, 0);
	
	mNetCDF_r->readVariableChunk("Nus", Nus, startElem, countElem);
	
	int totNuProc = 0;
	for (int i = 0; i<Nus.size(); i++) totNuProc+=Nus[i];

	
	std::vector<size_t> startElemNu, countElemNu;
	
	startElemNu.push_back(mStartElemNuFwd);
	startElemNu.push_back(0);
	startElemNu.push_back(0);
	startElemNu.push_back(0);
	
	countElemNu.push_back(mCountElemNuFwd);
	countElemNu.push_back(12);
	countElemNu.push_back(nPntEdge);
	countElemNu.push_back(nPntEdge);

	// fill with 0
	materials.assign(totNuProc, zero_ar12_RMatPP);
	vec_ar2_RMatPP vp(totNuProc,zero_ar2_RMatPP);
	
	mNetCDF_r->readVariableChunk("material_fields", materials, startElemNu, countElemNu);

	
	mNetCDF_r->close();
	
}