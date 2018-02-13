//KernerIO.cpp

#include "KernerIO.h"
#include "NetCDF_Writer.h"
#include "NetCDF_Reader.h"
#include "Parameters.h"
#include "XMPI.h"

void KernerIO::initialize(int dumpTimeKernels, int numFilters, int temp_startElem, int temp_countElem, int totSteps) {
	
	mTotSteps = totSteps;
	mStartElem = temp_startElem;
	mCountElem = temp_countElem;
	
	mNetCDF_w = new NetCDF_Writer();
	mNetCDF_r = new NetCDF_Reader();
	
	std::string fname_ker = Parameters::sOutputDirectory + "/kernels/kernels_db.nc4";
	std::string fname_wvf = Parameters::sOutputDirectory + "/wavefields/wavefield_db_fwd.nc4";

	if (XMPI::root()) {
		
		// get length of dimensions of fwd field : kernel will have the same for nus  
		mNetCDF_r->open(fname_wvf);
		
		std::vector<size_t> dims;
		mNetCDF_r->getVarDimensions("displacement_wavefield", dims);
		
		mNetCDF_r->close();
		
		//create dims and vars in kernel file 
		mNetCDF_w->open(fname_ker, true);
		
		std::vector<size_t> dims_notime; // for time integrated kernels 
		std::vector<size_t> dims_time; // for time kernels, optional. only allow to store one kernel in time

		dims_notime.push_back( dims[0] );
		//dims_notime.push_back( numFilters );
		dims_notime.push_back( dims[1] );
		dims_notime.push_back( 12 ); //real and imag parts of elastic radial ani kernels. 
		dims_notime.push_back( nPntEdge );
		dims_notime.push_back( nPntEdge );

		dims_time.push_back( dims[0] );
		dims_time.push_back( dims[1] );
		dims_time.push_back( 2 ); //real and imag 
		dims_time.push_back( nPntEdge );
		dims_time.push_back( nPntEdge );
		
		mNetCDF_w->defModeOn();
		mNetCDF_w->defineVariable<Real>("Kernels", dims_notime);
		if (dumpTimeKernels) mNetCDF_w->defineVariable<Real>("Kernels_time", dims_time);
		
		mNetCDF_w->defModeOff();
		mNetCDF_w->close();
		
	}

	
	
}

void KernerIO::finalize() {
	

	delete mNetCDF_w;
	delete mNetCDF_r;
	
}

void KernerIO::dumpToFile(vec_vec_ar12_RMatPP &kernels, int numFilters) {
	
	std::vector<size_t> startKernels, countKernels;
	startKernels.push_back(0);
	startKernels.push_back(mStartElemNu);
	startKernels.push_back(0);
	startKernels.push_back(0);
	startKernels.push_back(0);
	
	countKernels.push_back(1); //have to write one by one because of netcdf issue
	countKernels.push_back(mCountElemNu);
	countKernels.push_back(12);
	countKernels.push_back(nPntEdge);
	countKernels.push_back(nPntEdge); 
	
	std::string fname = Parameters::sOutputDirectory + "/kernels/kernels_db.nc4";
	mNetCDF_w->openParallel(fname);
	
	for (int ifilt = 0; ifilt < mTotSteps ; ifilt ++) {
		
		mNetCDF_w->writeVariableChunk("Kernels", kernels[ifilt], startKernels, countKernels);
		startKernels[0]++;
		
	}
	
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
	
	mCountElemNu = temp_countElemNu[XMPI::rank()];
	mStartElemNu = temp_startElemNu;
	
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
	vec_ar6_RMatPP initBuf(mCountElemNu, zero_ar6_RMatPP);
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
	
	startElemNu.push_back(mStartElemNu);
	startElemNu.push_back(0);
	startElemNu.push_back(0);
	startElemNu.push_back(0);
	
	countElemNu.push_back(mCountElemNu);
	countElemNu.push_back(12);
	countElemNu.push_back(nPntEdge);
	countElemNu.push_back(nPntEdge);

	// fill with 0
	materials.assign(totNuProc, zero_ar12_RMatPP);
	vec_ar2_RMatPP vp(totNuProc,zero_ar2_RMatPP);
	
	mNetCDF_r->readVariableChunk("material_fields", materials, startElemNu, countElemNu);

	
	mNetCDF_r->close();
	
}