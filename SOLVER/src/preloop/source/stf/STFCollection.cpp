// STFCollection.cpp
// source time function collection 

#include "STFCollection.h"
#include "SourceTimeFunction.h"
#include "Domain.h"
#include "XMPI.h"
#include "Parameters.h"
#include <fstream>
#include "ErfSTF.h"
#include "GaussSTF.h"
#include "RickerSTF.h"
#include "NetCDF_Reader.h"

STFCollection::STFCollection(double hdur, double duration, std::string mstf, double dt, int enforceMaxSteps, std::string offaxis_file, bool kernels) {
	
	if (hdur < 5. * dt) {
		hdur = 5. * dt;
	}
	double decay = 1.628;
	
	if (!kernels) { //if no kernels, only axis source
		
		STF *stf;
		if (boost::iequals(mstf,"erf")) {
			// Heaviside
			stf = new ErfSTF(dt, duration, hdur, decay);
		} else if (boost::iequals(mstf,"gauss")) {
			// Gaussian
			stf = new GaussSTF(dt, duration, hdur, decay);
		} else if (boost::iequals(mstf,"ricker")){
			// Ricker
			stf = new RickerSTF(dt, duration, hdur, decay);
		} else {
			throw std::runtime_error("STF::buildInparam || Unknown stf type: " + mstf);
		}
		
		// max total steps
		int maxTotalSteps = INT_MAX;
		if (enforceMaxSteps > 0) {
			maxTotalSteps = enforceMaxSteps;
		}
		if (stf->mSTFs.size() > maxTotalSteps) {
			stf->mSTFs.resize(maxTotalSteps);
			stf->mSTFp.resize(maxTotalSteps);
			stf->mSTFz.resize(maxTotalSteps);

		}
		
		mSTFs.push_back(stf);
		
	} else { //if kernels, create STFs for off axis sources, and no axis source. 
		
		int num_sources = 0;
		
		// get number of receivers, i.e. off axis sources 
		// get names of receivers to read the seismogram 
		std::vector<std::string> name, network;

		if (!boost::iequals(offaxis_file, "none")) {
			std::string input_file = Parameters::sInputDirectory + "/" + offaxis_file;
			if (XMPI::root()) {
				std::fstream fs(input_file, std::fstream::in);
				if (!fs) {
					throw std::runtime_error("STFCollection::STFCollection || " 
						"Error opening off_axis sources file " + input_file + ".");
				}
				std::string line;
				while (getline(fs, line)) {
					try {
						std::vector<std::string> strs = Parameters::splitString(line, "\t ");
						if (strs.size() < 5 || strs.size() > 6) {
							continue;
						}
						name.push_back(strs[0]);
						network.push_back(strs[1]);
						num_sources++;
					} catch(std::exception) {
						// simply ignore invalid lines
						continue;
					}
				}
				fs.close();
			}
			XMPI::bcast(name);
			XMPI::bcast(network);
			XMPI::bcast(num_sources);

		}
		
		// get downsampling factor from seismograms to wavefields 
		// get both dumping intervals, ratio is downsampling.
		int ratio = 0; 
		std::string fname_seismo = Parameters::sOutputDirectory + "/stations/axisem3d_synthetics.nc";
		std::string fname_wvf = Parameters::sOutputDirectory + "/wavefields/wavefield_db_fwd.nc4";
		
		NetCDF_Reader ncr = NetCDF_Reader();

		if (XMPI::root()) {
			
			

			
			int interval_seismo = 0;
			int interval_wvf = 0;
			
			ncr.open(fname_seismo);
			ncr.getAttribute("", "record_interval", interval_seismo);
			ncr.close();
			
			ncr.open(fname_wvf);
			ncr.getAttribute("", "record_interval", interval_wvf);
			ncr.close();
			
			ratio = interval_wvf / interval_seismo;
			
			
		}
		XMPI::bcast(ratio);
		
		// for each off axis source, load time reversed seismogram to act as STF 
		
		ncr.openParallel(fname_seismo);

		for (int i = 0; i < num_sources; i++) { // TODO : STF FROM SEISMOGRAMS 
			
			std::vector<size_t> dims;
			std::string key = network[i] + "." + name[i] + ".RTZ";
			//get dims of seismogram 
			ncr.getVarDimensions(key, dims);
			//define trace
			RMatX3 trace(dims[0], dims[1]);
			//read trace 
			ncr.read2D(key, trace);

			STF *stf = new GaussSTF(dt, duration, hdur, decay);
			
			// max total steps
			int maxTotalSteps = INT_MAX;
			if (enforceMaxSteps > 0) {
				maxTotalSteps = enforceMaxSteps;
			}
			if (stf->mSTFs.size() > maxTotalSteps) {
				stf->mSTFs.resize(maxTotalSteps);
				stf->mSTFp.resize(maxTotalSteps);
				stf->mSTFz.resize(maxTotalSteps);			}
			
			mSTFs.push_back(stf);
		}
	}
	

	
}

STFCollection::~STFCollection() {
	
	for (const auto &stf : mSTFs) {
		delete stf;
	}
}

void STFCollection::release(Domain &domain) const {
	
	for (int i = 0; i < mSTFs.size(); i++) {
    	mSTFs[i]->release(domain);
	}
}

void STFCollection::buildInparam(STFCollection *&stf, const Parameters &par, double dt, int verbose) {
    if (stf) {
        delete stf;
    }
    std::string cmtfile = Parameters::sInputDirectory + "/" + par.getValue<std::string>("SOURCE_FILE");
	
    double hdur = par.getValue<double>("SOURCE_STF_HALF_DURATION");
	double duration = par.getValue<double>("TIME_RECORD_LENGTH");
	std::string mstf = par.getValue<std::string>("SOURCE_TIME_FUNCTION");
	int enforceMaxSteps = par.getValue<int>("DEVELOP_MAX_TIME_STEPS");
	bool kernels = par.getValue<bool>("COMPUTE_KERNELS");
	std::string offaxis_file = par.getValue<std::string>("OUT_STATIONS_FILE");

	stf = new STFCollection(hdur, duration, mstf, dt, enforceMaxSteps, offaxis_file, kernels);
	
	if (verbose) {
		
		XMPI::cout<<stf->verbose();
		
	}
}

std::string STFCollection::verbose() {
	
	std::stringstream ss;
    ss << "\n========================= STFs ========================" << std::endl;
    ss << "  Number of STFs   =   " << mSTFs.size() << std::endl;
    if (mSTFs.size() > 0) {
        ss << "  STFs List: " << std::endl;
        ss << "    " << mSTFs[0]->verbose() << std::endl;
    }
    if (mSTFs.size() > 1) {
        ss << "    " << mSTFs[mSTFs.size() - 1]->verbose() << std::endl;
    }
    ss << "========================= STFs ========================\n" << std::endl;
    return ss.str();
}
