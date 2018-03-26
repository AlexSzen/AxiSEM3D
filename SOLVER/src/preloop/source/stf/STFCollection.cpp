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
#include "SeismogramSTF.h"
#include "NetCDF_Reader.h"
#include "eigenc.h"
#include "PreloopFFTW_time.h"
#include "Processor.h"

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
					} catch(std::exception) {
						// simply ignore invalid lines
						continue;
					}
				}
				fs.close();
			}
			XMPI::bcast(name);
			XMPI::bcast(network);

		}
		
		

		std::string fname_seismo = Parameters::sOutputDirectory + "/stations/axisem3d_synthetics_fwd.nc";
		std::string fname_adjoint_inputs = Parameters::sInputDirectory + "/adjoint_stations.nc4";
		
		
		
		NetCDF_Reader ncr = NetCDF_Reader(); // loads forward seismograms 
		NetCDF_Reader ncr_params = NetCDF_Reader(); // loads params for adjoint sources 
		
		ncr.openParallel(fname_seismo);
		ncr_params.openParallel(fname_adjoint_inputs);
		
		// number of adjoint sources 
		int num_sources = 0;
		ncr_params.getAttribute("", "num_sources", num_sources);						
		if (num_sources > network.size()) throw std::runtime_error("STFCollection : Number of adjoint sources bigger than number of stations in STATIONS file.");
		
		// type of tomography 
		std::string tomo;
		ncr_params.getAttributeString("", "measurement_type", tomo);

		// read time points
		std::vector<size_t> dim_time;
		ncr.getVarDimensions("time_points", dim_time);
		RColX time_points(dim_time[0]);
		ncr.read1D("time_points", time_points);
		// init fftw.  
		PreloopFFTW_time::initialize((int) dim_time[0]);
		
		// init processor 
		Processor::initialize((int) dim_time[0], time_points);
		
		// create filters
		int num_filters = 0;
		ncr_params.getAttribute("", "num_filters", num_filters);
		RMatXX filt_params = RMatXX(num_filters, 2); // temporary implementation, only with log gabor 
		std::vector<std::string> filter_names;

		for (int ifilt = 0; ifilt < num_filters; ifilt++) {
			
			std::string ifilt_name = "band_" + std::to_string(ifilt); 
			RDColX ifilt_params(2); // generated in python so need to read with double.
			ncr_params.read1D(ifilt_name,ifilt_params);
			filt_params.row(ifilt) = ifilt_params.cast<Real>(); 


		}
		
		Processor::createFilters(filt_params);
		
		// for each off axis source, load seismogram to act as STF 
		for (int i = 0; i < num_sources; i++) { // TODO : STF FROM SEISMOGRAMS 
			
			std::vector<size_t> dims, dims_params;
			std::string key = network[i] + "." + name[i] + ".SPZ";
			std::string key_params = network[i] + "." + name[i]; 
			
			//get dims of seismogram 
			ncr.getVarDimensions(key, dims);
			//read trace 
			RMatX3 trace(dims[0], dims[1]);
			ncr.read2D(key, trace);
			

			
			//get dims of params 
			ncr_params.getVarDimensions(key_params, dims_params);
			//read params
			RDMatXX adjoint_params(dims_params[0],dims_params[1]);
			ncr_params.read2D(key_params, adjoint_params);

			IColX filter_types = IColX(adjoint_params.rows());
			
			for (int i_measurement = 0; i_measurement < adjoint_params.rows(); i_measurement++) {
				int filter_type;
				ncr_params.getAttribute(key_params, "filter_" + std::to_string(i_measurement), filter_type);
				filter_types(i_measurement) = filter_type; // gives index of filter to use.
			} 
						
			STF *stf = new SeismogramSTF(trace, dt, duration, hdur, decay, adjoint_params, filter_types);
		//	STF *stf = new GaussSTF(dt, duration, hdur, decay);
			
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
		
		// finalize fftw 
		PreloopFFTW_time::finalize();
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
