// Source collection 

#include "SourceCollection.h"
#include "Parameters.h"
#include "Source.h"
#include "Earthquake.h"
#include "PointForce.h"
#include "OffAxisPointForce.h"
#include "XMPI.h"
#include <fstream>
#include <sstream>
#include "Source.h"
#include "NetCDF_Reader.h"



SourceCollection::SourceCollection(std::string axis_file, std::string axis_type, std::string offaxis_file, bool kernels) {
	
	
	// axis source first, if not calculating kernels 	
	if (!kernels) {
	    if (boost::iequals(axis_type, "earthquake")) {
	        std::string cmtfile = Parameters::sInputDirectory + "/" + axis_file;
	        double depth, lat, lon;
	        double Mrr, Mtt, Mpp, Mrt, Mrp, Mtp;
	        if (XMPI::root()) {
	            std::fstream fs(cmtfile, std::fstream::in);
	            if (!fs) {
	                throw std::runtime_error("Source::buildInparam || "
	                    "Error opening CMT data file: ||" + cmtfile);
	            }
	            std::string junk;
	            std::getline(fs, junk);
	            std::getline(fs, junk);
	            std::getline(fs, junk);
	            std::getline(fs, junk);
	            fs >> junk >> lat;
	            fs >> junk >> lon;
	            fs >> junk >> depth;
	            fs >> junk >> Mrr;
	            fs >> junk >> Mtt;
	            fs >> junk >> Mpp;
	            fs >> junk >> Mrt;
	            fs >> junk >> Mrp;
	            fs >> junk >> Mtp;
	            depth *= 1e3;
	            Mrr *= 1e-7;
	            Mtt *= 1e-7;
	            Mpp *= 1e-7;
	            Mrt *= 1e-7;
	            Mrp *= 1e-7;
	            Mtp *= 1e-7;
	            fs.close();
	        }
	        XMPI::bcast(depth);
	        XMPI::bcast(lat);
	        XMPI::bcast(lon);
	        XMPI::bcast(Mrr);
	        XMPI::bcast(Mtt);
	        XMPI::bcast(Mpp);
	        XMPI::bcast(Mrt);
	        XMPI::bcast(Mrp);
	        XMPI::bcast(Mtp);
	        mSources.push_back(new Earthquake(depth, lat, lon, Mrr, Mtt, Mpp, Mrt, Mrp, Mtp));
			
			mLatitude = lat;
			mLongitude = lon; 
			mDepth = depth;
			
	    } else if (boost::iequals(axis_type, "point_force")) {
	        // point force
	        std::string pointffile = Parameters::sInputDirectory + "/" + axis_file;
	        double depth, lat, lon;
	        double f1, f2, f3;
	        if (XMPI::root()) {
	            std::fstream fs(pointffile, std::fstream::in);
	            if (!fs) {
	                throw std::runtime_error("Source::buildInparam || "
	                    "Error opening point force data file: ||" + pointffile);
	            }
	            std::string junk;
	            std::getline(fs, junk);
	            std::getline(fs, junk);
	            std::getline(fs, junk);
	            std::getline(fs, junk);
	            fs >> junk >> lat;
	            fs >> junk >> lon;
	            fs >> junk >> depth;
	            fs >> junk >> f1;
	            fs >> junk >> f2;
	            fs >> junk >> f3;
	            depth *= 1e3;
	            fs.close();
	        }
	        XMPI::bcast(depth);
	        XMPI::bcast(lat);
	        XMPI::bcast(lon);
	        XMPI::bcast(f1);
	        XMPI::bcast(f2);
	        XMPI::bcast(f3);
	        mSources.push_back(new PointForce(depth, lat, lon, f1, f2, f3));
			
			mLatitude = lat;
			mLongitude = lon; 
			mDepth = depth;
			
	    } else {
	         throw std::runtime_error("Source::buildInparam || Unknown source type: " + axis_type);
	    }
	} else { // if kernels, use off axis only 
		
		std::vector<std::string> name, network;
		std::vector<double> theta, phi, depth;
		if (!boost::iequals(offaxis_file, "none")) {
			std::string input_file = Parameters::sInputDirectory + "/" + offaxis_file;
			if (XMPI::root()) {
				std::fstream fs(input_file, std::fstream::in);
				if (!fs) {
					throw std::runtime_error("SourceCollection::SourceCollection || " 
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
						theta.push_back(boost::lexical_cast<double>(strs[2]));
						phi.push_back(boost::lexical_cast<double>(strs[3]));
						depth.push_back(boost::lexical_cast<double>(strs[strs.size() - 1]));
					} catch(std::exception) {
						// simply ignore invalid lines
						continue;
					}
				}
				fs.close();
			}
			XMPI::bcast(name);
			XMPI::bcast(network);
			XMPI::bcast(theta);
			XMPI::bcast(phi);
			XMPI::bcast(depth);
		}
		
		//get coords of forward source 
		double srcDepth = 0;
		double srcLat = 0;
		double srcLon = 0;
		
		if (XMPI::root()) {
			std::string fname_fwd = Parameters::sOutputDirectory + "/wavefields/wavefield_db_fwd.nc4";
			NetCDF_Reader nc_reader = NetCDF_Reader();
			nc_reader.open(fname_fwd);
			nc_reader.getAttribute("", "source_latitude", srcLat);
			nc_reader.getAttribute("", "source_longitude", srcLon);
			nc_reader.getAttribute("", "source_depth", srcDepth);
			nc_reader.close();

		}
		XMPI::bcast(srcDepth);
		XMPI::bcast(srcLat);
		XMPI::bcast(srcLon);
		
		mLatitude = srcLat;
		mLongitude = srcLon;
		mDepth = srcDepth;

		// create sources
		std::vector<std::string> recKeys;
		for (int i = 0; i < name.size(); i++) {
			// check duplicated
			std::string key = network[i] + "." + name[i];
			if (std::find(recKeys.begin(), recKeys.end(), key) != recKeys.end()) {
	
				// error
				throw std::runtime_error("SourceCollection::SourceCollection duplicated off axis sources ");
				
			}
			recKeys.push_back(key);
			// add source
			mSources.push_back(new OffAxisPointForce(depth[i], theta[i], phi[i], srcLat, srcLon, srcDepth));

		}
		
	}
	
}

SourceCollection::~SourceCollection() {
	
	for (const auto &src : mSources) {
		delete src;
	}
	
}

void SourceCollection::release(Domain &domain, const Mesh &mesh) const {
		
	for (int isource = 0; isource < mSources.size(); isource++) {
		mSources[isource]->release(domain, mesh, isource);
	}

	
}

std::string SourceCollection::verbose() {
	
	std::stringstream ss;
    ss << "\n========================= Sources ========================" << std::endl;
    ss << "  Number of Sources   =   " << mSources.size() << std::endl;
    if (mSources.size() > 0) {
        ss << "  Sources List: " << std::endl;
        ss << "    " << mSources[0]->verbose() << std::endl;
    }
    if (mSources.size() > 1) {
        ss << "    " << mSources[mSources.size() - 1]->verbose() << std::endl;
    }
    ss << "========================= Sources ========================\n" << std::endl;
    return ss.str();
}

void SourceCollection::buildInparam(SourceCollection *&src, const Parameters &par, int verbose) {
	
	std::string axis_src_type = par.getValue<std::string>("SOURCE_TYPE");
	std::string axis_src_file = par.getValue<std::string>("SOURCE_FILE");	
	std::string offaxis_src_file = par.getValue<std::string>("OUT_STATIONS_FILE");
	bool kernels = par.getValue<bool>("COMPUTE_KERNELS");
	
	src = new SourceCollection(axis_src_file, axis_src_type, offaxis_src_file, kernels);
	
	if (verbose) {
		
		XMPI::cout<<src->verbose();
		
	}
	
}