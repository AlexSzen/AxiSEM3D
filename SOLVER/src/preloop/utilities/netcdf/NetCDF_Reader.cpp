// NetCDF_Reader.h
// created by Kuangdai on 17-May-2017 
// NetCDF Reader

#include "NetCDF_Reader.h"
#include <netcdf.h>

#ifdef _USE_PARALLEL_NETCDF
    #include <netcdf_par.h>
#endif

void NetCDF_Reader::open(const std::string &fname) {
    close();
    mFileName = fname;
    if (nc_open(fname.c_str(), NC_NETCDF4, &mFileID) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::open || "
            "Error opening NetCDF file: || " + fname);
    }
}

void NetCDF_Reader::openParallel(const std::string &fname) {
    #ifdef _USE_PARALLEL_NETCDF
        close();
        mFileName = fname;
        if (nc_open_par(fname.c_str(), NC_MPIIO | NC_NETCDF4, 
            MPI_COMM_WORLD, MPI_INFO_NULL, &mFileID) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Reader::openParallel || "
                "Error opening NetCDF file: || " + fname);
        }
    #else
        throw std::runtime_error("NetCDF_Reader::openParallel || "
            "Parallel NetCDF is disabled in CMakeLists.txt");
    #endif
}

void NetCDF_Reader::close() {
    if (isOpen()) {
        netcdfError(nc_close(mFileID), "nc_close");
        mFileID = -1;
        mFileName = "";
    }
}

void NetCDF_Reader::getAttributeString(const std::string &vname, const std::string &attname, std::string &attvalue) const {
	//// Problem with dynamic allocation of chars. 
	/// If dynamically allocated char array, some lengths of attributes cause problems for some reason (e.g. if shorter than the name of the attribute)
	/// Workaround by declaring compile time array big enough to hold any string (e.g. 1000). 
	
	int varid = -1;
	int varloc = -1;
	if (vname == "") {
		varid = NC_GLOBAL;
		varloc = mFileID;
	} else {
		varid = inquireVariable(vname);
		varloc = mFileID;
	}
	
	size_t attlen = 0;
	netcdfError(nc_inq_attlen(varloc, varid, attname.c_str(), &attlen), "nc_inq_attlen");
	char cstr[1000]; // assume will never go above 1000 chars string... 
	if (nc_get_att_text(varloc, varid, attname.c_str(), cstr) != NC_NOERR) {
		throw std::runtime_error("NetCDF_Reader::getAttributeString || "
			"Error getting attribute from variable, variable: " + vname + ", attribute: " + attname  
			+ " || NetCDF file: " + mFileName);
	}
	attvalue = std::string(cstr);
		
}

void NetCDF_Reader::readString(const std::string &vname, std::vector<std::string> &data) const {
    // access variable
    int var_id = -1;
    if (nc_inq_varid(mFileID, vname.c_str(), &var_id) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::readString || "
            "Error finding variable: " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get ndims
    int var_ndims = -1;
    netcdfError(nc_inq_varndims(mFileID, var_id, &var_ndims), "nc_inq_varndims");
    if (var_ndims != 2) {
        throw std::runtime_error("NetCDF_Reader::readString || "
            "Number of dimensions is not 2, Variable = " + vname + " || NetCDF file: " + mFileName);
    }
    
    // get dim length
    std::vector<int> var_dimids(var_ndims, -1);
    std::vector<size_t> dims = std::vector<size_t>(var_ndims, 0);
    netcdfError(nc_inq_vardimid(mFileID, var_id, var_dimids.data()), "nc_inq_vardimid");
    for (int i = 0; i < var_ndims; i++) {
        netcdfError(nc_inq_dimlen(mFileID, var_dimids[i], &dims[i]), "nc_inq_dimlen");
    }
    
    // size
    int numString = dims[0];
    int lenString = dims[1];
    char *cstr = new char[numString * lenString];
    netcdfError(nc_get_var_text(mFileID, var_id, cstr), "nc_get_var_text");
    for (int i = 0; i < numString; i++) {
        data.push_back(std::string(&cstr[i * lenString]));
    }
    
    delete [] cstr;
}

void NetCDF_Reader::getVarDimensions(const std::string &vname, std::vector<size_t> &dims) const {
	
	// access variable
	int var_id = -1;
	if (nc_inq_varid(mFileID, vname.c_str(), &var_id) != NC_NOERR) {
		throw std::runtime_error("NetCDF_Reader::getVarDimensions || "
			"Error finding variable: " + vname + " || NetCDF file: " + mFileName);
	}
	
	// get ndims
	int var_ndims = -1;
	netcdfError(nc_inq_varndims(mFileID, var_id, &var_ndims), "nc_inq_varndims");
	
	// get dim length
	std::vector<int> var_dimids(var_ndims, -1);
	netcdfError(nc_inq_vardimid(mFileID, var_id, var_dimids.data()), "nc_inq_vardimid");
	dims.resize(var_ndims);
	size_t total_len = 1;
	for (int i = 0; i < var_ndims; i++) {
		netcdfError(nc_inq_dimlen(mFileID, var_dimids[i], &dims[i]), "nc_inq_dimlen");
	}
}

int NetCDF_Reader::inquireVariable(const std::string &vname) const {
    int var_id = -1;
    if (nc_inq_varid(mFileID, vname.c_str(), &var_id) != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::inquireVariable || "
            "Error finding variable: " + vname + " || NetCDF file: " + mFileName);
    }
    return var_id;
}

bool NetCDF_Reader::isNetCDF(const std::string &fname) {
    int fid, stat;
    stat = nc_open(fname.c_str(), 0, &fid);
    nc_close(fid);
    return stat == NC_NOERR;
}

void NetCDF_Reader::netcdfError(const int retval, const std::string &func_name) const {
    if (retval != NC_NOERR) {
        throw std::runtime_error("NetCDF_Reader::netcdfError || "
            "Error in NetCDF function: " + func_name + " || NetCDF file: " + mFileName);
    }
}

#include "NetCDF_ReaderAscii.h"
bool NetCDF_Reader::checkNetCDF_isAscii(const std::string &fname) {
    if (isNetCDF(fname)) {
        return false; 
    }
    
    if (NetCDF_ReaderAscii::isNetCDFAscii(fname)) {
        return true;
    }
    
    throw std::runtime_error("NetCDF_ReaderAscii::checkNetCDFAscii || "
        "Error identifying NetCDF or NetCDF-alternative ascii file: || " + fname);
}
