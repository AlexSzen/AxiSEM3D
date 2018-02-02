// NetCDF_Reader.h
// created by Kuangdai on 17-May-2017 
// NetCDF Reader

#pragma once

#include <netcdf.h>
#include <vector>
#include <string>
#include <stdexcept>

class NetCDF_Reader {
public:    
    // file
    void open(const std::string &fname);
    void openParallel(const std::string &fname);
    void close();
    bool isOpen() const {return mFileName != "";};
    
    // double
    template<class Container>
    void readMetaData(const std::string &vname, Container &data, std::vector<size_t> &dims) const {
        // access variable
        int var_id = -1;
        if (nc_inq_varid(mFileID, vname.c_str(), &var_id) != NC_NOERR) {
            throw std::runtime_error("NetCDF_Reader::readMetaData || "
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
            total_len *= dims[i];
        }
        
        // get data
        data.resize(total_len);
        netcdfError(nc_get_var(mFileID, var_id, data.data()), "nc_get_var");
    };
    
    template<class Container>
    void read1D(const std::string &vname, Container &data) const {
        // read meta data
        std::vector<size_t> dims;
        readMetaData(vname, data, dims);
        
        // check ndims
        int var_ndims = dims.size();
        if (!(var_ndims == 1 || (var_ndims == 2 && dims[0] == 1))) {
            throw std::runtime_error("NetCDF_Reader::read1D || "
                "Variable is not 1D, Variable = " + vname + " || NetCDF file: " + mFileName);
        }
    };
    
    template<class Container>
    void read2D(const std::string &vname, Container &data) const {
        // read meta data
        std::vector<size_t> dims;
        std::vector<typename Container::Scalar> mdata;
        readMetaData(vname, mdata, dims);
        
        // check ndims
        int var_ndims = dims.size();
        if (var_ndims != 2) {
            throw std::runtime_error("NetCDF_Reader::read2D || "
                "Variable is not 2D, Variable = " + vname + " || NetCDF file: " + mFileName);
        }
        
        // Container can be both RowMajor and ColMajor
        int pos = 0;
        data = Container::Zero(dims[0], dims[1]);
        for (int j = 0; j < dims[0]; j++) {
            for (int k = 0; k < dims[1]; k++) {
                data(j, k) = mdata[pos++];
            }
        }
    };
	
	// read chunk of a variable
	// NOTE: if the Container is an Eigen Matrix, it has to be Eigen::RowMajor
	template<class Container>   
	void readVariableChunk(const std::string &vname, Container &data, 
		std::vector<size_t> &start, std::vector<size_t> &count) const {
		int varid = inquireVariable(vname);
		if (data.size() > 0) {
			if (nc_get_vara(mFileID, varid, start.data(), count.data(), data.data()) != NC_NOERR) {
				throw std::runtime_error("NetCDF_Reader::readVariableChunk || "
					"Error reading variable, variable: " + vname + " || NetCDF file: " + mFileName);
			}
		}
	};
	
	template<class base_type>
	void getAttribute(const std::string &vname, 
		const std::string &attname, base_type &attvalue) const {
		int varid = -1;
		int varloc = -1;
		if (vname == "") {
			varid = NC_GLOBAL;
			varloc = mFileID;
		} else {
			varid = inquireVariable(vname);
			varloc = mFileID;
		}
		if (nc_get_att(varloc, varid, attname.c_str(), &attvalue) != NC_NOERR) {
			throw std::runtime_error("NetCDF_Reader::getAttribute || "
				"Error getting attribute from variable, variable: " + vname 
				+ ", attribute: " + attname  
				+ " || NetCDF file: " + mFileName);
		}
	};
    
	// get lengths of dimensions for a variable
	void getVarDimensions(const std::string &vname, std::vector<size_t> &dims) const;
	
	int inquireVariable(const std::string &vname) const;	
    // string
    void readString(const std::string &vname, std::vector<std::string> &data) const;
    
    // check
    static bool isNetCDF(const std::string &fname);
    
    // determine netcdf or ascii
    static bool checkNetCDF_isAscii(const std::string &fname);
    
private:    
    // error handler
    void netcdfError(const int retval, const std::string &func_name) const;
    
private:
    int mFileID = -1;
    std::string mFileName = "";
};

