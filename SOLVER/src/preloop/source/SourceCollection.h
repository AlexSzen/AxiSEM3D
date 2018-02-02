// SourceCollection.h


#pragma once

#include "eigenp.h"
#include "eigenc.h"

class Quad;
class Mesh;
class Domain;
class Parameters;
class Source; 

class SourceCollection {
public:
    SourceCollection(std::string axis_file, std::string axis_type, std::string offaxis_file, bool kernels);
    
    ~SourceCollection();
    
    void release(Domain &domain, const Mesh &mesh) const;
    
    std::string verbose();
  
    static void buildInparam(SourceCollection *&src, const Parameters &par, int verbose);
	
	double getLatitude() const {return mLatitude;};
	double getLongitude() const {return mLongitude;};
	double getDepth() const {return mDepth;};
    
private:

	std::vector<Source *> mSources;
	
	// In forward simulation this is actual coords of source.
	// In adjoint, this is coordinate of forward source, so that domain can be rotated correctly.
	double mDepth;
	double mLatitude;
	double mLongitude;
	
};

