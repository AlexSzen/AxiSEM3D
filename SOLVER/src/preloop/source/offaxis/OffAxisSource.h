// OffAxisSource.h
// created by Kuangdai on 11-Nov-2017 
// base class of off-axis source

#pragma once

#include "eigenp.h"
#include "eigenc.h"

class Quad;
class Mesh;
class Domain;
class Parameters;

class OffAxisSource {
public:
    OffAxisSource(double depth, double lat, double lon,
		double srcLat, double srcLon, double srcDep);
    
    virtual ~OffAxisSource() {};
    
    void release(Domain &domain, const Mesh &mesh) const;
    
    virtual std::string verbose() const = 0;
        
    double getLatitude() const {return mLatitude;};
    double getLongitude() const {return mLongitude;};
    double getDepth() const {return mDepth;};
	double getThetaSrc() const {return mThetaSrc;};
	double getPhiSrc() const {return mPhiSrc;};
    
    static void buildInparam(std::vector<OffAxisSource> *&offsrc, 
		const Parameters &par, int verbose);
    
protected:
        
    double mDepth;
    double mLatitude;
    double mLongitude;
	// theta and phi in source-centered coordinate system
	double mThetaSrc;
	double mPhiSrc;
        
private:
    bool locate(const Mesh &mesh, int &locTag, 
		RDColP &interpFactXii, RDColP &interpFactEta) const;
};

