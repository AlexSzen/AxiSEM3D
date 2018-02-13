// Source.h
// base class of Source
// modified for on and off axis sources 

#pragma once

#include "eigenp.h"
#include "eigenc.h"

class Quad;
class Mesh;
class Domain;
class Parameters;

class Source {
public:
	// on axis
    Source(double depth = 0., double lat = 90., double lon = 0.);
    // off axis 
	Source(double depth, double lat, double lon,
		double srcLat, double srcLon, double srcDep);
		
    virtual ~Source() {};
    
    void release(Domain &domain, const Mesh &mesh, int isource) const;
    
    virtual std::string verbose() const = 0;
        
    double getLatitude() const {return mLatitude;};
    double getLongitude() const {return mLongitude;};
    double getDepth() const {return mDepth;};
	
	//off axis : theta and phi of forward source
	double getThetaSrc() const {return mThetaSrc;};
	double getPhiSrc() const {return mPhiSrc;};
    
protected:
	
	// on and off axis	
	virtual void computeSourceFourier(const Quad &myQuad, 
		const RDColP &interpFactZ, 
		const RDColP &interpFactXii,
		const RDColP &interpFactEta,
		double phi,
		arPP_CMatX3 &fouriers) const = 0;
        
    double mDepth;
    double mLatitude;
    double mLongitude;
	
	// theta and phi in source-centered coordinate system, for off axis source
	double mThetaSrc = 0.;
	double mPhiSrc = 0.;
	
	// on/off axis 
	bool mAxial;

        
private:
	// on axis
	bool locate(const Mesh &mesh, int &locTag, RDColP &interpFactZ) const;
	//off axis
	bool locate(const Mesh &mesh, int &locTag,
		RDColP &interpFactXii, RDColP &interpFactEta) const;
	
};

