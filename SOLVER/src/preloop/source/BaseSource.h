// Source.h
//  base class of Source
// we only consider point Sources located on the axis

#pragma once

#include "eigenp.h"
#include "eigenc.h"

class Quad;
class Mesh;
class Domain;
class Parameters;

class Source {
public:
    Source(double depth = 0., double lat = 90., double lon = 0.);
    
    virtual ~Source() {};
    
    void release(Domain &domain, const Mesh &mesh) const;
    
    virtual std::string verbose() const = 0;
        
    double getLatitude() const {return mLatitude;};
    double getLongitude() const {return mLongitude;};
    double getDepth() const {return mDepth;};
	double getThetaSrc() const {return mThetaSrc;};
	double getPhiSrc() const {return mPhiSrc;};
    
protected:
	// on axis
	virtual void computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
        arPP_CMatX3 &fouriers) const = 0;
	
	// off axis	
	virtual void computeSourceFourier(const Quad &myQuad, 
		const RDColP &interpFactXii,
		const RDColP &interpFactEta,
		double phi,
		vec_arPP_CMatX3 &fouriers) const = 0;
        
    double mDepth;
    double mLatitude;
    double mLongitude;
	
	// theta and phi in source-centered coordinate system, for off axis source
	double mThetaSrc;
	double mPhiSrc;
	
	// on/off axis 
	bool mAxial;
        
private:
	// on axis
	bool locate(const Mesh &mesh, int &locTag, RDColP &interpFactZ) const;
	bool locate(const Mesh &mesh, int &locTag, 
		RDColP &interpFactXii, RDColP &interpFactEta) const;
	
};

