// AxisSource.h
// created by Kuangdai on 8-May-2016 
// base class of AxisSource
// we only consider point AxisSources located on the axis

#pragma once

#include "eigenp.h"
#include "eigenc.h"

class Quad;
class Mesh;
class Domain;
class Parameters;

class AxisSource : public Source {
public:
    AxisSource(double depth = 0., double lat = 90., double lon = 0.);
    
    virtual ~AxisSource() {};
    
    void release(Domain &domain, const Mesh &mesh) const;
    
    virtual std::string verbose() const = 0;
        
    double getLatitude() const {return mLatitude;};
    double getLongitude() const {return mLongitude;};
    double getDepth() const {return mDepth;};
    
    
protected:

    double mDepth;
    double mLatitude;
    double mLongitude;
        
private:
    bool locate(const Mesh &mesh, int &locTag, RDColP &interpFactZ) const;
};

