// OffAxisSourceCollection.h
// created by Kuangdai on 11-Nov-2017 
// base class of off-axis source

#pragma once

#include "eigenp.h"
#include "eigenc.h"

class Quad;
class Mesh;
class Domain;
class Parameters;

class OffAxisSourceCollection {
public:
    OffAxisSourceCollection(double depth, double lat, double lon,
		double srcLat, double srcLon, double srcDep);
    
    virtual ~OffAxisSourceCollection() {};
    
    void release(Domain &domain, const Mesh &mesh) const;
    
    virtual std::string verbose() const = 0;
        

    
    static void buildInparam(std::vector<OffAxisSource> *&offsrc, 
		const Parameters &par, int verbose);
    
private:
	
	std::vector<OffAxisSource *> mSources;

};

