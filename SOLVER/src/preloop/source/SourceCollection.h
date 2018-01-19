// SourceCollection.h
// created by Kuangdai on 11-Nov-2017 
// base class of off-axis source

#pragma once

#include "eigenp.h"
#include "eigenc.h"

class Quad;
class Mesh;
class Domain;
class Parameters;

class SourceCollection {
public:
    SourceCollection(double depth, double lat, double lon,
		double srcLat, double srcLon, double srcDep);
    
    virtual ~SourceCollection() {};
    
    void release(Domain &domain, const Mesh &mesh) const;
    
    virtual std::string verbose() const = 0;
        

    
    static void buildInparam(const Parameters &par, int verbose);
    
private:
	
	Source *mSource;
	std::vector<OffAxisSource *> mOffAxisSources;
	
};

