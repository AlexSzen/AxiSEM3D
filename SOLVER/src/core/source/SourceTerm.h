// SourceTerm.h
// created by Kuangdai on 5-Apr-2016 
// base class of sources in the solver
// NOTE: there are many different types of seismic sources, such as earthquakes,
// explosions, surface forces (ocean), and point forces (adjoint). But to the 
// solver, all those sources are seen as a combination of force vectors on solid 
// points, or of scalars on fluid points if any.
// So, if we only consider sources in the solid domain, only one class is enough
// to implement all source types. 


#pragma once

#include "eigenc.h"

class Element;

class SourceTerm {
public:
    SourceTerm(Element *element, const arPP_CMatX3 &force, int isource);
    
    // apply source at a time step 
    void apply(Real stfs, Real stfp, Real stfz);
    
    void setElement(Element *element) {mElement = element;};
	
	int getNumSource() {return mISource;}; // get the source number to know which stf to apply
    
private:
    
    Element *mElement;
    arPP_CMatX3 mForce;
    arPP_CMatX3 mForceXSTF;
	
	// because we can have many sources with many stfs,
	// to which source does this term belong 
	int mISource;
};