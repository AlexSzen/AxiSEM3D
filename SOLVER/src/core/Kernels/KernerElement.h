// elements of kerner that perform all element wise operations
// to obtain kernels

#pragma once
 
#include "eigenc.h"
#include "Element.h"

class KernerElement {
	
public:
	KernerElement(const Element *elem);
	void computeKernels(bool dumpTimeKernels);
	void feedTimeKernelsBuffer();
	
	void setForwardDisp(const vec_vec_ar3_CMatPP disp) {mForwardDisp = disp;}; //not passed by ref. we actually make a copy because we then clear the global field.
	void setBackwardDisp(const vec_vec_ar3_CMatPP disp) {mBackwardDisp = disp;}; 
	void setNuForward(const int nu) {mNuForward = nu;};
	void setMaterials(const vec_ar6_CMatPP mat) {mMaterials = mat;};
	
	const int getNuForward() const {return mNuForward;};
	const int getNuBackward() const {return mElement->getMaxNu();} ;
	
private:
	const Element *mElement;
	int mNuForward;
	
	// fields 
	vec_vec_ar3_CMatPP mForwardDisp;
	vec_vec_ar3_CMatPP mBackwardDisp;
	
	// material fields. numbering is rho, vph, vpv, vsh, vsv, eta. 
	vec_ar6_CMatPP mMaterials;
	
	// base kernels (time integrated). 
	// numbering is rho, lambda, mu, a, b, c. 
	vec_ar6_CMatPP mBaseKernels;
	
	// physical kernels (time integrated)
	// numbering is rho, vsh, vsv, vph, vpv, eta.
	vec_ar6_CMatPP mPhysicalKernels;
};