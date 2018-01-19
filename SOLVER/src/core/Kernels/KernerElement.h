// elements of kerner that perform all element wise operations
// to obtain kernels

#pragma once
 
#include "eigenc.h"
#include "Element.h"

class KernerElement {
	
public:
	KernerElement(const Element *elem);
	void computeKernels(bool dumpTimeKernels);
	void feedKernels(vec_vec_ar12_RMatPP &physKernels, int nuLine, int nuElem);
	void clearKernels() {mPhysicalKernels.clear(); mBaseKernels.clear();};
	void test();
	
	void setForwardDisp(const vec_vec_ar3_CMatPP disp) {mForwardDisp = disp;}; //not passed by ref. we actually make a copy because we then clear the global field.
	void setBackwardDisp(const vec_vec_ar3_CMatPP disp) {mBackwardDisp = disp;}; 
	void setNuForward(const int nu) {mNuForward = nu;};
	void setNrForward(const int nr) {mNrForward = nr; mNyquist = (int)(mNrForward % 2 == 0);};
	void setMaterials(const vec_ar6_CMatPP mat) {mMaterials = mat;};
	void setTimeAndFreqSize(int totSteps) {mTimeSize = totSteps; mFreqSize = totSteps / 2 + 1;}
	
	const int getNuForward() const {return mNuForward;};
	const int getNuBackward() const {return mElement->getMaxNu();} ;
	

	
private:
			
	const Element *mElement;
	int mNuForward;
	int mNrForward;
	int mNyquist;
	
	// time and freq len 
	int mTimeSize, mFreqSize;
	
	// fields 
	vec_vec_ar3_CMatPP mForwardDisp;
	vec_vec_ar3_CMatPP mBackwardDisp;
	
	// material fields. numbering is rho, vph, vpv, vsh, vsv, eta. 
	vec_ar6_CMatPP mMaterials;
	
	// base kernels (time integrated). 
	// numbering is rho, lambda, mu, a, b, c. 
	vec_vec_ar9_CMatPP mBaseKernels;
	
	// physical kernels (time integrated) for each filter
	// numbering is rho, vsh, vsv, vph, vpv, eta.
	vec_vec_ar6_CMatPP mPhysicalKernels;
};