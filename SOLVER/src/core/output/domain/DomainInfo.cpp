// DomainInfo.cpp 

#include "DomainInfo.h"
#include "Element.h"

DomainInfo::DomainInfo(const Element* elem): mElement(elem) {
	
}

void DomainInfo::feedBuffer(int bufferLineTime, int bufferLineNu, vec_vec_ar6_RMatPP& bufferDisp) {
	
	const vec_ar3_CMatPP &disp = mElement->getDisp();
	
	for (int inu = 0; inu <= getMaxNu(); inu++) {
		
		bufferDisp[bufferLineTime][bufferLineNu + inu][0] = disp[inu][0].real();
		bufferDisp[bufferLineTime][bufferLineNu + inu][1] = disp[inu][0].imag();
		bufferDisp[bufferLineTime][bufferLineNu + inu][2] = disp[inu][1].real();
		bufferDisp[bufferLineTime][bufferLineNu + inu][3] = disp[inu][1].imag();
		bufferDisp[bufferLineTime][bufferLineNu + inu][4] = disp[inu][2].real();
		bufferDisp[bufferLineTime][bufferLineNu + inu][5] = disp[inu][2].imag();

	}
	
}

int DomainInfo::getMaxNu() const {
	return mElement->getMaxNu();
}