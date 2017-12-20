// DomainInfo.h
// created by alex.
// gets disp for each element to feed a domain wide bufferSize

#pragma once 

#include "eigenc.h"

class Element;

class DomainInfo {
	
public:
	
	DomainInfo(const Element* elem);
	~DomainInfo() {};
	
	void feedBuffer(int bufferLineTime, int bufferLineNu, vec_vec_ar6_RMatPP& bufferDisp);
	
	int getMaxNu() const;
	
	
private:
	const Element* mElement;
	
	
};