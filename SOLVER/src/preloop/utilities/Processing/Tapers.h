// create various tapers 

#pragma once 

#include "eigenc.h"
#include "eigenp.h"

class Tapers {
	
public:
	
	static void cosineTaper(RColX &taper, int len) {
		
		
		
		taper.resize(len);
		int cut = (5 * len)/100; //number of timesteps on which to taper : 5% of tot length
		
		for (int i = 0; i < len; i++) {
			taper(i) = one;
		}
		
		for (int i=0;i<cut;i++) {
			taper(i) = cos( (3*pi/2)  + (pi/2)*( (Real) i/(cut-1) ) );
			taper(len-1-i) = taper(i);
		}
		
	};

	
};