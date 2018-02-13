// create filters 

#pragma once 

#include "eigenc.h"
#include <cmath>

class Filters{
	
public:
	
	static void logGabor(RMatXX &filters, const RColX &freq, Real center, Real sigma, int indFilt) {
			

			filters(indFilt,0) = 0.;
			
			for (int i = 1; i < freq.size(); i++) {
				
				filters(indFilt,i) = exp( (-pow(log(freq(i)/center),2.) )/(2*pow(log(sigma),2.)) );
				
			}
			
		
		
	}
	
};