// create various tapers 

#pragma once 

#include "eigenc.h"
#include "eigenp.h"
#include <iostream>
#include "XMPI.h"
class Tapers {
	
public:
	
	static void cosineTaper(RColX &taper, int len) { // taper all trace to avoid fft problems 
		
		
		
		taper.resize(len);
		int cut = (5 * len)/100; //number of timesteps on which to taper : 5% of tot length
		if (cut < 2) cut = 2; 
		for (int i = 0; i < len; i++) {
			taper(i) = one;
		}
		
		for (int i=0;i<cut;i++) {
			taper(i) = cos( (3*pi/2)  + (pi/2)*( (Real) i/(cut-1) ) );
			taper(len-1-i) = taper(i);
		}
		
	};
	
	static void cosineTaper(RColX &taper, const RColX &t, Real begWin, Real endWin) { // taper in a given time window 

		taper.resize(t.size());
		taper.setZero();
		
		//find indices of beg and start of time window and size of window
		int ind_start = 0, ind_end = 0;
		
		for (int it = 0; it < t.size(); it++) {
			
			if (it != 0) { // start right before and end right after time window 
				if (t(it - 1) <= begWin && t(it) >= begWin)
					ind_start = it - 1;
				if (t(it - 1) <= endWin && t(it) >= endWin)
					ind_end = it;	
			}				
								
		}
		
		int winlen = ind_end - ind_start;
		
		for (int i = 0; i < winlen; i++) {
			taper(ind_start + i) = one;
		}
		
		// cosine taper 5% of window length before and after. make sure we dont cause undefined behavior 
		int cut = (5 * winlen) / 100;
		if (cut < 2) cut = 2;
							
		for (int i = 0; i < cut; i++) {
			taper(ind_start + i) = cos( (3*pi/2)  + (pi/2)*( (Real) i/(cut-1) ) );
			taper(ind_end - i) = taper(ind_start + i);
		}
		

	};

	
};