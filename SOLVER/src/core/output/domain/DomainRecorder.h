// Alex 
// Recorder for domain wide wavefields

#pragma once 

#include 

class DomainRecorder{
	
public:
	
	DomainRecorder();
	~DomainRecorder();
	
	// before time loop
    void initialize();

    // after time loop
    void finalize();

    // record at a time step
    void record(int tstep, Real t);

    // dump to netcdf
    void dumpToFile();

private:
	
	DomainIO* mIOinv; //dumps wavefields for inversion
	DomainIO* mIOani; //dumps wavefields for animation
	
	
	
	
	
	
	
};