// Alex 
// Recorder for domain wide wavefields
// rewritten for kuangdai's implementation

#pragma once 

#include "eigenc.h"

class DomainIO;
class DomainInfo;
class Element;

class DomainRecorder{

friend class Kerner; //to allow usage of mBufferDisp. 

public:
	
	DomainRecorder(int totalRecordSteps, int recordInterval, int bufferSize, bool write);
	~DomainRecorder();
	
	void addElement(const Element* elem);
	// before time loop
    void initialize();

    // after time loop
    void finalize();

    // record at a time step
    void record(int tstep, Real t);

    // dump to netcdf
    void dumpToFile();
	
private:
	
	DomainIO* mIO; //dumps wavefields for inversion
	
	std::vector<DomainInfo> mDomainInfo;
	
	// interval
	int mTotalRecordSteps;
	int mRecordInterval;

	// buffer
	int mBufferSize;
	int mBufferLineNu;
	int mBufferLineTime;
	vec_vec_ar6_RMatPP mBufferDisp;
	RColX mBufferTime;
	
	// write file or just keep in ram 
	bool mWrite;
	
	
	
	
	
	
	
	
};