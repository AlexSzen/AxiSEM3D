// kernel computer (kerner)

#pragma once 

#include <vector>
#include <string>
#include "eigenc.h"
#include "MultilevelTimer.h"

class DomainRecorder;
class KernerElement;
class KernerIO;
//class Processor;


class Kerner {
	
public:
	Kerner(bool dumpTimeKernels, int totSteps, int bufferSize, int recInterval, int maxStep);
	~Kerner();
	
	void initialize();
	void finalize();
	void dumpToFile();
	
	void setDomainRecorder(DomainRecorder *recorderDM) {mDomainRecorder = recorderDM;};
	void addKernerElement(KernerElement *kerElem) {mKerElements.push_back(kerElem);};
	//int getSize() {return mKerElements.size();};
	void computeKernels( int tstep );
	const int getNumElements() const {return mKerElements.size();};
	
private:
	
	void distributeFwdWvfToElements();
	void distributeBwdWvfToElements();
	void distributeMaterialToElements();
	
	
	DomainRecorder *mDomainRecorder;
	KernerIO *mIO;
	std::vector<KernerElement*> mKerElements;
	vec_vec_ar6_RMatPP mPhysicalKernels; // real and imag part of kernels for all filters and whole domain.
	std::vector<int> mNusKernel, mNrsKernel;
	// time kernels options 
	bool mDumpTimeKernels;
	
	int mMaxNr; //used for fftw init 
	int mTotSteps; // also for fftw (and others)
	int mBufferSize;
	int mRecordInterval;
	int mMaxStep; //max step of time loop
};