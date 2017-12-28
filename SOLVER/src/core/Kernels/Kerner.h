// kernel computer (kerner)

#pragma once 

#include <vector>
#include <string>
class DomainRecorder;
class KernerElement;
class KernerIO;
//class Processor;


class Kerner {
	
public:
	Kerner(bool dumpTimeKernels, int numKernels, std::vector<std::string> kerTypes, int totSteps, int maxNr);
	~Kerner();
	
	void initialize();
	void finalize();
	
	void setDomainRecorder(DomainRecorder *recorderDM) {mDomainRecorder = recorderDM;};
	void addKernerElement(KernerElement *kerElem) {mKerElements.push_back(kerElem);};
	//int getSize() {return mKerElements.size();};
	void computeKernels();
	const int getNumElements() const {return mKerElements.size();};
	
private:
	
	void distributeFwdWvfToElements();
	void distributeBwdWvfToElements();
	void distributeMaterialToElements();
	
	DomainRecorder *mDomainRecorder;
	KernerIO *mIO;
//	Processor *mPrc;
	std::vector<KernerElement*> mKerElements;
	
	// time kernels options 
	bool mDumpTimeKernels;

	// integrated kernels 
	int mNumKernels;
	std::vector<std::string> mKerTypes;
	
	int mMaxNr; //used for fftw init 
	int mTotSteps; // also for fftw (and others)
};