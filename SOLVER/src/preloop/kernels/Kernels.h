// builds params for kerner and signal proc 

#pragma once 

#include <vector>
#include <string>
#include "global.h"

class Parameters;
class Domain;
class Mesh;

class Kernels {

public:
	
	Kernels(bool computeKer); // if computeKer = 0 we call this constructor 
	~Kernels();
	
	std::string verbose();
	
	static void buildInparam(Kernels *&kernels, const Parameters &par, int totalStepsSTF, int verbose);
	void release(Domain &domain, const Mesh &mesh);
	
private:
	
	bool mComputeKernels;
	bool mDumpTimeKernels;
	//filtering 
	std::string mFilter;
	std::vector<Real> mFiltParams;
	
	//tapering 
	std::string mTaper;
	
	// kernels to compute 
	int mNumKernels;
	std::vector<std::string> mKerTypes;
	
	double mRmin = 0., mRmax = 7.e6;
	double mThetaMin = 0., mThetaMax = 180.;
	
	int mTotalStepsKernels;
	
	
};