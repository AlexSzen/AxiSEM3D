// kernels.cpp
//build params for kernel computation 

#include "Kernels.h"
#include "Domain.h"
#include "Mesh.h"
#include "Element.h"
#include "Parameters.h"
#include "XMPI.h"
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "Kerner.h"
#include "KernerElement.h"

Kernels::Kernels(bool computeKer): mComputeKernels(computeKer) {
	
}


Kernels::~Kernels() {
	
}

void Kernels::buildInparam(Kernels *&kernels, const Parameters &par, int totalStepsSTF,  int verbose) {
	
	bool computeKer = par.getValue<bool>("COMPUTE_KERNELS");

	kernels = new Kernels(computeKer);
	if (!computeKer) return;

	std::string filt = par.getValue<std::string>("FILTER_TYPE");
	std::vector<Real> filtParams;
	
	if (boost::iequals(filt, "log_gabor")) {
			
			Real center = par.getValue<Real>("CENTER_FREQ");
			filtParams.push_back(center);
			
			Real sigma = par.getValue<Real>("SIGMA");
			filtParams.push_back(sigma);
		
		
	} else throw std::runtime_error("Parameters :: filter " + filt + " not implemented");
	
	kernels->mFilter = filt;
	kernels->mFiltParams = filtParams;
	
	kernels->mTaper = par.getValue<std::string>("TAPER_TYPE");
	
	kernels->mNumKernels = par.getSize("KERNEL_TYPES");
	for (int i = 0; i < kernels->mNumKernels; i++) {
		kernels->mKerTypes.push_back(par.getValue<std::string>("KERNEL_TYPES", i));
	}
			
	kernels->mRmin = par.getValue<double>("RMIN")*1e3;
	kernels->mRmax = par.getValue<double>("RMAX")*1e3;
	kernels->mThetaMin = par.getValue<double>("THETA_MIN")*degree;
	kernels->mThetaMax = par.getValue<double>("THETA_MAX")*degree;
	kernels->mDumpTimeKernels = par.getValue<bool>("OUTPUT_TIME_KERNELS");

	int recInterval = par.getValue<int>("WAVEFIELD_RECORD_INTERVAL");
	if (recInterval <= 0) {
		recInterval = 1;
	}
	
	kernels->mTotalStepsKernels = totalStepsSTF / recInterval;
	if (totalStepsSTF % recInterval > 0) {
		kernels->mTotalStepsKernels += 1;
	}
	
	if (verbose) XMPI::cout << kernels->verbose();

}

void Kernels::release(Domain &domain, const Mesh &mesh) {
	
	if (!mComputeKernels) return;
	
	Kerner *kerner = new Kerner(mDumpTimeKernels, mNumKernels, mKerTypes, mTotalStepsKernels, mesh.getMaxNr());  
	
	for (int ielem = 0; ielem < domain.getNumElements(); ielem ++) {
		Element *elem = domain.getElement(ielem);
		if ( elem->needDumping(mRmin,mRmax,mThetaMin,mThetaMax) ) {
			KernerElement *kerElem = new KernerElement(elem);
			kerner->addKernerElement(kerElem);
		}
	}
	kerner->setDomainRecorder(domain.getDomainRecorder());
	domain.setKerner(kerner);
	
}

std::string Kernels::verbose() {
	
	std::stringstream ss;
	ss << "\n========================= Kernels ========================" << std::endl;
	ss << "Kernels to be computed : " ; 
	for (int i = 0; i<mNumKernels; i++) {
		ss << mKerTypes[i] << " ";
	}
	ss << std::endl;
	ss << "Filter : " << mFilter << std::endl;
	ss << "Taper : " << mTaper << std::endl;
	ss << "========================= Kernels ========================\n" << std::endl;
	return ss.str();
}


