

#include "Kerner.h"
#include "KernerElement.h"
#include "KernerIO.h"
#include "XMPI.h"
#include "DomainRecorder.h"
#include "eigenc.h"
#include "PreloopFFTW.h" //just for lucky number 
#include "KernerFFTW_N3.h"
#include "KernerFFTW_N6.h"
#include "KernerFFTW_N9.h"
#include "Processor.h"

Kerner::Kerner(bool dumpTimeKernels, int totSteps, int bufferSize, int recInterval, int maxStep): 
mDumpTimeKernels(dumpTimeKernels), mBufferSize(bufferSize), mRecordInterval(recInterval), mMaxStep (maxStep),
mTotSteps(totSteps)
 {
	


}

Kerner::~Kerner() {
	for (const auto &e: mKerElements) {delete e;}
	delete mIO;

}

void Kerner::initialize() {
	
	int startElem;
	std::vector<int> countElem(XMPI::nproc(), 0);
	XMPI::gather(mKerElements.size(), countElem, true);
	startElem = std::accumulate(countElem.begin(), countElem.begin()+XMPI::rank(),0);
	int totElems = XMPI::sum((int) mKerElements.size());
	
	mIO = new KernerIO(mDumpTimeKernels, startElem, countElem[XMPI::rank()], mTotSteps, mBufferSize);
	
	// distribute Nus  
	distributeNus();
	//distribute materials
	distributeMaterialToElements();
	
	
	// gather nus and initialize workspace of elements
	int totNu = 0; 
	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {

		int nuMax = mKerElements[ielem]->getNuMax();
		totNu += nuMax + 1;
		mNusKernel.push_back(nuMax + 1);
		mNrsKernel.push_back(mKerElements[ielem]->getNrMax());
	}
			
	//init kernel 
	vec_ar6_RMatPP initKernels(totNu, zero_ar6_RMatPP);
	if (mDumpTimeKernels) {
		mPhysicalKernels.assign(mBufferSize, initKernels);
	} else {
		mPhysicalKernels.assign(1, initKernels);
	}
	
	int startElemNu;
	std::vector<int> countElemNu(XMPI::nproc(),0);
	XMPI::gather(totNu, countElemNu, true);
	startElemNu = std::accumulate(countElemNu.begin(), countElemNu.begin() + XMPI::rank(), 0);
	int totTotNu = XMPI::sum(totNu);
	
	mIO->initialize(totTotNu, startElemNu, countElemNu[XMPI::rank()], totElems, mNusKernel, mNrsKernel);

	
}

void Kerner::finalize() {
	
	mIO->finalize();

}

void Kerner::computeKernels( int tstep ) {

	
	if ((tstep == mMaxStep) && (mTotSteps % mBufferSize != 0) ) {
		mBufferSize = mTotSteps % mBufferSize;
	}
	
	if ((tstep % (mBufferSize * mRecordInterval) == 0) || (tstep == mMaxStep) ) {
		
		// distribute forward
		distributeFwdWvfToElements();
		// distribute backward each time we compute the kernels 
		distributeBwdWvfToElements();

		int nuLine = 0;
		
		
		for (int ielem = 0; ielem < mKerElements.size(); ielem++) {


			KernerElement *kerElem = mKerElements[ielem];
			int nuMax = mKerElements[ielem]->getNuMax();
			
			if (tstep == mMaxStep) {
				kerElem->setBufferSize(mBufferSize);
			}
						
			kerElem->computeKernels();
			kerElem->feedKernels(mPhysicalKernels, nuLine, nuMax, mDumpTimeKernels, tstep);
			kerElem->clearKernels();
			nuLine += nuMax + 1;

		}
		if (mDumpTimeKernels && tstep != mMaxStep) dumpToFile();
	}

}

void Kerner::dumpToFile() {
	mIO->dumpToFile(mPhysicalKernels, mBufferSize);
}

void Kerner::distributeFwdWvfToElements() {

	vec_vec_ar6_RMatPP forward_disp;

	mIO->loadWavefield(forward_disp, mBufferSize);
	
	int nuOffset = 0;
	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {
		
		int nuFwd = mKerElements[ielem]->getNuForward();
		int nuMax = mKerElements[ielem]->getNuMax(); 
	
		vec_ar3_CMatPP initDispElem(nuMax + 1, zero_ar3_CMatPP);
		vec_vec_ar3_CMatPP dispElem(mBufferSize, initDispElem); 

		for (int it = 0; it < mBufferSize; it++) {

			for (int inu = 0; inu <= nuFwd; inu++) {
				
				dispElem[it][inu][0] = forward_disp[it][nuOffset + inu][0] + ii * forward_disp[it][nuOffset + inu][1];
				dispElem[it][inu][1] = forward_disp[it][nuOffset + inu][2] + ii * forward_disp[it][nuOffset + inu][3];
				dispElem[it][inu][2] = forward_disp[it][nuOffset + inu][4] + ii * forward_disp[it][nuOffset + inu][5];

			}					
		}
		
		nuOffset += nuFwd + 1;	
		mKerElements[ielem]->setForwardDisp(dispElem);
	

	}
	
	
}

void Kerner::distributeBwdWvfToElements() {
	
	int nuOffset = 0;

	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {
		
		KernerElement *kerElem = mKerElements[ielem];
		int nuBwd = kerElem->getNuBackward();
		int nuMax = kerElem->getNuMax();


		vec_ar3_CMatPP initDispElem(nuMax + 1, zero_ar3_CMatPP); 
		vec_vec_ar3_CMatPP dispElem(mBufferSize, initDispElem); 
		
		for (int it = 0; it < mBufferSize; it++) {
			
			for (int inu = 0; inu <= nuBwd; inu ++) {
				
				dispElem[it][inu][0] = mDomainRecorder->mBufferDisp[it][nuOffset + inu][0] + ii * mDomainRecorder->mBufferDisp[it][nuOffset + inu][1];
				dispElem[it][inu][1] = mDomainRecorder->mBufferDisp[it][nuOffset + inu][2] + ii * mDomainRecorder->mBufferDisp[it][nuOffset + inu][3];
				dispElem[it][inu][2] = mDomainRecorder->mBufferDisp[it][nuOffset + inu][4] + ii * mDomainRecorder->mBufferDisp[it][nuOffset + inu][5];

			}
			
		}
		
		nuOffset+= nuBwd + 1;	
		kerElem->setBackwardDisp(dispElem);


	}

}

void Kerner::distributeMaterialToElements() {

	vec_ar12_RMatPP materials;	//order is real and imag of rho, vph, vpv, vsh, vsv, eta.
	mIO->loadMaterial(materials);
	
	int nuOffset = 0;

	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {
		
		int nuFwd = mKerElements[ielem]->getNuForward();
		int nuMax = mKerElements[ielem]->getNuMax(); 
		
		vec_ar6_CMatPP materialsElem(nuMax + 1, zero_ar6_CMatPP);
		
		for (int inu = 0; inu <= nuFwd; inu ++) {
			
			materialsElem[inu][0] = materials[nuOffset + inu][0] + ii * materials[nuOffset + inu][1];
			materialsElem[inu][1] = materials[nuOffset + inu][2] + ii * materials[nuOffset + inu][3];
			materialsElem[inu][2] = materials[nuOffset + inu][4] + ii * materials[nuOffset + inu][5];
			materialsElem[inu][3] = materials[nuOffset + inu][6] + ii * materials[nuOffset + inu][7];
			materialsElem[inu][4] = materials[nuOffset + inu][8] + ii * materials[nuOffset + inu][9];
			materialsElem[inu][5] = materials[nuOffset + inu][10] + ii * materials[nuOffset + inu][11];

			
		}
		
		mKerElements[ielem]->setMaterials(materialsElem);
		nuOffset+=nuFwd + 1;	

	}
	
}


void Kerner::distributeNus() {
	
	std::vector<int> NusFwd;
	std::vector<int> NrsFwd;
	
	mIO->loadNus(NusFwd);
	mIO->loadNrs(NrsFwd);
	
	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {
		
		int nuBwd = mKerElements[ielem]->getNuBackward();
		int nuFwd = NusFwd[ielem]-1;
		int nuMax = nuBwd > nuFwd ? nuBwd : nuFwd; // we use max nu between fwd and bwd 
		int nrBwd = mKerElements[ielem]->getNrBackward();
		int nrFwd = NrsFwd[ielem];
		int nrMax = nrBwd > nrFwd ? nrBwd : nrFwd; // we use max nr between fwd and bwd 
		
		mKerElements[ielem]->setNuForward(nuFwd);
		mKerElements[ielem]->setNrForward(nrFwd);
		mKerElements[ielem]->setNuMax(nuMax);
		mKerElements[ielem]->setNrMax(nrMax);
		
		
	}
	
}