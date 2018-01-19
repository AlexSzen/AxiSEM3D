

#include "Kerner.h"
#include "KernerElement.h"
#include "KernerIO.h"
#include "XMPI.h"
#include "DomainRecorder.h"
#include "eigenc.h"
#include "PreloopFFTW.h" //just for lucky number 
#include "Processor.h"

Kerner::Kerner(bool dumpTimeKernels, int numKernels, std::vector<std::string> kerTypes, int totSteps, int maxNr, const RMatX2 filtParams, Real begWin, Real endWin): 
mDumpTimeKernels(dumpTimeKernels), mNumKernels(numKernels), mKerTypes(kerTypes), mMaxNr(maxNr), mFiltParams(filtParams), mBegWin(begWin), mEndWin(endWin) {
	
	mIO = new KernerIO();
	mTotSteps = PreloopFFTW::nextLuckyNumber(2 * totSteps + 1, false); //for fft we need to padd the wavefields with 0
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
	
	mIO->initialize(mDumpTimeKernels, mFiltParams.rows(), startElem, countElem[XMPI::rank()], mDomainRecorder->mTotalRecordSteps);

	Processor::initialize(mTotSteps, mDomainRecorder->mBufferTime, mFiltParams, mBegWin, mEndWin);	
	
}

void Kerner::finalize() {
	
	mIO->finalize();
	Processor::finalize();
}

void Kerner::computeKernels() {
	
	distributeBwdWvfToElements();
	distributeFwdWvfToElements();

	distributeMaterialToElements();
		
	// gather nus
	int totNu = 0; 
	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {
		totNu += mKerElements[ielem]->getNuForward() + 1;
	}
			
	//init kernel 
	vec_ar12_RMatPP initKernels(totNu, zero_ar12_RMatPP);
	//mPhysicalKernels.assign(Processor::sNumFilters, initKernels);
	mPhysicalKernels.assign(mTotSteps, initKernels);

	int nuLine = 0;
	
	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {
		
		KernerElement *kerElem = mKerElements[ielem];
		int nuElem = kerElem->getNuForward() + 1;

		kerElem->computeKernels(mDumpTimeKernels);
		kerElem->feedKernels(mPhysicalKernels, nuLine, nuElem);
		kerElem->clearKernels();
		nuLine += nuElem;
	}
	mIO->dumpToFile(mPhysicalKernels, Processor::sNumFilters);

}


void Kerner::distributeFwdWvfToElements() {
	
	vec_vec_ar6_RMatPP forward_disp;
	std::vector<int> Nus;
	std::vector<int> Nrs;
	int totSteps = mDomainRecorder->mTotalRecordSteps; //this is tot steps of the wvf, then we padd for fft  
	
	
	mIO->loadWavefield(forward_disp, Nus, Nrs);
	int nuOffset = 0;

	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {
		
		
		vec_ar3_CMatPP initDispElem(Nus[ielem], zero_ar3_CMatPP);
		vec_vec_ar3_CMatPP dispElem(mTotSteps, initDispElem); //here we pad for fft

		for (int it = 0; it < totSteps; it++) {
			
			
			for (int inu = 0; inu < Nus[ielem]; inu++) {
				
				dispElem[it][inu][0] = forward_disp[it][nuOffset + inu][0] + ii * forward_disp[it][nuOffset + inu][1];
				dispElem[it][inu][1] = forward_disp[it][nuOffset + inu][2] + ii * forward_disp[it][nuOffset + inu][3];
				dispElem[it][inu][2] = forward_disp[it][nuOffset + inu][4] + ii * forward_disp[it][nuOffset + inu][5];

			}					
		}
		
		nuOffset+=Nus[ielem];	
		
		// tapering 	
		Processor::taper(dispElem);
		
		// set displacement 
		mKerElements[ielem]->setForwardDisp(dispElem);
		mKerElements[ielem]->setNuForward(Nus[ielem]-1);
		mKerElements[ielem]->setNrForward(Nrs[ielem]);
	}
	
	
}

void Kerner::distributeBwdWvfToElements() {
	
	int totSteps = mDomainRecorder->mTotalRecordSteps; // these are those of the wavefield
	
	
	std::vector<int> NusFwd;	
	mIO->loadNus(NusFwd);
	int nuOffset = 0;

	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {
		
		KernerElement *kerElem = mKerElements[ielem];
		vec_ar3_CMatPP initDispElem(NusFwd[ielem], zero_ar3_CMatPP); // we assume NuFwd>=NuBwd, so bwd disp is padded to NuFwd for FFTs sake
		vec_vec_ar3_CMatPP dispElem(mTotSteps, initDispElem); //here we pad for fft 
		
		for (int it = 0; it < totSteps; it++) {
			
			for (int inu = 0; inu <= kerElem->getNuBackward(); inu ++) {
				
				dispElem[it][inu][0] = mDomainRecorder->mBufferDisp[it][nuOffset + inu][0] + ii * mDomainRecorder->mBufferDisp[it][nuOffset + inu][1];
				dispElem[it][inu][1] = mDomainRecorder->mBufferDisp[it][nuOffset + inu][2] + ii * mDomainRecorder->mBufferDisp[it][nuOffset + inu][3];
				dispElem[it][inu][2] = mDomainRecorder->mBufferDisp[it][nuOffset + inu][4] + ii * mDomainRecorder->mBufferDisp[it][nuOffset + inu][5];

			}
			
		}
		
		nuOffset+=kerElem->getNuBackward() + 1;	
		
		Processor::taper(dispElem);
		
		kerElem->setBackwardDisp(dispElem);
		kerElem->setTimeAndFreqSize(mTotSteps);

	}
		
	mDomainRecorder->mBufferDisp.clear(); //it will already be dumped if needed, so we clear it to free up RAM 
	
}

void Kerner::distributeMaterialToElements() {

	std::vector<int> NusFwd;
	vec_ar12_RMatPP materials;	//order is real and imag of rho, vph, vpv, vsh, vsv, eta.
	mIO->loadMaterial(materials, NusFwd);
	
	int nuOffset = 0;

	for (int ielem = 0; ielem < mKerElements.size(); ielem++) {
		
		vec_ar6_CMatPP materialsElem(NusFwd[ielem], zero_ar6_CMatPP);
		
		for (int inu = 0; inu < NusFwd[ielem]; inu ++) {
			
			materialsElem[inu][0] = materials[nuOffset + inu][0] + ii * materials[nuOffset + inu][1];
			materialsElem[inu][1] = materials[nuOffset + inu][2] + ii * materials[nuOffset + inu][3];
			materialsElem[inu][2] = materials[nuOffset + inu][4] + ii * materials[nuOffset + inu][5];
			materialsElem[inu][3] = materials[nuOffset + inu][6] + ii * materials[nuOffset + inu][7];
			materialsElem[inu][4] = materials[nuOffset + inu][8] + ii * materials[nuOffset + inu][9];
			materialsElem[inu][5] = materials[nuOffset + inu][10] + ii * materials[nuOffset + inu][11];

			
		}
		
		mKerElements[ielem]->setMaterials(materialsElem);
		nuOffset+=NusFwd[ielem];	

	}
	
}