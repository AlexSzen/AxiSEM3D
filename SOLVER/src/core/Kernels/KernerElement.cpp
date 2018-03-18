
#include "KernerElement.h"
#include "FieldFFT.h"
#include "SolverFFTW_N3.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"
#include "Gradient.h"
#include "Processor.h"
#include <iostream>
KernerElement::KernerElement(const Element* elem, int bufferSize, Real deltaT): mElement(elem), mBufferSize(bufferSize), mDeltaT(deltaT) {
	
	
}

void KernerElement::computeKernels() {
	
	///////// workspace 
	vec_ar9_CMatPP strainFwdF_it(mNuMax + 1, zero_ar9_CMatPP); //fourier strain for one timestep. defined here to avoid reallocation 
	vec_ar9_CMatPP strainBwdF_it(mNuMax + 1, zero_ar9_CMatPP); //fourier strain for one timestep. defined here to avoid reallocation 
	
	vec_ar9_RMatPP strainFwdR_it(mNrMax, zero_ar9_RMatPP); //physical strain for one timestep. defined here to avoid reallocation 
	vec_ar9_RMatPP strainBwdR_it(mNrMax, zero_ar9_RMatPP); //physical strain for one timestep. defined here to avoid reallocation 
	vec_ar3_RMatPP dispFwdR_it(mNrMax, zero_ar3_RMatPP);
	vec_ar3_RMatPP dispBwdR_it(mNrMax, zero_ar3_RMatPP);
	
	vec_ar9_RMatPP initstrainRT(mBufferSize, zero_ar9_RMatPP); 
	vec_ar9_RMatPP initstrainTR(mNrMax, zero_ar9_RMatPP); 
	vec_vec_ar9_RMatPP strainForwardRT(mNrMax, initstrainRT); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainForwardTR(mBufferSize, initstrainTR); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainBackwardRT(mNrMax, initstrainRT); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainBackwardTR(mBufferSize, initstrainTR); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP baseKernels(mBufferSize, initstrainTR); // all time steps for each physical slice 

	
	vec_ar6_RMatPP materialsR(mNrMax, zero_ar6_RMatPP);
	
	mBaseKernels.assign(mBufferSize,strainFwdF_it);
	int nrBwd = getNrBackward();
	int nyquistBwd = (int) (nrBwd % 2 == 0);
	// transform material to real domain 
	FieldFFT::transformF2P(mMaterials, mNrMax);
	RMatXN6 unstructured_mat = SolverFFTW_N6::getC2R_RMat(mNrMax);
	FieldFFT::makeStruct<vec_ar6_RMatPP, RMatXN6>(materialsR, unstructured_mat, mNrMax - 1);
	
	///// get strain and go to physical azimuthal domain 
	for (int it = 0; it < mBufferSize; it ++) {
		// disp to strain 	
		mElement->mGradient->computeGrad9(mForwardDisp[mBufferSize - 1 - it], strainFwdF_it, mNuMax, mNyquistFwd); //time reverse
		mElement->mGradient->computeGrad9(mBackwardDisp[it], strainBwdF_it, mNuMax, nyquistBwd); 

		// fft F2P
		FieldFFT::transformF2P(strainFwdF_it, mNrMax);
		RMatXN9 unstructured = SolverFFTW_N9::getC2R_RMat(mNrMax);
		FieldFFT::makeStruct<vec_ar9_RMatPP, RMatXN9>(strainFwdR_it, unstructured, mNrMax - 1);

		FieldFFT::transformF2P(strainBwdF_it, mNrMax);
		unstructured = SolverFFTW_N9::getC2R_RMat(mNrMax);
		FieldFFT::makeStruct<vec_ar9_RMatPP, RMatXN9>(strainBwdR_it, unstructured, mNrMax - 1);
		
		// fft F2P
		FieldFFT::transformF2P(mForwardDisp[it], mNrMax);
		RMatXN3 unstructured_disp = SolverFFTW_N3::getC2R_RMat(mNrMax);
		FieldFFT::makeStruct<vec_ar3_RMatPP, RMatXN3>(dispFwdR_it, unstructured_disp, mNrMax - 1);

		FieldFFT::transformF2P(mBackwardDisp[it], mNrMax);
		unstructured_disp = SolverFFTW_N3::getC2R_RMat(mNrMax);
		FieldFFT::makeStruct<vec_ar3_RMatPP, RMatXN3>(dispBwdR_it, unstructured_disp, mNrMax - 1);
		
		for (int inu = 0; inu < mNrMax; inu ++) { // only trial vp 
			baseKernels[it][inu][2] = two *  (strainFwdR_it[inu][0] + strainFwdR_it[inu][4] + strainFwdR_it[inu][8]).schur((strainBwdR_it[inu][0] +
										strainBwdR_it[inu][4] + strainBwdR_it[inu][8])).schur(materialsR[inu][0]).schur(materialsR[inu][2]);
			
		}


	}

	/////// go back to fourier azimuthal domain.
	for (int it = 0; it < mBufferSize; it ++) {
		
		//fft P2F
		RMatXN9 &unstruc = SolverFFTW_N9::getR2C_RMat(mNrMax);
		FieldFFT::makeFlat<vec_ar9_RMatPP, RMatXN9>(baseKernels[it], unstruc, mNrMax - 1);
		FieldFFT::transformP2F(mBaseKernels[it], mNrMax);

	}
	
	/////// clear member variables to free up RAM.
	mBackwardDisp.clear();
	mForwardDisp.clear();

	
}



void KernerElement::feedKernels(vec_vec_ar6_RMatPP &physKernels, int nuLine, int nuElem, bool dumpTimeKernels, int tstep) { 

	if (dumpTimeKernels) { //movie
		
		if (mBufferKernels.size() == 0) {
			mBufferKernels.assign(nuElem + 1, zero_ar9_CMatPP);
		}
	
		for (int inu = 0; inu <= nuElem; inu ++) {
				
			physKernels[0][nuLine + inu][0] = mDeltaT * mBaseKernels[0][inu][0].real() + mBufferKernels[inu][0].real();
			physKernels[0][nuLine + inu][1] = mDeltaT * mBaseKernels[0][inu][0].imag() + mBufferKernels[inu][0].imag();
			physKernels[0][nuLine + inu][2] = mDeltaT * mBaseKernels[0][inu][1].real() + mBufferKernels[inu][1].real();
			physKernels[0][nuLine + inu][3] = mDeltaT * mBaseKernels[0][inu][1].imag() + mBufferKernels[inu][1].imag();
			physKernels[0][nuLine + inu][4] = mDeltaT * mBaseKernels[0][inu][2].real() + mBufferKernels[inu][2].real();
			physKernels[0][nuLine + inu][5] = mDeltaT * mBaseKernels[0][inu][2].imag() + mBufferKernels[inu][2].imag();
		
			for (int it = 1; it < mBufferSize; it ++) {
		
				
				physKernels[it][nuLine + inu][0] += physKernels[it-1][nuLine + inu][0] + mDeltaT * mBaseKernels[it][inu][0].real();
				physKernels[it][nuLine + inu][1] += physKernels[it-1][nuLine + inu][1] + mDeltaT * mBaseKernels[it][inu][0].imag();
				physKernels[it][nuLine + inu][2] += physKernels[it-1][nuLine + inu][2] + mDeltaT * mBaseKernels[it][inu][1].real();
				physKernels[it][nuLine + inu][3] += physKernels[it-1][nuLine + inu][3] + mDeltaT * mBaseKernels[it][inu][1].imag();
				physKernels[it][nuLine + inu][4] += physKernels[it-1][nuLine + inu][4] + mDeltaT * mBaseKernels[it][inu][2].real();
				physKernels[it][nuLine + inu][5] += physKernels[it-1][nuLine + inu][5] + mDeltaT * mBaseKernels[it][inu][2].imag();

				
			}
			
		
		
		mBufferKernels[inu][0].real() = physKernels[mBufferSize - 1][nuLine + inu][0];
		mBufferKernels[inu][0].imag() = physKernels[mBufferSize - 1][nuLine + inu][1];
		mBufferKernels[inu][1].real() = physKernels[mBufferSize - 1][nuLine + inu][2];
		mBufferKernels[inu][1].imag() = physKernels[mBufferSize - 1][nuLine + inu][3];
		mBufferKernels[inu][2].real() = physKernels[mBufferSize - 1][nuLine + inu][4];
		mBufferKernels[inu][2].imag() = physKernels[mBufferSize - 1][nuLine + inu][5];
	}

	} else {
		for (int it = 0; it < mBufferSize; it ++) {
			for (int inu = 0; inu <= nuElem; inu ++) {
				
				physKernels[0][nuLine + inu][0] += mDeltaT * mBaseKernels[it][inu][0].real();
				physKernels[0][nuLine + inu][1] += mDeltaT * mBaseKernels[it][inu][0].imag();
				physKernels[0][nuLine + inu][2] += mDeltaT * mBaseKernels[it][inu][1].real();
				physKernels[0][nuLine + inu][3] += mDeltaT * mBaseKernels[it][inu][1].imag();
				physKernels[0][nuLine + inu][4] += mDeltaT * mBaseKernels[it][inu][2].real();
				physKernels[0][nuLine + inu][5] += mDeltaT * mBaseKernels[it][inu][2].imag();

				
			}
			
		}		
	}
}

