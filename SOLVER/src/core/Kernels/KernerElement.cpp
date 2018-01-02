
#include "KernerElement.h"
#include "FieldFFT.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"
#include "Gradient.h"
#include "Processor.h"
#include <iostream>
KernerElement::KernerElement(const Element* elem): mElement(elem) {
	
	
}

void KernerElement::computeKernels(bool dumpTimeKernels) {
	
	
	///////// workspace 	
	vec_ar6_CMatPP baseKernelsFr_inu(mFreqSize, zero_ar6_CMatPP);
	vec_ar6_CMatPP baseKernelsF_it(mNuForward + 1, zero_ar6_CMatPP);
	vec_ar6_RMatPP baseKernelsT_inu(mTimeSize, zero_ar6_RMatPP);
	vec_ar6_RMatPP baseKernelsR_it(mNrForward, zero_ar6_RMatPP);	
	
	vec_vec_ar6_CMatPP baseKernelsRFr(mNrForward, baseKernelsFr_inu);
	vec_vec_ar6_CMatPP baseKernelsTF(mTimeSize, baseKernelsF_it);
	vec_vec_ar6_RMatPP baseKernelsTR(mTimeSize, baseKernelsR_it);
	
	mBaseKernels.assign(Processor::sNumFilters, baseKernelsF_it); // integrated kernels for all filters 
	mPhysicalKernels.assign(Processor::sNumFilters, baseKernelsF_it);
	
	
	// compute frequency domain kernels
	computeBaseKernelsRFr(baseKernelsRFr);
	
	
	// for each filter, compute time domain kernels, integrate in time, and compute physical kernels 
	for (int ifilt = 0; ifilt < Processor::sNumFilters; ifilt++) {
	
		computeBaseKernelsF(baseKernelsRFr, baseKernelsTF, baseKernelsTR, baseKernelsT_inu, baseKernelsR_it, ifilt);
		computePhysicalKernels();
	
	}
	
	// clear member variables as they are now unused
	mMaterials.clear(); 
	mForwardDisp.clear();
	mBackwardDisp.clear();
	
}



void KernerElement::computeBaseKernelsF(vec_vec_ar6_CMatPP &kerRFr, vec_vec_ar6_CMatPP &kerTF, vec_vec_ar6_RMatPP &kerTR, vec_ar6_RMatPP &baseKernelsT_inu, vec_ar6_RMatPP &baseKernelsR_it, int ifilt) {
	
	
	for (int inu = 0; inu < mNrForward; inu++) {
		
		//filter 			
		Processor::filter<vec_ar6_CMatPP>(kerRFr[inu], ifilt);
		
		// fft F2T
		Processor::transformF2T(kerRFr[inu], baseKernelsT_inu);
		
		for (int it = 0; it < mTimeSize; it++) {
			
			kerTR[it][inu] = baseKernelsT_inu[inu];
			
		}
	}
	
	// fft P2F
	for (int it = 0; it < mTimeSize; it ++) {
		
		RMatXN6 &unstructured = SolverFFTW_N6::getR2C_RMat(mNrForward);
		FieldFFT::makeFlat<vec_ar6_RMatPP, RMatXN6>(kerTR[it], unstructured, mNrForward - 1);
		FieldFFT::transformP2F(kerTF[it], mNrForward);
		
	}
	
	Processor::timeWindow<vec_vec_ar6_CMatPP, vec_ar6_CMatPP>(kerTF, mBaseKernels[ifilt]);
	
			
	
}

void KernerElement::computePhysicalKernels() {
	
	// following Fichtner p. 169
	
}

void KernerElement::computeBaseKernelsRFr(vec_vec_ar6_CMatPP &kerRFr) {
	
	///////// workspace 
	vec_ar9_CMatPP strainFwdF_it(mNuForward + 1, zero_ar9_CMatPP); //fourier strain for one timestep. defined here to avoid reallocation 
	vec_ar9_CMatPP strainBwdF_it(mNuForward + 1, zero_ar9_CMatPP); //fourier strain for one timestep. defined here to avoid reallocation 

	vec_ar9_RMatPP strainFwdR_it(mNrForward, zero_ar9_RMatPP); //physical strain for one timestep. defined here to avoid reallocation 
	vec_ar9_RMatPP strainBwdR_it(mNrForward, zero_ar9_RMatPP); //physical strain for one timestep. defined here to avoid reallocation 

	vec_ar9_RMatPP initstrainRT(mTimeSize, zero_ar9_RMatPP); 
	vec_vec_ar9_RMatPP strainForwardRT(mNrForward, initstrainRT); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainBackwardRT(mNrForward, initstrainRT);
	
	vec_ar9_CMatPP strainFwdFr_inu(mFreqSize, zero_ar9_CMatPP); //frequency strain for one nr. defined here to avoid reallocation 
	vec_ar9_CMatPP strainBwdFr_inu(mFreqSize, zero_ar9_CMatPP); 
	
	vec_ar6_CMatPP baseKernelsFr_inu(mFreqSize, zero_ar6_CMatPP);

	
	///// get strain and go to physical domain 
	for (int it = 0; it < mTimeSize; it ++) {
		
		// disp to strain 	
		mElement->mGradient->computeGrad9(mForwardDisp[it], strainFwdF_it, mNuForward, mNyquist); 
		mElement->mGradient->computeGrad9(mBackwardDisp[it], strainBwdF_it, mNuForward, mNyquist); 
		
		// rotate from spz to rtp 
		//
		
		// fft F2P
		FieldFFT::transformF2P(strainFwdF_it, mNrForward);
		RMatXN9 unstructured = SolverFFTW_N9::getC2R_RMat(mNrForward);
		FieldFFT::makeStruct<vec_ar9_RMatPP, RMatXN9>(strainFwdR_it, unstructured, mNrForward - 1);
		
		FieldFFT::transformF2P(strainBwdF_it, mNrForward);
		unstructured = SolverFFTW_N9::getC2R_RMat(mNrForward);
		FieldFFT::makeStruct<vec_ar9_RMatPP, RMatXN9>(strainBwdR_it, unstructured, mNrForward - 1);
		
		for (int inu = 0; inu < mNrForward; inu++) {
			strainBackwardRT[inu][it] = strainBwdR_it[inu];
			strainForwardRT[inu][it] = strainFwdR_it[inu];
		}
					
	}
	
	////// go to frequency domain and perform convolution for each slice 
	for (int inu = 0; inu < mNrForward; inu++) {
		
		//fft T2F
		Processor::transformT2F(strainForwardRT[inu], strainFwdFr_inu);
		Processor::transformT2F(strainBackwardRT[inu], strainBwdFr_inu);

		/////////////////// convolve to get base kernels, following Fichtner p.169 
		
		///// Rho 
		// not implemented yet 
		
		///// Lambda 
		IColX indsFwdLambda(3);
		IColX indsBwdLambda(3);
		indsFwdLambda << 0, 4, 8; //trace 
		indsBwdLambda << 0, 4, 8;
		Real prefLambda = 1.;
		Processor::sumAndConvolve<vec_ar9_CMatPP>(strainFwdFr_inu, strainBwdFr_inu, baseKernelsFr_inu, indsFwdLambda, indsBwdLambda, prefLambda , 1);
		
		///// Mu 
		IColX indsFwdMu(1);
		IColX indsBwdMu(1);
		Real prefMu = 2.;
		for (int ic = 0; ic < 9 /* tensor product */; ic++) {
			indsFwdMu(0) = ic;
			indsBwdMu(0) = ic;
			Processor::sumAndConvolve<vec_ar9_CMatPP>(strainFwdFr_inu, strainBwdFr_inu, baseKernelsFr_inu, indsFwdMu, indsBwdMu, prefMu, 2);
		}
		
		///// a 
		IColX indsFwdA(2);
		IColX indsBwdA(2);
		indsFwdA << 4, 8; //tt + pp 
		indsBwdA << 4, 8;
		Real prefA = 1.;
		Processor::sumAndConvolve<vec_ar9_CMatPP>(strainFwdFr_inu, strainBwdFr_inu, baseKernelsFr_inu, indsFwdA, indsBwdA, prefA, 3);

		///// b 
		IColX indsFwdB(1);
		IColX indsBwdB(1);
		Real prefB = 4.;
		indsFwdB(0) = 3; // rt
		indsBwdB(0) = 3; 
		Processor::sumAndConvolve<vec_ar9_CMatPP>(strainFwdFr_inu, strainBwdFr_inu, baseKernelsFr_inu, indsFwdB, indsBwdB, prefB, 4);
		indsFwdB(0) = 2; // rp
		indsBwdB(0) = 2; 
		Processor::sumAndConvolve<vec_ar9_CMatPP>(strainFwdFr_inu, strainBwdFr_inu, baseKernelsFr_inu, indsFwdB, indsBwdB, prefB, 4);

		///// c 
		IColX indsFwdC(1);
		IColX indsBwdC(2);
		Real prefC = 1.;
		indsFwdC << 1; // rr
		indsBwdC << 4, 8; //pp + tt 
		Processor::sumAndConvolve<vec_ar9_CMatPP>(strainFwdFr_inu, strainBwdFr_inu, baseKernelsFr_inu, indsFwdC, indsBwdC, prefC, 5);
		indsFwdC.resize(2);
		indsBwdC.resize(1);
		indsFwdC << 4, 8;
		indsBwdC << 1;
		Processor::sumAndConvolve<vec_ar9_CMatPP>(strainFwdFr_inu, strainBwdFr_inu, baseKernelsFr_inu, indsFwdC, indsBwdC, prefC, 5);
		
		
		// assign to slice 
		kerRFr[inu] = baseKernelsFr_inu;
		
		
		
	}
	

	
}

void KernerElement::feedTimeKernelsBuffer() {
	
	
}