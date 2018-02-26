
#include "KernerElement.h"
#include "FieldFFT.h"
#include "SolverFFTW_N3.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"
#include "Gradient.h"
#include "Processor.h"
#include <iostream>
KernerElement::KernerElement(const Element* elem): mElement(elem) {
	
	
}

void KernerElement::computeKernels(bool dumpTimeKernels) {
	
	
	///////// workspace 
	vec_ar9_CMatPP strainFwdF_it(mNuForward + 1, zero_ar9_CMatPP); //fourier strain for one timestep. defined here to avoid reallocation 
	vec_ar9_CMatPP strainBwdF_it(mNuForward + 1, zero_ar9_CMatPP); //fourier strain for one timestep. defined here to avoid reallocation 

	vec_ar9_RMatPP strainFwdR_it(mNrForward, zero_ar9_RMatPP); //physical strain for one timestep. defined here to avoid reallocation 
	vec_ar9_RMatPP strainBwdR_it(mNrForward, zero_ar9_RMatPP); //physical strain for one timestep. defined here to avoid reallocation 
	vec_ar3_RMatPP dispFwdR_it(mNrForward, zero_ar3_RMatPP);
	vec_ar3_RMatPP dispBwdR_it(mNrForward, zero_ar3_RMatPP);
	
	vec_ar9_RMatPP initstrainRT(mTimeSize, zero_ar9_RMatPP); 
	vec_ar9_RMatPP initstrainTR(mNrForward, zero_ar9_RMatPP); 
	vec_vec_ar9_RMatPP strainForwardRT(mNrForward, initstrainRT); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainForwardTR(mTimeSize, initstrainTR); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainBackwardRT(mNrForward, initstrainRT); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainBackwardTR(mTimeSize, initstrainTR); // all time steps for each physical slice 
	
	vec_ar3_RMatPP initdispRT(mTimeSize, zero_ar3_RMatPP); 
	vec_vec_ar3_RMatPP dispForwardRT(mNrForward, initdispRT);
	vec_vec_ar3_RMatPP dispBackwardRT(mNrForward, initdispRT);

	
	vec_ar9_CMatPP strainFwdFr_inu(mFreqSize, zero_ar9_CMatPP); //frequency strain for one nr. defined here to avoid reallocation 
	vec_ar9_CMatPP strainBwdFr_inu(mFreqSize, zero_ar9_CMatPP);  
	vec_ar3_CMatPP dispFwdFr_inu(mFreqSize, zero_ar3_CMatPP); 
	vec_ar3_CMatPP dispBwdFr_inu(mFreqSize, zero_ar3_CMatPP); 
	vec_ar3_CMatPP velFwdFr_inu(mFreqSize, zero_ar3_CMatPP); 
	vec_ar3_CMatPP velBwdFr_inu(mFreqSize, zero_ar3_CMatPP); 
	
	vec_ar9_CMatPP baseKernelsFr_inu(mFreqSize, zero_ar9_CMatPP);
	vec_ar9_CMatPP baseKernelsF_it(mNuForward + 1, zero_ar9_CMatPP);
	
	vec_ar6_RMatPP materialsR(mNrForward, zero_ar6_RMatPP);
	
	mBaseKernels.assign(mTimeSize,strainFwdF_it);

	//test material properties
	/*for (int it = 0; it<mTimeSize;it++)
		for (int inu = 0; inu<=mNuForward;inu++)
			mBaseKernels[it][inu][0] = mMaterials[inu][2];*/

	// transform material to real domain 
	FieldFFT::transformF2P(mMaterials, mNrForward);
	RMatXN6 unstructured_mat = SolverFFTW_N6::getC2R_RMat(mNrForward);
	FieldFFT::makeStruct<vec_ar6_RMatPP, RMatXN6>(materialsR, unstructured_mat, mNrForward - 1);
	
	///// get strain and go to physical azimuthal domain 
	for (int it = 0; it < mTimeSize; it ++) {
		// disp to strain 	
		mElement->mGradient->computeGrad9(mForwardDisp[it], strainFwdF_it, mNuForward, mNyquistFwd); 
		mElement->mGradient->computeGrad9(mBackwardDisp[it], strainBwdF_it, mNuForward, mNyquistFwd); 

		// fft F2P
		FieldFFT::transformF2P(strainFwdF_it, mNrForward);
		RMatXN9 unstructured = SolverFFTW_N9::getC2R_RMat(mNrForward);
		FieldFFT::makeStruct<vec_ar9_RMatPP, RMatXN9>(strainFwdR_it, unstructured, mNrForward - 1);

		FieldFFT::transformF2P(strainBwdF_it, mNrForward);
		unstructured = SolverFFTW_N9::getC2R_RMat(mNrForward);
		FieldFFT::makeStruct<vec_ar9_RMatPP, RMatXN9>(strainBwdR_it, unstructured, mNrForward - 1);
		
		// fft F2P
		FieldFFT::transformF2P(mForwardDisp[it], mNrForward);
		RMatXN3 unstructured_disp = SolverFFTW_N3::getC2R_RMat(mNrForward);
		FieldFFT::makeStruct<vec_ar3_RMatPP, RMatXN3>(dispFwdR_it, unstructured_disp, mNrForward - 1);

		FieldFFT::transformF2P(mBackwardDisp[it], mNrForward);
		unstructured_disp = SolverFFTW_N3::getC2R_RMat(mNrForward);
		FieldFFT::makeStruct<vec_ar3_RMatPP, RMatXN3>(dispBwdR_it, unstructured_disp, mNrForward - 1);
		
		for (int inu = 0; inu < mNrForward; inu ++) {
			strainForwardRT[inu][it] = strainFwdR_it[inu];
			strainBackwardRT[inu][it] = strainBwdR_it[inu];
			dispForwardRT[inu][it] = dispFwdR_it[inu];
			dispBackwardRT[inu][it] = dispBwdR_it[inu]; 
			
		}


	}
	
	
	/////// for each azimuthal slice, perform fft time to freq. 
	/////// compute base then physical kernels, and go back to time domain 
	for (int inu = 0; inu < mNrForward; inu++) {
		
		//fft T2F
		Processor::transformT2F(strainForwardRT[inu], strainFwdFr_inu);
		Processor::transformT2F(strainBackwardRT[inu], strainBwdFr_inu);
		Processor::transformT2F(dispForwardRT[inu], velFwdFr_inu);
		Processor::transformT2F(dispBackwardRT[inu], velBwdFr_inu);
		
		//get velocity 
		Processor::derivate(velFwdFr_inu);
		Processor::derivate(velBwdFr_inu);

		
		//convolve. formulas from Ficthner p.166
		
		for (int ifr = 0; ifr<mFreqSize; ifr++) {
				
				// need to multiply kernels by number time steps before FFT f2t, 
				// because each wavefield has been divided during FFT t2f but kernels would only be multiplied once
				// K_rho0
				baseKernelsFr_inu[ifr][0] = mTimeSize * (- velFwdFr_inu[ifr][0].schur(velBwdFr_inu[ifr][0]) - 
											velFwdFr_inu[ifr][1].schur(velBwdFr_inu[ifr][1]) - 
											velFwdFr_inu[ifr][2].schur(velBwdFr_inu[ifr][2]));
			
				
				//K_lambda0
				baseKernelsFr_inu[ifr][1] = mTimeSize * ((strainFwdFr_inu[ifr][0] + strainFwdFr_inu[ifr][4] + strainFwdFr_inu[ifr][8]).schur((strainBwdFr_inu[ifr][0] +
											strainBwdFr_inu[ifr][4] + strainBwdFr_inu[ifr][8])));
				
				//K_mu0
				baseKernelsFr_inu[ifr][2] = mTimeSize * (strainFwdFr_inu[ifr][0].schur(strainBwdFr_inu[ifr][0]) +
											strainFwdFr_inu[ifr][1].schur(strainBwdFr_inu[ifr][1]) +
											strainFwdFr_inu[ifr][2].schur(strainBwdFr_inu[ifr][2]) +
											strainFwdFr_inu[ifr][3].schur(strainBwdFr_inu[ifr][3]) +
											strainFwdFr_inu[ifr][4].schur(strainBwdFr_inu[ifr][4]) +
											strainFwdFr_inu[ifr][5].schur(strainBwdFr_inu[ifr][5]) +
											strainFwdFr_inu[ifr][6].schur(strainBwdFr_inu[ifr][6]) +
											strainFwdFr_inu[ifr][7].schur(strainBwdFr_inu[ifr][7]) +
											strainFwdFr_inu[ifr][8].schur(strainBwdFr_inu[ifr][8]));
											
		
		
		}
		// for now only one filter.
		// eventually will loop through filters and time integrate kernels 
	//	Processor::filter(baseKernelsFr_inu, 0);
		
		 
		//reuse strainForwardRT and strainBackwardRT for kernels to avoid extra allocation
		Processor::transformF2T(baseKernelsFr_inu, strainForwardRT[inu]);
		
		// compute physical kernels. formulas Fichtner p.167
		
		// before going back to Fourier azimuth, need to divide by Nr,
		// because both wavefields have been multiplied during FFT but kernels would only be divided once 
		Real inv_nr = one / (Real) mNrForward;
		for (int it = 0; it < mTimeSize; it ++) {
			
			// K_rho
			strainBackwardTR[it][inu][0] = inv_nr * (strainForwardRT[inu][it][0] + (materialsR[inu][4].schur(materialsR[inu][4])).schur(strainForwardRT[inu][it][2]) +
			(materialsR[inu][2].schur(materialsR[inu][2]) - two * materialsR[inu][4].schur(materialsR[inu][4])).schur(strainForwardRT[inu][it][1]));
			
			//K_vp
			strainBackwardTR[it][inu][1] = inv_nr * (two * strainForwardRT[inu][it][1].schur(materialsR[inu][0].schur(materialsR[inu][2])));
			
			//K_vs
			strainBackwardTR[it][inu][2] = inv_nr * (two * strainForwardRT[inu][it][2].schur(materialsR[inu][0].schur(materialsR[inu][4])) - 
											(Real)4.0 * strainForwardRT[inu][it][1].schur(materialsR[inu][0].schur(materialsR[inu][4])));

		}
		 
	}
	
	/////// go back to fourier azimuthal domain.
	for (int it = 0; it < mTimeSize; it ++) {
		
		//fft P2F
		RMatXN9 &unstruc = SolverFFTW_N9::getR2C_RMat(mNrForward);
		FieldFFT::makeFlat<vec_ar9_RMatPP, RMatXN9>(strainBackwardTR[it], unstruc, mNrForward - 1);
		FieldFFT::transformP2F(mBaseKernels[it], mNrForward);

	}
	
	/////// clear member variables to free up RAM.
	mForwardDisp.clear();
	mBackwardDisp.clear();
	mMaterials.clear();


	
	
}

void KernerElement::computeKernels2() {
	
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
	Real inv_nr = one / (Real) mNrMax;
	for (int it = 0; it < mBufferSize; it ++) {
		// disp to strain 	
		mElement->mGradient->computeGrad9(mForwardDisp[mTimeSize - 1 - it - mTimeLine], strainFwdF_it, mNuMax, mNyquistFwd); //time reverse
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
			baseKernels[it][inu][2] = two * inv_nr * (strainFwdR_it[inu][0] + strainFwdR_it[inu][4] + strainFwdR_it[inu][8]).schur((strainBwdR_it[inu][0] +
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
	mTimeLine += mBufferSize;

	
}



void KernerElement::feedKernels(vec_vec_ar6_RMatPP &physKernels, int nuLine, int nuElem, bool dumpTimeKernels) {

	if (dumpTimeKernels) {
		for (int it = 0; it < mBufferSize; it ++) {
			for (int inu = 0; inu < nuElem; inu ++) {
				
				physKernels[it][nuLine + inu][0] = mBaseKernels[it][inu][0].real();
				physKernels[it][nuLine + inu][1] = mBaseKernels[it][inu][0].imag();
				physKernels[it][nuLine + inu][2] = mBaseKernels[it][inu][1].real();
				physKernels[it][nuLine + inu][3] = mBaseKernels[it][inu][1].imag();
				physKernels[it][nuLine + inu][4] = mBaseKernels[it][inu][2].real();
				physKernels[it][nuLine + inu][5] = mBaseKernels[it][inu][2].imag();

				
			}
			
		}
	} else {
		for (int it = 0; it < mBufferSize; it ++) {
			for (int inu = 0; inu < nuElem; inu ++) {
				
				physKernels[0][nuLine + inu][0] += mBaseKernels[it][inu][0].real();
				physKernels[0][nuLine + inu][1] += mBaseKernels[it][inu][0].imag();
				physKernels[0][nuLine + inu][2] += mBaseKernels[it][inu][1].real();
				physKernels[0][nuLine + inu][3] += mBaseKernels[it][inu][1].imag();
				physKernels[0][nuLine + inu][4] += mBaseKernels[it][inu][2].real();
				physKernels[0][nuLine + inu][5] += mBaseKernels[it][inu][2].imag();

				
			}
			
		}		
	}
}

