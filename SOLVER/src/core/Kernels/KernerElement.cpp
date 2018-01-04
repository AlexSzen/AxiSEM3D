
#include "KernerElement.h"
#include "FieldFFT.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"
#include "Gradient.h"
#include "Processor.h"
#include <iostream>
KernerElement::KernerElement(const Element* elem): mElement(elem) {
	
	
}

void KernerElement::test() {
	
	
	///////// workspace 
	vec_ar9_CMatPP strainFwdF_it(mNuForward + 1, zero_ar9_CMatPP); //fourier strain for one timestep. defined here to avoid reallocation 
	vec_ar9_CMatPP strainBwdF_it(mNuForward + 1, zero_ar9_CMatPP); //fourier strain for one timestep. defined here to avoid reallocation 

	vec_ar9_RMatPP strainFwdR_it(mNrForward, zero_ar9_RMatPP); //physical strain for one timestep. defined here to avoid reallocation 
	vec_ar9_RMatPP strainBwdR_it(mNrForward, zero_ar9_RMatPP); //physical strain for one timestep. defined here to avoid reallocation 

	vec_ar9_RMatPP initstrainRT(mTimeSize, zero_ar9_RMatPP); 
	vec_ar9_RMatPP initstrainTR(mNrForward, zero_ar9_RMatPP); 
	vec_vec_ar9_RMatPP strainForwardRT(mNrForward, initstrainRT); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainForwardTR(mTimeSize, initstrainTR); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainBackwardRT(mNrForward, initstrainRT); // all time steps for each physical slice 
	vec_vec_ar9_RMatPP strainBackwardTR(mTimeSize, initstrainTR); // all time steps for each physical slice 
	
	vec_ar9_CMatPP strainFwdFr_inu(mFreqSize, zero_ar9_CMatPP); //frequency strain for one nr. defined here to avoid reallocation 
	vec_ar9_CMatPP strainBwdFr_inu(mFreqSize, zero_ar9_CMatPP); //frequency strain for one nr. defined here to avoid reallocation 

	vec_ar9_CMatPP baseKernelsFr_inu(mFreqSize, zero_ar9_CMatPP);
	vec_ar9_CMatPP baseKernelsF_it(mNuForward + 1, zero_ar9_CMatPP);
	
	
	mBaseKernels.assign(mTimeSize,strainFwdF_it);
	
	///// get strain and go to physical domain 
	/*for (int it = 0; it < mTimeSize; it ++) {
		
		// disp to strain 	
		mElement->mGradient->computeGrad9(mForwardDisp[it], strainFwdF_it, mNuForward, mNyquist); 
		mElement->mGradient->computeGrad9(mBackwardDisp[it], strainBwdF_it, mNuForward, mNyquist); 

		// fft F2P
		FieldFFT::transformF2P(strainFwdF_it, mNrForward);
		RMatXN9 unstructured = SolverFFTW_N9::getC2R_RMat(mNrForward);
		FieldFFT::makeStruct<vec_ar9_RMatPP, RMatXN9>(strainFwdR_it, unstructured, mNrForward - 1);
		
		FieldFFT::transformF2P(strainBwdF_it, mNrForward);
		unstructured = SolverFFTW_N9::getC2R_RMat(mNrForward);
		FieldFFT::makeStruct<vec_ar9_RMatPP, RMatXN9>(strainBwdR_it, unstructured, mNrForward - 1);
		
		for (int inu = 0; inu < mNrForward; inu ++) {
			
			strainForwardRT[inu][it] = strainFwdR_it[inu];
			strainBackwardRT[inu][it] = strainBwdR_it[inu];

			
		}


	}
	
	for (int inu = 0; inu < mNrForward; inu++) {
		
		//fft T2F
		Processor::transformT2F(strainForwardRT[inu], strainFwdFr_inu);
		Processor::transformT2F(strainBackwardRT[inu], strainBwdFr_inu);
		//fft F2T
		//Processor::transformF2T(strainFwdFr_inu, strainForwardRT[inu]);
		
		//fft T2F

		//fft F2T
		//Processor::transformF2T(strainBwdFr_inu, strainBackwardRT[inu]);
		
		//convolve 
		for (int ifr = 0; ifr<mFreqSize; ifr++) {
				baseKernelsFr_inu[ifr][0] = (strainFwdFr_inu[ifr][0] + strainFwdFr_inu[ifr][4] + strainFwdFr_inu[ifr][8]).schur((strainBwdFr_inu[ifr][0] + strainBwdFr_inu[ifr][4] + strainBwdFr_inu[ifr][8]));
		}
		Processor::filter(baseKernelsFr_inu, 0);
		
		//fft F2T 
		Processor::transformF2T(baseKernelsFr_inu, strainBackwardRT[inu]);
		for (int it = 0; it < mTimeSize; it ++) {
			//strainForwardTR[it][inu] = strainForwardRT[inu][it];
			strainBackwardTR[it][inu] = strainBackwardRT[inu][it];
		}
		 
	}
	
	for (int it = 0; it < mTimeSize; it ++) {
		
		//fft P2F
	//	RMatXN9 &unstructured = SolverFFTW_N9::getR2C_RMat(mNrForward);
		//FieldFFT::makeFlat<vec_ar9_RMatPP, RMatXN9>(strainForwardTR[it], unstructured, mNrForward - 1);
		//FieldFFT::transformP2F(mBaseKernels[it], mNrForward);
		RMatXN9 &unstruc = SolverFFTW_N9::getR2C_RMat(mNrForward);
		FieldFFT::makeFlat<vec_ar9_RMatPP, RMatXN9>(strainBackwardTR[it], unstruc, mNrForward - 1);
		FieldFFT::transformP2F(mBaseKernels[it], mNrForward);
//		for (int inu = 0; inu <=mNuForward; inu++) {
//			mBaseKernels[it][inu][0] = mBaseKernels[it][inu][0] + mBaseKernels[it][inu][4] + mBaseKernels[it][inu][8];
//		}
	}
	*/
	
	//test material. 
	for (int inu = 0; inu < mNuForward +1 ;inu++)
		mBaseKernels[0][inu][0] = mMaterials[inu][1];

	
	
}

void KernerElement::computeKernels(bool dumpTimeKernels) {
	
	
	///////// workspace 	
	vec_ar6_CMatPP baseKernelsFr_inu(mFreqSize, zero_ar6_CMatPP);
	vec_ar6_CMatPP baseKernelsF_it(mNuForward + 1, zero_ar6_CMatPP);
	vec_ar6_RMatPP baseKernelsT_inu(mTimeSize, zero_ar6_RMatPP);
	vec_ar6_RMatPP baseKernelsR_it(mNrForward, zero_ar6_RMatPP);	
	
	vec_vec_ar6_CMatPP baseKernelsRFr(mNrForward, baseKernelsFr_inu);
//	vec_vec_ar6_CMatPP baseKernelsTF(mTimeSize, baseKernelsF_it);
//	mBaseKernels.assign(mTimeSize,baseKernelsF_it);
	vec_vec_ar6_RMatPP baseKernelsTR(mTimeSize, baseKernelsR_it);
	
	//mPhysicalKernels.assign(Processor::sNumFilters, baseKernelsF_it);
	
	
	// compute frequency domain kernels on slices 
	computeBaseKernelsRFr(baseKernelsRFr);
	
	
	// for each filter, compute time domain kernels, integrate in time, and compute physical kernels 
	//for (int ifilt = 0; ifilt < Processor::sNumFilters; ifilt++) {
	
	//computeBaseKernelsTF(baseKernelsRFr, mBaseKernels, baseKernelsTR, baseKernelsT_inu, baseKernelsR_it, 0);
//		computePhysicalKernels(ifilt);
	
//	}
	

	// clear member variables as they are now unused
	mMaterials.clear(); 
//	mForwardDisp.clear();
	//mBackwardDisp.clear();
	
}



void KernerElement::computeBaseKernelsTF(vec_vec_ar6_CMatPP &kerRFr, vec_vec_ar6_CMatPP &kerTF, vec_vec_ar6_RMatPP &kerTR, vec_ar6_RMatPP &baseKernelsT_inu, vec_ar6_RMatPP &baseKernelsR_it, int ifilt) {
	
	//mBaseKernels.assign(mNuForward + 1, zero_ar6_CMatPP);


	for (int inu = 0; inu < mNrForward; inu++) {
		
		//vec_ar6_RMatPP baseKernelsT_inu(mTimeSize, zero_ar6_RMatPP);
	//	vec_ar6_RMatPP baseKernelsR_it(mNrForward, zero_ar6_RMatPP);	
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
	
	//Processor::timeWindow<vec_vec_ar6_CMatPP, vec_ar6_CMatPP>(kerTF, mBaseKernels);

			
	
}

/*void KernerElement::computePhysicalKernels(int ifilt) {
	
	// following Fichtner p. 169
	
	// material numbering : rho, vph, vpv, vsh, vsv, eta
	// base kernel numbering : rho, lambda, mu, a, b, c
	// phyical kernels numbering : rho, vph, vpv, vsh, vsv, eta 
	
	vec_ar6_RMatPP materialsR(mNrForward, zero_ar6_RMatPP);
	vec_ar6_RMatPP physicalKernelsR(mNrForward, zero_ar6_RMatPP);
	vec_ar6_RMatPP baseKernelsR(mNrForward, zero_ar6_RMatPP);
	
	// fft F2P materials 
	FieldFFT::transformF2P(mMaterials, mNrForward);
	RMatXN6 unstructured = SolverFFTW_N6::getC2R_RMat(mNrForward);
	FieldFFT::makeStruct<vec_ar6_RMatPP, RMatXN6>(materialsR, unstructured, mNrForward - 1);
	
	FieldFFT::transformF2P(mBaseKernels, mNrForward);
	unstructured = SolverFFTW_N6::getC2R_RMat(mNrForward);
	FieldFFT::makeStruct<vec_ar6_RMatPP, RMatXN6>(baseKernelsR, unstructured, mNrForward - 1);
	
	/////// physical kernels 
	for (int inu = 0; inu < mNrForward; inu++) {
		
		// rho 
		// not implemented 
		
		// vph 
		physicalKernelsR[inu][1] = two * materialsR[inu][0] * materialsR[inu][1] * (baseKernelsR[inu][3] + materialsR[inu][5] * baseKernelsR[inu][5]);
		
		// vpv 
		physicalKernelsR[inu][2] = materialsR[inu][3];//two * materialsR[inu][0] * materialsR[inu][2] * (baseKernelsR[inu][1] - baseKernelsR[inu][3] - baseKernelsR[inu][5]);
		
		// vsh 
		physicalKernelsR[inu][3] = materialsR[inu][3];//two * materialsR[inu][3] * 
		//(baseKernelsR[inu][2] -  two * baseKernelsR[inu][1] - baseKernelsR[inu][4] + two * ( (one - materialsR[inu][5].array()).matrix()) * baseKernelsR[inu][4]);
		
		// vsv 
		physicalKernelsR[inu][4] = materialsR[inu][3];//two * materialsR[inu][0] * materialsR[inu][4] * baseKernelsR[inu][4];
		
		// eta 
		physicalKernelsR[inu][5] = materialsR[inu][3];//(materialsR[inu][1] + materialsR[inu][3]) * baseKernelsR[inu][5]; 
		
		
	}
	
	RMatXN6 &unstruct = SolverFFTW_N6::getR2C_RMat(mNrForward);
	FieldFFT::makeFlat<vec_ar6_RMatPP, RMatXN6>(physicalKernelsR, unstruct, mNrForward - 1);
	FieldFFT::transformP2F(mPhysicalKernels[ifilt], mNrForward);
	

	mBaseKernels.clear(); //we clear them so we can reuse them in next filter. 
	
	
}*/

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

void KernerElement::feedKernels(vec_vec_ar12_RMatPP &physKernels, int nuLine, int nuElem) {
	
	
	/*for (int ifilt = 0; ifilt < Processor::sNumFilters; ifilt ++) {
		for (int inu = 0; inu < nuElem; inu ++) {
			
			physKernels[ifilt][nuLine + inu][0] = mPhysicalKernels[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][1] = mPhysicalKernels[ifilt][inu][0].imag();
			physKernels[ifilt][nuLine + inu][2] = mPhysicalKernels[ifilt][inu][1].real();
			physKernels[ifilt][nuLine + inu][3] = mPhysicalKernels[ifilt][inu][1].imag();
			physKernels[ifilt][nuLine + inu][4] = mPhysicalKernels[ifilt][inu][2].real();
			physKernels[ifilt][nuLine + inu][5] = mPhysicalKernels[ifilt][inu][2].imag();
			physKernels[ifilt][nuLine + inu][6] = mPhysicalKernels[ifilt][inu][3].real();
			physKernels[ifilt][nuLine + inu][7] = mPhysicalKernels[ifilt][inu][3].imag();
			physKernels[ifilt][nuLine + inu][8] = mPhysicalKernels[ifilt][inu][4].real();
			physKernels[ifilt][nuLine + inu][9] = mPhysicalKernels[ifilt][inu][4].imag();
			physKernels[ifilt][nuLine + inu][10] = mPhysicalKernels[ifilt][inu][5].real();
			physKernels[ifilt][nuLine + inu][11] = mPhysicalKernels[ifilt][inu][5].imag();

			
		}
		
		
	}*/
	
	for (int ifilt = 0; ifilt < mTimeSize; ifilt ++) {
		for (int inu = 0; inu < nuElem; inu ++) {
			
			physKernels[ifilt][nuLine + inu][0] = mBaseKernels[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][1] = mBaseKernels[ifilt][inu][0].imag();
			physKernels[ifilt][nuLine + inu][2] = mBaseKernels[ifilt][inu][1].real();
			physKernels[ifilt][nuLine + inu][3] = mBaseKernels[ifilt][inu][1].imag();
			physKernels[ifilt][nuLine + inu][4] = mBaseKernels[ifilt][inu][2].real();
			physKernels[ifilt][nuLine + inu][5] = mBaseKernels[ifilt][inu][2].imag();
			physKernels[ifilt][nuLine + inu][6] = mBaseKernels[ifilt][inu][3].real();
			physKernels[ifilt][nuLine + inu][7] = mBaseKernels[ifilt][inu][3].imag();
			physKernels[ifilt][nuLine + inu][8] = mBaseKernels[ifilt][inu][4].real();
			physKernels[ifilt][nuLine + inu][9] = mBaseKernels[ifilt][inu][4].imag();
			physKernels[ifilt][nuLine + inu][10] = mBaseKernels[ifilt][inu][5].real();
			physKernels[ifilt][nuLine + inu][11] = mBaseKernels[ifilt][inu][5].imag();

			
		}
		
		
	}
	
	/*for (int ifilt = 0; ifilt < mTimeSize; ifilt ++) {
		for (int inu = 0; inu < nuElem; inu ++) {
			
			physKernels[ifilt][nuLine + inu][0] = mForwardDisp[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][1] = mForwardDisp[ifilt][inu][0].imag();
			physKernels[ifilt][nuLine + inu][2] = mForwardDisp[ifilt][inu][1].real();
			physKernels[ifilt][nuLine + inu][3] = mForwardDisp[ifilt][inu][1].imag();
			physKernels[ifilt][nuLine + inu][4] = mForwardDisp[ifilt][inu][2].real();
			physKernels[ifilt][nuLine + inu][5] = mForwardDisp[ifilt][inu][2].imag();
			physKernels[ifilt][nuLine + inu][6] = mForwardDisp[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][7] = mForwardDisp[ifilt][inu][0].imag();
			physKernels[ifilt][nuLine + inu][8] = mForwardDisp[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][9] = mForwardDisp[ifilt][inu][0].imag();
			physKernels[ifilt][nuLine + inu][10] = mForwardDisp[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][11] = mForwardDisp[ifilt][inu][0].imag();

			
		}
		
		
	}*/
	
	
/*	for (int ifilt = 0; ifilt < mTimeSize; ifilt ++) {
		for (int inu = 0; inu < nuElem; inu ++) {
			
			physKernels[ifilt][nuLine + inu][0] = mBackwardDisp[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][1] = mBackwardDisp[ifilt][inu][0].imag();
			physKernels[ifilt][nuLine + inu][2] = mBackwardDisp[ifilt][inu][1].real();
			physKernels[ifilt][nuLine + inu][3] = mBackwardDisp[ifilt][inu][1].imag();
			physKernels[ifilt][nuLine + inu][4] = mBackwardDisp[ifilt][inu][2].real();
			physKernels[ifilt][nuLine + inu][5] = mBackwardDisp[ifilt][inu][2].imag();
			physKernels[ifilt][nuLine + inu][6] = mBackwardDisp[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][7] = mBackwardDisp[ifilt][inu][0].imag();
			physKernels[ifilt][nuLine + inu][8] = mBackwardDisp[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][9] = mBackwardDisp[ifilt][inu][0].imag();
			physKernels[ifilt][nuLine + inu][10] = mBackwardDisp[ifilt][inu][0].real();
			physKernels[ifilt][nuLine + inu][11] = mBackwardDisp[ifilt][inu][0].imag();

			
		}
		
		
	}*/


}

