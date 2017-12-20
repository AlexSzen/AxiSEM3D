// Acoustic1D.h
// created by Kuangdai on 23-Apr-2016 
// 1D acoustic
// In fluid it is potential formulation, so here strainToStress is actually potential to displacement.
// It is kept as strainToStress for homogeneity with other classes.
// However, in strainToStress the resulting displacement is premultiplied with integral factor to 
// then compute stifness.
// Hence potToDisp which doesnt premultiply and thus si;ply gives displacement. 
#pragma once

#include "Acoustic.h"
#include "eigenc.h"

class Acoustic1D: public Acoustic {
public:
    // constructor
    Acoustic1D(const RMatPP &KFluid, const RMatPP &KFluid_noquad): mKStruct(KFluid), mKStructNoQuad(KFluid_noquad) {};
    
    // STEP 2: strain ==> stress
    void strainToStress(FluidResponse &response) const;
	
	// computes displacement from potential in fluid.
	void potToDisp(FluidResponse &response) const;

    
    // verbose
    std::string verbose() const {return "Acoustic1D";};
    
    // 1D or Fourier space
    bool is1D() const {return true;};
    
private:
    RMatPP mKStruct; 
	RMatPP mKStructNoQuad; 

};
