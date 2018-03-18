// OffAxisPointForce.h
// created by Kuangdai on 11-Nov-2017
// off-axis point-force source

#include "OffAxisPointForce.h"
#include "Quad.h"
#include "SpectralConstants.h"
#include "XMath.h"
#include <sstream>
#include "Source.h"
#include "Relabelling.h"
#include <iostream>
OffAxisPointForce::OffAxisPointForce(double depth, double lat, double lon,
	double srcLat, double srcLon, double srcDep): Source(depth, lat, lon, srcLat, srcLon, srcDep)
 {
    // nothing
}

void OffAxisPointForce::computeSourceFourier(const Quad &myQuad, 
	const RDColP &interpFactZ, 
	const RDColP &interpFactXii,
	const RDColP &interpFactEta,
	double phi,
	arPP_CMatX3 &fouriers) const {

	// Fourier order
	int nu = myQuad.getNu();
    // set zero
	for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        fouriers[ipnt] = CMatX3::Zero(nu + 1, 3);
    }

    // particle relabelling
    RDRowN JPRT;
    if (myQuad.hasRelabelling()) {
        const RDMatXN &JJ = myQuad.getRelabelling().getStiffJacobian();
		JPRT = XMath::computeFourierAtPhi(JJ, phi);
    } else {
		JPRT = RDRowN::Ones();
	}
	//double eps = 1e-9; //we want the gaussian(numax) = eps 
	double a = 0.02;// sqrt( - (4./pow(nu,2.)) * log(eps) ); // a = variance / sqrt(2) 
    // compute source pointwise
	for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            double fact = interpFactXii(ipol) * interpFactEta(jpol);
			for (int beta = 0; beta <= nu; beta++) {
				double gauss_fact = exp( - pow(a, 2.) * pow(beta, 2.) / 4); //gaussian approximate delta in freq domain    
				for (int idim = 0; idim < 3; idim++) {
					fouriers[ipnt](beta, 2) = Complex( 
						(1. / (2. * pi)) *  fact * gauss_fact * exp(beta * phi * iid) * JPRT(ipnt));
					
				}
			}
		}
    }
	

	
}

std::string OffAxisPointForce::verbose() const {
	std::stringstream ss;
    ss << "\n================ Source ================" << std::endl;
    ss << "  Type         =   " << "Adjoint Source" << std::endl;
    ss << "  Latitude     =   " << mLatitude << std::endl;
    ss << "  Longitude    =   " << mLongitude << std::endl;
    ss << "  Depth (km)   =   " << mDepth / 1e3 << std::endl;
    ss << "================ Source ================\n" << std::endl;
	return ss.str();
}
