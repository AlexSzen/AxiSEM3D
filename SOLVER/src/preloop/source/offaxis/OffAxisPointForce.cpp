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
    // compute source pointwise
	for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            double fact = interpFactXii(ipol) * interpFactEta(jpol);
			for (int beta = 0; beta <= nu; beta++) {
				for (int idim = 0; idim < 3; idim++) {
					if (idim == 2 && beta == 0) {// only monopole vertical force 
						fouriers[ipnt](beta, idim) = Complex(
							1e20 * fact * exp(beta * phi * iid) * JPRT(ipnt));
					}
				}
			}
		}
    }
	
}

std::string OffAxisPointForce::verbose() const {
	// TODO: should be verbosed collectively like receivers
	return "";
}
