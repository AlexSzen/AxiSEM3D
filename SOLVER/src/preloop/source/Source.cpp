// Source.cpp
// base class of source
// sources can be on or off axis 

#include "Source.h"
#include "Quad.h"
#include "Domain.h"
#include "Element.h"
#include "SourceTerm.h"
#include "SourceTerm.h"
#include "Mesh.h"
#include "XMath.h"
#include "SpectralConstants.h"
#include "XMPI.h"
#include "MultilevelTimer.h"
#include "Geodesy.h"

Source::Source(double depth, double lat, double lon):
mDepth(depth), mLatitude(lat), mLongitude(lon) {
    // handle singularity at poles
    if (std::abs(mLatitude - 90.) < tinyDouble) {
        mLatitude = 90.;
        mLongitude = 0.;
    }
    if (std::abs(mLatitude + 90.) < tinyDouble) {
        mLatitude = -90.;
        mLongitude = 0.;
    }
	
	mAxial = true;

}

Source::Source(double depth, double lat, double lon,
	double srcLat, double srcLon, double srcDep):
mDepth(depth), mLatitude(lat), mLongitude(lon) {
    // handle singularity at poles
    if (std::abs(mLatitude - 90.) < tinyDouble) {
        mLatitude = 90.;
        mLongitude = 0.;
    }
    if (std::abs(mLatitude + 90.) < tinyDouble) {
        mLatitude = -90.;
        mLongitude = 0.;
    }
	// compute theta and phi in source-centered coordinate system
	RDCol3 rtpG, rtpS;
	rtpG(0) = 1.;
	rtpG(1) = Geodesy::lat2Theta_d(mLatitude, mDepth);
	rtpG(2) = Geodesy::lon2Phi(mLongitude);
	rtpS = Geodesy::rotateGlob2Src(rtpG, srcLat, srcLon, srcDep);
    mThetaSrc = rtpS(1);
    mPhiSrc = rtpS(2);
	
	mAxial = false;
	
}


void Source::release(Domain &domain, const Mesh &mesh, int isource) const {
    MultilevelTimer::begin("Locate Source", 2);
	
	if (mAxial) { // on axis 
	    // locate local
	    int myrank = XMPI::nproc();
	    int locTag;
	    RDColP interpFactZ, interpFactXii, interpFactEta;
	    if (locate(mesh, locTag, interpFactZ)) {
	        myrank = XMPI::rank();
	    }

	    // min recRank
	    int myrank_min = XMPI::min(myrank);
	    if (myrank_min == XMPI::nproc()) {
	        throw std::runtime_error("Source::release || Error locating source.");
	    }
	    MultilevelTimer::end("Locate Source", 2);

	    MultilevelTimer::begin("Compute Source", 2);
	    // release to me
	    if (myrank_min == XMPI::rank()) {
	        // compute source term
	        arPP_CMatX3 fouriers;
	        const Quad *myQuad = mesh.getQuad(locTag);
	        computeSourceFourier(*myQuad, interpFactZ, interpFactXii, interpFactEta, mPhiSrc, fouriers);
	        // add to domain
	        Element *myElem = domain.getElement(myQuad->getElementTag());
	        domain.addSourceTerm(new SourceTerm(myElem, fouriers, isource));
	    }
	    MultilevelTimer::end("Compute Source", 2);
	} else { // off axis
		// locate local
		int myrank = XMPI::nproc();
		int locTag;
		RDColP interpFactZ, interpFactXii, interpFactEta;
		if (locate(mesh, locTag, interpFactXii, interpFactEta)) {
			myrank = XMPI::rank();
		}

		// min recRank
		int myrank_min = XMPI::min(myrank);
		if (myrank_min == XMPI::nproc()) {
			throw std::runtime_error("OffAxisSource::release || Error locating off-axis source.");
		}
		MultilevelTimer::end("Locate Off-axis Source", 2);

		MultilevelTimer::begin("Compute Off-axis Source", 2);
		// release to me
		if (myrank_min == XMPI::rank()) {
			// compute OffAxisSource term
			arPP_CMatX3 fouriers;
			const Quad *myQuad = mesh.getQuad(locTag);
			computeSourceFourier(*myQuad, interpFactZ, interpFactXii, interpFactEta, mPhiSrc, fouriers);
			// add to domain
			Element *myElem = domain.getElement(myQuad->getElementTag());
	        domain.addSourceTerm(new SourceTerm(myElem, fouriers, isource));
		}
		MultilevelTimer::end("Compute Off-axis Source", 2);
	}
}

bool Source::locate(const Mesh &mesh, int &locTag, RDColP &interpFactZ) const {
    MultilevelTimer::begin("R Source", 3);
    RDCol2 srcCrds = RDCol2::Zero();
    srcCrds(1) = mesh.computeRadiusRef(mDepth, mLatitude, mLongitude);
    MultilevelTimer::end("R Source", 3);

    // check range of subdomain
    if (srcCrds(0) > mesh.sMax() + tinySingle || srcCrds(0) < mesh.sMin() - tinySingle) {
        return false;
    }
    if (srcCrds(1) > mesh.zMax() + tinySingle || srcCrds(1) < mesh.zMin() - tinySingle) {
        return false;
    }
    // find host element
    RDCol2 srcXiEta;
    for (int iloc = 0; iloc < mesh.getNumQuads(); iloc++) {
        const Quad *quad = mesh.getQuad(iloc);
        if (!quad->isAxial() || quad->isFluid() || !quad->nearMe(srcCrds(0), srcCrds(1))) {
            continue;
        }
        if (quad->invMapping(srcCrds, srcXiEta)) {
            if (std::abs(srcXiEta(1)) <= 1.000001) {
                if (std::abs(srcXiEta(0) + 1.) > tinySingle) {
                    throw std::runtime_error("Source::locate || Bad source location.");
                }
                locTag = iloc;
                XMath::interpLagrange(srcXiEta(1), nPntEdge,
                    SpectralConstants::getP_GLL().data(), interpFactZ.data());
                return true;
            }
        }
    }
    return false;
}


bool Source::locate(const Mesh &mesh, int &locTag, 
	RDColP &interpFactXii, RDColP &interpFactEta) const {
    MultilevelTimer::begin("R OffAxisSource", 3);
    RDCol2 srcCrds = RDCol2::Zero();
    double r = mesh.computeRadiusRef(mDepth, mLatitude, mLongitude);
	srcCrds(0) = r * sin(mThetaSrc);
	srcCrds(1) = r * cos(mThetaSrc);
    MultilevelTimer::end("R OffAxisSource", 3);

    // check range of subdomain
    if (srcCrds(0) > mesh.sMax() + tinySingle || srcCrds(0) < mesh.sMin() - tinySingle) {
        return false;
    }
    if (srcCrds(1) > mesh.zMax() + tinySingle || srcCrds(1) < mesh.zMin() - tinySingle) {
        return false;
    }
    // find host element
    RDCol2 srcXiEta;
    for (int iloc = 0; iloc < mesh.getNumQuads(); iloc++) {
        const Quad *quad = mesh.getQuad(iloc);
        if (quad->isFluid() || !quad->nearMe(srcCrds(0), srcCrds(1))) {
            continue;
        }
        if (quad->invMapping(srcCrds, srcXiEta)) {
            if (std::abs(srcXiEta(0)) <= 1.000001 && std::abs(srcXiEta(1)) <= 1.000001) {
                locTag = iloc;
				XMath::interpLagrange(srcXiEta(0), nPntEdge,
                    (quad->isAxial() ? SpectralConstants::getP_GLJ().data()
								     : SpectralConstants::getP_GLL().data()), 
					interpFactXii.data());
                XMath::interpLagrange(srcXiEta(1), nPntEdge,
                    SpectralConstants::getP_GLL().data(), interpFactEta.data());
                return true;
            }
        }
    }
    return false;
}



