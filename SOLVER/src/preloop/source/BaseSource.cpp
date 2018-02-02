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
	    RDColP interpFactZ;
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
	        computeSourceFourier(*myQuad, interpFactZ, fouriers);
	        // add to domain
	        Element *myElem = domain.getElement(myQuad->getElementTag());
	        domain.addSourceTerm(new SourceTerm(myElem, fouriers, isource));
	    }
	    MultilevelTimer::end("Compute Source", 2);
	} else { // off axis
		// locate local
		int myrank = XMPI::nproc();
		int locTag;
		RDColP interpFactXii, interpFactEta;
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
			vec_arPP_CMatX3 fouriers;
			const Quad *myQuad = mesh.getQuad(locTag);
			computeSourceFourier(*myQuad, interpFactXii, interpFactEta, mPhiSrc, fouriers);
			// add to domain
			Element *myElem = domain.getElement(myQuad->getElementTag());
		   // TODO: implement this in core
		   //       the existing class SourceTerm is not general enough
		   //       because the source time funciton varies for each source
	         domain.addSourceTerm(new SourceTerm(myElem, fouriers, isource));
		}
		MultilevelTimer::end("Compute Off-axis Source", 2);
	}
}



