// AxisSource.cpp
// created by Kuangdai on 8-May-2016
// base class of AxisSource
// we only consider point AxisSources located on the axis

#include "AxisSource.h"
#include "Quad.h"
#include "Domain.h"
#include "Element.h"
#include "AxisSourceTerm.h"
#include "Mesh.h"
#include "XMath.h"
#include "SpectralConstants.h"
#include "XMPI.h"
#include "MultilevelTimer.h"

AxisSource::AxisSource(double depth, double lat, double lon):
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
}

void AxisSource::release(Domain &domain, const Mesh &mesh) const {
    MultilevelTimer::begin("Locate AxisSource", 2);
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
        throw std::runtime_error("AxisSource::release || Error locating AxisSource.");
    }
    MultilevelTimer::end("Locate AxisSource", 2);

    MultilevelTimer::begin("Compute AxisSource", 2);
    // release to me
    if (myrank_min == XMPI::rank()) {
        // compute AxisSource term
        arPP_CMatX3 fouriers;
        const Quad *myQuad = mesh.getQuad(locTag);
        computeAxisSourceFourier(*myQuad, interpFactZ, fouriers);
        // add to domain
        Element *myElem = domain.getElement(myQuad->getElementTag());
        domain.addAxisSourceTerm(new AxisSourceTerm(myElem, fouriers));
    }
    MultilevelTimer::end("Compute AxisSource", 2);
}

bool AxisSource::locate(const Mesh &mesh, int &locTag, RDColP &interpFactZ) const {
    MultilevelTimer::begin("R AxisSource", 3);
    RDCol2 srcCrds = RDCol2::Zero();
    srcCrds(1) = mesh.computeRadiusRef(mDepth, mLatitude, mLongitude);
    MultilevelTimer::end("R AxisSource", 3);

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
                    throw std::runtime_error("AxisSource::locate || Bad AxisSource location.");
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

    
    if (verbose) {
        XMPI::cout << src->verbose();
    }
}
