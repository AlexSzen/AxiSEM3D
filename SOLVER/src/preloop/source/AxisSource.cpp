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

#include "Parameters.h"
#include "Earthquake.h"
#include "PointForce.h"
#include "NullAxisSource.h"
#include <fstream>
#include <boost/algorithm/string.hpp>

void AxisSource::buildInparam(AxisSource *&src, const Parameters &par, int verbose) {
    if (src) {
        delete src;
    }

    // null AxisSource
    if (par.getValue<bool>("DEVELOP_NON_Source_MODE")) {
        src = new NullAxisSource();
        if (verbose) {
            XMPI::cout << src->verbose();
        }
        return;
    }

    std::string src_type = par.getValue<std::string>("Source_TYPE");
    std::string src_file = par.getValue<std::string>("Source_FILE");

    if (boost::iequals(src_type, "earthquake")) {
        std::string cmtfile = Parameters::sInputDirectory + "/" + src_file;
        double depth, lat, lon;
        double Mrr, Mtt, Mpp, Mrt, Mrp, Mtp;
        if (XMPI::root()) {
            std::fstream fs(cmtfile, std::fstream::in);
            if (!fs) {
                throw std::runtime_error("AxisSource::buildInparam || "
                    "Error opening CMT data file: ||" + cmtfile);
            }
            std::string junk;
            std::getline(fs, junk);
            std::getline(fs, junk);
            std::getline(fs, junk);
            std::getline(fs, junk);
            fs >> junk >> lat;
            fs >> junk >> lon;
            fs >> junk >> depth;
            fs >> junk >> Mrr;
            fs >> junk >> Mtt;
            fs >> junk >> Mpp;
            fs >> junk >> Mrt;
            fs >> junk >> Mrp;
            fs >> junk >> Mtp;
            depth *= 1e3;
            Mrr *= 1e-7;
            Mtt *= 1e-7;
            Mpp *= 1e-7;
            Mrt *= 1e-7;
            Mrp *= 1e-7;
            Mtp *= 1e-7;
            fs.close();
        }
        XMPI::bcast(depth);
        XMPI::bcast(lat);
        XMPI::bcast(lon);
        XMPI::bcast(Mrr);
        XMPI::bcast(Mtt);
        XMPI::bcast(Mpp);
        XMPI::bcast(Mrt);
        XMPI::bcast(Mrp);
        XMPI::bcast(Mtp);
        src = new Earthquake(depth, lat, lon, Mrr, Mtt, Mpp, Mrt, Mrp, Mtp);
    } else if (boost::iequals(src_type, "point_force")) {
        // point force
        std::string pointffile = Parameters::sInputDirectory + "/" + src_file;
        double depth, lat, lon;
        double f1, f2, f3;
        if (XMPI::root()) {
            std::fstream fs(pointffile, std::fstream::in);
            if (!fs) {
                throw std::runtime_error("AxisSource::buildInparam || "
                    "Error opening point force data file: ||" + pointffile);
            }
            std::string junk;
            std::getline(fs, junk);
            std::getline(fs, junk);
            std::getline(fs, junk);
            std::getline(fs, junk);
            fs >> junk >> lat;
            fs >> junk >> lon;
            fs >> junk >> depth;
            fs >> junk >> f1;
            fs >> junk >> f2;
            fs >> junk >> f3;
            depth *= 1e3;
            fs.close();
        }
        XMPI::bcast(depth);
        XMPI::bcast(lat);
        XMPI::bcast(lon);
        XMPI::bcast(f1);
        XMPI::bcast(f2);
        XMPI::bcast(f3);
        src = new PointForce(depth, lat, lon, f1, f2, f3);
    } else {
         throw std::runtime_error("AxisSource::buildInparam || Unknown AxisSource type: " + src_type);
    }
    
    if (verbose) {
        XMPI::cout << src->verbose();
    }
}
