// Element.cpp
// created by Kuangdai on 27-Mar-2016 
// base class of AxiSEM3D spectral elements

#include "Element.h"
#include "Point.h"

#include "Gradient.h"
#include "PRT.h"

Element::Element(Gradient *grad, PRT *prt, const std::array<Point *, nPntElem> &points):
mGradient(grad), mPRT(prt), mHasPRT(prt != 0) {
    mMaxNr = mMaxNu = -1;
    for (int i = 0; i < nPntElem; i++) {
        mPoints[i] = points[i];
        mMaxNr = std::max(mMaxNr, mPoints[i]->getNr());
        mMaxNu = std::max(mMaxNu, mPoints[i]->getNu());
    }
    if (mHasPRT) {
        mPRT->checkCompatibility(mMaxNr);
    }
}

Element::~Element() {
    delete mGradient;
    if (mHasPRT) {
        delete mPRT;
    }
}

void Element::addSourceTerm(const arPP_CMatX3 &source) const {
    for (int i = 0; i < nPntElem; i++) {
        mPoints[i]->addToStiff(source[i]);
    }
}

bool Element::axial() const {
    return mPoints[0]->axial();
}

std::string Element::costSignature() const {
    std::stringstream ss;
    ss << verbose() << "$DimAzimuth=" << mMaxNr << "$Axial=" << (axial() ? "T" : "F");
    return ss.str();
}

#include "Geodesy.h"
RDMatPP Element::formThetaMat() const {
    RDMatPP theta;
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &sz = mPoints[ipnt++]->getCoords();
            theta(ipol, jpol) = Geodesy::theta(sz);
        }
    }
    return theta;
}

RDMatXX Element::getCoordsOnSide(int side) const {
	int ipol0 = 0, ipol1 = 0, jpol0 = 0, jpol1 = 0;
	if (side == 0) {
		ipol0 = 0;
		ipol1 = nPol;
		jpol0 = jpol1 = 0;
	} else if (side == 1) {
		ipol0 = ipol1 = nPol;
		jpol0 = 0;
		jpol1 = nPol;
	} else if (side == 2) {
		ipol0 = 0;
		ipol1 = nPol;
		jpol0 = jpol1 = nPol;
	} else {
		ipol0 = ipol1 = 0;
		jpol0 = 0;
		jpol1 = nPol;
	}
	RDMatXX sz = RDMatXX::Zero(2, nPntEdge);
	int ipntedge = 0;
	for (int ipol = ipol0; ipol <= ipol1; ipol++) {
		for (int jpol = jpol0; jpol <= jpol1; jpol++) {
			int ipnt = ipol * nPntEdge + jpol;
			sz.col(ipntedge) = mPoints[ipnt]->getCoords();
			ipntedge++;
		}
	}
	return sz;
}

vec_RMatPP Element::getCoordsPoints() const {
	
	vec_RMatPP elemCoords(2,RMatPP::Zero());
	
	for (int ipol = 0; ipol <= nPol; ipol++) {
		for (int jpol = 0; jpol <= nPol; jpol ++) {
			int ipnt = nPntEdge * ipol + jpol;
			RDCol2 coords = mPoints[ipnt]->getCoords();
			elemCoords[0](ipol,jpol) = coords(0);
			elemCoords[1](ipol,jpol) = coords(1);
		}
	}
	return elemCoords;
	
	
}

#include "XMath.h"
bool Element::needDumping(double rmin, double rmax,double tmin,double tmax) {
	RDCol2 sz= mPoints[12]->getCoords(); //midpoint
	RDCol2 rt = Geodesy::rtheta(sz);
	if (rmin<rt(0) && rmax>rt(0) && tmin<rt(1) && tmax>rt(1)) return true;
	return false; 
		
}

