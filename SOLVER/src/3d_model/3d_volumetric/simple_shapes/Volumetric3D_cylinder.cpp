// Volumetric3D_cylinder.cpp
// created by Kuangdai on 19-Oct-2016 
// a cylinder-shaped heterogeneity

#include "Volumetric3D_cylinder.h"
#include "XMath.h"
#include <sstream>

void Volumetric3D_cylinder::initialize(const std::vector<double> &params) {
    // need at least 6 parameters to make a cylinder
    if (params.size() < 10) throw std::runtime_error("Volumetric3D_cylinder::initialize || "
        "Not enough parameters to initialize a Volumetric3D_cylinder object.");
    
    // initialize location    
    mR1 = params[0] * 1e3;
    mTheta1 = params[1] * degree;
    mPhi1 = params[2] * degree;
    
    mR2 = params[3] * 1e3;
    mTheta2 = params[4] * degree;
    mPhi2 = params[5] * degree;
    
    // initialize Gaussian parameters
    mMaxAxis = params[6];
    mHWHM_lateral = params[7] * 1e3;
    mHWHM_top_bot = params[8] * 1e3;
    
    // initialize reference type
    if (params[9] < 0.5) 
        mReferenceType = ReferenceTypes::Absolute;
    else if (params[9] < 1.5) 
        mReferenceType = ReferenceTypes::Reference1D;
    else if (params[9] < 2.5)
        mReferenceType = ReferenceTypes::Reference3D;
    else 
        mReferenceType = ReferenceTypes::ReferenceDiff;    
        
    try {
        int ipar = 10;
        mChangeVp = (params.at(ipar++) > 0.);
        mChangeVs = (params.at(ipar++) > 0.);
        mChangeRho = (params.at(ipar++) > 0.);
    } catch (std::out_of_range) {
        // nothing
    }    
}

bool Volumetric3D_cylinder::get3dProperties(double r, double theta, double phi, double rElemCenter,
    double &dvpv, double &dvph, double &dvsv, double &dvsh, double &drho) const {
    // zero results
    dvpv = dvph = dvsv = dvsh = drho = 0.;
    
    // distance from point to axis
    RDCol3 rtpPoint1, rtpPoint2, rtpTarget;
    rtpPoint1(0) = mR1;
    rtpPoint1(1) = mTheta1;
    rtpPoint1(2) = mPhi1;
    rtpPoint2(0) = mR2;
    rtpPoint2(1) = mTheta2;
    rtpPoint2(2) = mPhi2;
    rtpTarget(0) = r;
    rtpTarget(1) = theta;
    rtpTarget(2) = phi;
    const RDCol3 &xyzPoint1 = XMath::toCartesian(rtpPoint1);
    const RDCol3 &xyzPoint2 = XMath::toCartesian(rtpPoint2);
    const RDCol3 &xyzTarget = XMath::toCartesian(rtpTarget);
    const RDCol3 &xyzDiff1 = xyzTarget - xyzPoint1;
    const RDCol3 &xyzDiff2 = xyzTarget - xyzPoint2; 
    double length = (xyzPoint1 - xyzPoint2).norm();
    double distToLine = xyzDiff1.cross(xyzDiff2).norm() / length;
    
    // outside range
    if (distToLine > 4. * mHWHM_lateral) return false;
    
    // compute Gaussian lateral
    double stddev = mHWHM_lateral / sqrt(2. * log(2.));
    double gaussian_lateral = mMaxAxis * exp(-distToLine * distToLine / (stddev * stddev * 2.));
    
    // distance from point to ends
    double d1 = xyzDiff1.norm();
    double d2 = xyzDiff2.norm();
    double dmax = std::max(d1, d2);
    double distToTopBot = sqrt(dmax * dmax - distToLine * distToLine) - length;
    
    // outside range
    if (distToTopBot > 4. * mHWHM_top_bot) return false;
    
    double gaussian_topbot = gaussian_lateral;
    if (distToTopBot > 0.) {
        // compute Gaussian top-bottom
        stddev = mHWHM_top_bot / sqrt(2. * log(2.));
        gaussian_topbot = gaussian_lateral * exp(-gaussian_topbot * gaussian_topbot / (stddev * stddev * 2.)); 
    }
    
    // set perturbations    
    if (mChangeVp) dvpv = dvph = gaussian_topbot;
    if (mChangeVs) dvsv = dvsh = gaussian_topbot;
    if (mChangeRho) drho = gaussian_topbot;
    
    return true;
}

std::string Volumetric3D_cylinder::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name            =   cylinder" << std::endl;
    ss << "  Radius_1 / km         =   " << mR1 / 1e3 << std::endl;
    ss << "  Theta_1 / degree      =   " << mTheta1 / degree << std::endl;
    ss << "  Phi_1 / degree        =   " << mPhi1 / degree << std::endl;
    ss << "  Radius_2 / km         =   " << mR2 / 1e3 << std::endl;
    ss << "  Theta_2 / degree      =   " << mTheta2 / degree << std::endl;
    ss << "  Phi_2 / degree        =   " << mPhi2 / degree << std::endl;
    ss << "  Maximum on Axis       =   " << mMaxAxis << std::endl;
    ss << "  HWHM lateral / km     =   " << mHWHM_lateral / 1e3 << std::endl;
    ss << "  HWHM top & bot / km   =   " << mHWHM_top_bot / 1e3 << std::endl;
    ss << "  Reference Type        =   " << ReferenceTypesString[mReferenceType] << std::endl;
    ss << "  Affect VP             =   " << (mChangeVp ? "YES" : "NO") << std::endl;
    ss << "  Affect VS             =   " << (mChangeVs ? "YES" : "NO") << std::endl;
    ss << "  Affect Density        =   " << (mChangeRho ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}
