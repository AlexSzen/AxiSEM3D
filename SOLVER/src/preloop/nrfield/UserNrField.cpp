// UserNrField.cpp
// created by Kuangdai on 11-Jun-2017 
// user-defined nr integer field

#include "UserNrField.h"
#include <sstream>
#include "Geodesy.h"

UserNrField::UserNrField(bool useLucky, const std::vector<double> &params): 
NrField(useLucky), mParameters(params) {
    // nothing
}

int UserNrField::getNrAtPoint(const RDCol2 &coords) const {
    // input: coordinates
    // s, z, r, theta, depth
    double s = coords(0);
    double z = coords(1);
    double r, theta;
    Geodesy::rtheta(coords, r, theta);
    double depth = Geodesy::getROuter() - r;
    
    // output: Fourier Expansion order nu
    int nu = 2;
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    // TODO: Compute nu from s, z, r, theta, depth and mParameters
    // NOTE: Only edit within this box!
    //       If this is the first time you are looking at this part
    //       of code, the following example implements a Nu field
    //       that only depends on radius (depth), with a very large Nu
    //       assigned near the core-mantle boundary
   /* 
    double rcmb = 3480e3;
    double r_low_broad = rcmb - 500e3;
    double r_low_narrow = rcmb - 100e3;
    double r_upp_narrow = rcmb + 200e3;
    double r_upp_broad = rcmb + 1000e3;
    int nu_narrow = 1000;
    int nu_broad = 100; 
    
    if (r <= r_low_broad) {
        nu = nu_broad;
    } else if (r <= r_low_narrow) {
        nu = (nu_narrow - nu_broad) / (r_low_narrow - r_low_broad) * (r - r_low_broad) + nu_broad;
    } else if (r <= r_upp_narrow) {
        nu = nu_narrow;
    } else if (r <= r_upp_broad) {
        nu = (nu_narrow - nu_broad) / (r_upp_narrow - r_upp_broad) * (r - r_upp_broad) + nu_broad;
    } else {
        nu = nu_broad;
    }
	*/
	// max nu in box around adjoint source 
	double r_outer = 6.371e6;
	double r_lower_start = r_outer - 100e3; // between r_lower_start and r_outer nu=nu_max
	double r_lower_end = r_outer - 500e3; // below r_lower_end nu = nu_min
	double theta_adj = 0.351219; //source centered theta in rad 
	double d_theta1 = 2.0 * pi / 180.; //2 degrees dtheta for first band 
	double d_theta2 = 5.0 * pi /180.; //5 degrees dtheta for second band 
	
    int nu_max = 200;
	int nu_min = 10;
	
	nu = nu_min;
	
	if ( r >= r_lower_end ) { // max nu in crust 
		
		nu = nu_max;

	}
	
	
    // NOTE: Only edit within this box!
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    // get nr from nu
    int nr = 2 * nu + 1;
    if (nr <= 0) {
        throw std::runtime_error("UserNrField::getNrAtPoint || Non-positive Nr.");
    }
    return nr;
}

std::string UserNrField::verbose() const {
    std::stringstream ss;
    ss << "\n================= Fourier Expansion Order ==================" << std::endl;
    ss << "  Type                     =   User-defined" << std::endl;
    if (mParameters.size() > 0) {
        ss << "  Parameter List           =   ";
        for (int ipar = 0; ipar < mParameters.size(); ipar++) {
            ss << mParameters[ipar] << " ";
        }
        ss << std::endl;
    }
    ss << "  Use FFTW Lucky Numbers   =   " << (mUseLuckyNumber ? "YES" : "NO") << std::endl;
    ss << "================= Fourier Expansion Order ==================\n" << std::endl;
    return ss.str();
}

