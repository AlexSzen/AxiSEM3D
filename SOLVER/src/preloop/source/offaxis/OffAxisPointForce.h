// OffAxisPointForce.h
// created by Kuangdai on 11-Nov-2017
// off-axis point-force source

#pragma once
#include "Source.h"

class OffAxisPointForce: public Source {
public:
    
    OffAxisPointForce(double depth, double lat, double lon,
		double srcLat, double srcLon, double srcDep);

    std::string verbose() const;

protected:
	void computeSourceFourier(const Quad &myQuad, 
		const RDColP &interpFactZ, 
		const RDColP &interpFactXii,
		const RDColP &interpFactEta,
		double phi,
		arPP_CMatX3 &fouriers) const;

private:
    // store: cylindrical components (q_s, q_phi, q_z)

};
