// DomainRecorder.cpp
// created by Alex, rewritten for new output syntax.

#include "DomainRecorder.h"
#include "DomainIO.h"
#include "DomainInfo.h"

DomainRecorder::DomainRecorder(int totalRecordSteps, int recordInterval, int bufferSize, bool write,
double srcLat, double srcLon, double srcDep): 
mTotalRecordSteps(totalRecordSteps),
mRecordInterval(recordInterval), mBufferSize(bufferSize), mWrite(write),
mSrcLat(srcLat), mSrcLon(srcLon), mSrcDep(srcDep) {
	mBufferLineTime = 0;
	mBufferLineNu = 0;
	mIO = new DomainIO();
}


DomainRecorder::~DomainRecorder() {
	delete mIO;
}

void DomainRecorder::addElement(const Element* elem) {
	mDomainInfo.push_back(DomainInfo(elem));
}

void DomainRecorder::initialize() {
	
	
	mBufferTime = RColX::Zero(mBufferSize);

	int totNu = 0;
	for (int ielem = 0; ielem < mDomainInfo.size(); ielem++) {
		totNu += mDomainInfo[ielem].getMaxNu() + 1;
	}
	vec_ar6_RMatPP initBuf(totNu, zero_ar6_RMatPP);
	mBufferDisp.assign(mBufferSize, initBuf);
	if (mWrite) {
		mIO->initialize(mTotalRecordSteps, mRecordInterval, mBufferSize, totNu, mSrcLat, mSrcLon, mSrcDep);
	}
}

void DomainRecorder::finalize() {
	if (mWrite) {
		mIO->finalize();
	}
}

void DomainRecorder::record(int tstep, Real t) {
	
	if (tstep % mRecordInterval != 0) {
		return;
	}
	
	// time
	mBufferTime(mBufferLineTime) = t;
	
	for (int ielem = 0; ielem < mDomainInfo.size(); ielem++) {
		mDomainInfo[ielem].feedBuffer(mBufferLineTime, mBufferLineNu, mBufferDisp);
		mBufferLineNu += mDomainInfo[ielem].getMaxNu() + 1;
	}
	
	mBufferLineNu = 0;
	mBufferLineTime++;
	
	// dump and clear buffer
	if (mBufferLineTime == mBufferSize) {
		dumpToFile();
	}
	
}

void DomainRecorder::dumpToFile() {
	if (mWrite) {
		mIO->dumpToFile(mBufferDisp, mBufferTime, mBufferLineTime);
	}
	mBufferLineTime = 0;
}