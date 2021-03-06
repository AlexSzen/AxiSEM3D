// Domain.h
// created by Kuangdai on 8-Apr-2016 
// computational domain

#pragma once
#include <vector>
#include "global.h"

#ifdef _MEASURE_TIMELOOP
    #include "MultilevelTimer.h"
#endif

class Point;
class Element;
class SolidFluidPoint;
class SourceTerm;
class SourceTimeFunction;
class PointwiseRecorder;
class SurfaceRecorder;
class DomainRecorder;
struct MessagingInfo;
struct MessagingBuffer;
struct LearnParameters;
class Kerner;

class Domain {
public:
    Domain();
    ~Domain();
    
    ////////////// methods before time loop //////////////
    // domain setup
    int addPoint(Point *point);
    int addElement(Element *elem);
    void addSourceTerm(SourceTerm *source) {mSourceTerms.push_back(source);};
    void addSTF(SourceTimeFunction *stf) {mSTFs.push_back(stf);};
    void setPointwiseRecorder(PointwiseRecorder *recorderPW) {mPointwiseRecorder = recorderPW;};
	void setSurfaceRecorder(SurfaceRecorder *recorderSF) {mSurfaceRecorder = recorderSF;};
	void setDomainRecorder(DomainRecorder *recorderDM) {mDomainRecorder = recorderDM;};
    void setMessaging(MessagingInfo *msgInfo, MessagingBuffer *msgBuffer) 
        {mMsgInfo = msgInfo; mMsgBuffer = msgBuffer;};
    void addSFPoint(SolidFluidPoint *SFPoint) {mSFPoints.push_back(SFPoint);};
    void setLearnParameters(LearnParameters *lpar) {mLearnPar = lpar;};
    void setKerner(Kerner *kerner) {mKerner = kerner;};    
    // get const components
    const SourceTimeFunction &getSTF() const {return *mSTFs[0];};
    int getNumPoints() const {return mPoints.size();};
    int getNumElements() const {return mElements.size();};
    
    // get pointer components
    Point *getPoint(int index) const {return mPoints[index];};
    Element *getElement(int index) const {return mElements[index];};
    // SourceTerm *getSourceTerm(int index) {return mSourceTerms[index];};

	// get recorder for kernels 
	DomainRecorder *getDomainRecorder() const {return mDomainRecorder;};
    
    // test domain 
    void test() const;
    
    // initialize displacement with tiny random numbers  
    void initDisplTinyRandom() const;
    
    ////////////// methods during time loop //////////////
    // element operations
    void computeStiff() const;
    void applySource(int tstep) const;
    
    // point operations
    void assembleStiff(int phase = 0) const; 
    void updateNewmark(Real dt) const;
    void coupleSolidFluid() const;
    
    // point-wise stations
    void initializeRecorders() const;
    void finalizeRecorders() const;
    void record(int tstep, Real t) const;
    void dumpLeft() const;
    
    // stability
    void checkStability(double dt, int tstep, double t) const;
    
    // statistics of element and point types
    std::string verbose() const;
    
    // cost measurement
    std::string reportCost() const;
    
    // wisdom
    void learnWisdom(int tstep) const;
    void dumpWisdom() const;
	
	// kernels 
	void initializeKerner();
	void finalizeKerner();
	void computeKernels( int tstep );
    
private:
    bool pointInPreviousRank(int myPointTag) const;
    
    // points
    std::vector<Point *> mPoints;
    // elements
    std::vector<Element *> mElements;
    // solid-fluid boundary
    std::vector<SolidFluidPoint *> mSFPoints;
    // source 
    std::vector<SourceTerm *> mSourceTerms;
    // source time function
    std::vector<SourceTimeFunction *> mSTFs;
    // point-wise stations
    PointwiseRecorder *mPointwiseRecorder = 0;
	// surface wavefield
    SurfaceRecorder *mSurfaceRecorder = 0;
	// domain wavefield
	DomainRecorder *mDomainRecorder = 0;
	// kernels 
	Kerner *mKerner = 0;
    // massaging 
    MessagingInfo *mMsgInfo = 0;
    MessagingBuffer *mMsgBuffer = 0;
    
    // timers 
    #ifdef _MEASURE_TIMELOOP
        MyBoostTimer *mTimerElemts;
        MyBoostTimer *mTimerPoints;
        MyBoostTimer *mTimerAssemb;
        MyBoostTimer *mTimerAsWait;
        MyBoostTimer *mTimerOthers;
    #endif
    
    // wisdom
    LearnParameters *mLearnPar;
};


