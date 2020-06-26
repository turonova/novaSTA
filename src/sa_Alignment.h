#pragma once

#include "parameterSetup.h"
#include <map>
#include <mpi.h>
#include <vector>
#include <complex>
#include "sa_Alignment.h"
#include "emFile.h"
#include "wedgeMask.h"
#include "volume.h"
#include "quaternion.h"
#include "fsc.h"
#include "motl.h"
#include "log.h"

using namespace std;

class sa_Alignment
{
public:

    static sa_Alignment* create(ParameterSetup& iParams, Log& iLog, int iProcessID, int iNumberOfProcesses, MPI_Comm iActiveComm);
    void run();

protected:

    virtual void align(unsigned int iteration) = 0;
    virtual void computePartialAverage(unsigned int iteration) = 0;

    static void splitIntoEvenOdd(ParameterSetup& iParams);
    void extractSubtomos();
    void createStartingReference();
    void computeFinalAverage(unsigned int iteration);

    void initVaribales();
    void prepareReference(unsigned int index);
    void prepareMotl();
    void prepareMasks();
    void clean(unsigned int iteration);
    void copySubtomos();

    void setRotations();
    void setBandPassFilter();

    bool runAlignment(unsigned int iteration);

    void computeFSC(vector<float>& evenVolume, vector<float>& oddVolume,unsigned int iteration);
    void combineEvenOdd(unsigned int iteration);

    void writeIntermediateOutputs(unsigned int iteration);

    float findSubpixelPeak(vector<float>& volume, vector<float>& shift);

    static string createName(string basename, unsigned int number, string typeExtension);
    string createName(string basename, unsigned int iteration, unsigned int number, string typeExtension);


    void updateConvergenceStats(vector<float>& oldShift, vector<float>& newShift, Quaternion& oldRotation, Quaternion& newRotation, float oldPhi, float newPhi);
    void computeFinalStats(unsigned int iteration);

    void preprocessMotl();
    void postprocessMotl(unsigned int iteration);

    float computeAddaptiveFilter();
    
    Log logger;
    ParameterSetup params;
    int processID;
    int numberOfProcesses;
    int globalProcessID;

    WedgeMask* wedgeMask;
    Particle* ref;
    FSC* fsc;

    Mask* mask;
    Mask* ccMask;

    size_t subtomosPerMotl;

    vector<size_t> dim;
    size_t dimx;
    size_t dimy;
    size_t dimz;
    size_t volSize;

    vector<float> convergenceStats;

    float volSignal;

    unsigned int realStartinIndex;
    unsigned int rotationAndFilterIndex;

    Motl* motl;

    vector<float> partialAverage;
    vector<float> partialWedgeMask;

    vector<float> volCenter;

    vector<float> bandPassMaskReduced;
    //vector<float> bandPassMask;

    float phiRange;
    float phiStep;

    size_t numberOfSubtomos;


    vector<Quaternion> rotations;
    vector<float> phiAngles;
    
    string subtomoName;

    MPI_Comm activeComm;

    bool binarizeWedgeMask;

    float currentFscValue;
	bool terminateEarly; 
};
