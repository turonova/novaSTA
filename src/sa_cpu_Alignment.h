#pragma once

#include <vector>
#include "sa_Alignment.h"
#include "emFile.h"
#include "wedgeMask.h"
#include "volume.h"
#include "quaternion.h"


using namespace std;

class sa_cpu_Alignment : public sa_Alignment
{
public:
    sa_cpu_Alignment(ParameterSetup& iParams, Log& iLog, int iProcessID, int iNumberOfProcesses, MPI_Comm iActiveComm,int iGlobalProcessID);
    
    ~sa_cpu_Alignment();

private:

    void align(unsigned int iteration);
    void computePartialAverage(unsigned int iteration);

    void prepareSubtomogram(Particle* originalSubtomo, Quaternion& rotation, vector<float>& shift);
    void computeSubtomoRefCC(Particle* subtomo, vector<float>& ccf);
    void computeRosemanCC(Particle* subtomo, vector<float>& ccf);

    void  divideByMeanAmplitude(vector<complex<double>>& data);

    void updateResults(Particle* subtomogram, vector<float> shifts, Quaternion quat);
};
