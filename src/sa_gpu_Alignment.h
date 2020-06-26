#pragma once

#include "sa_Alignment.h"

/*
#include <map>
#include <mpi.h>
#include <vector>
#include <complex>
#include "emFile.h"
#include "wedgeMask.h"
#include "volume.h"
#include "quaternion.h"
*/

#include <cuda_runtime.h>
#include <cufft.h>
#include "gpuRoutines.h"

class sa_gpu_Alignment : public sa_Alignment
{
public:
    sa_gpu_Alignment(ParameterSetup& iParams, Log& iLog, int iProcessID, int iNumberOfProcesses, MPI_Comm iActiveComm, int iGlobalProcessID);
    ~sa_gpu_Alignment();

private:

    void align(unsigned int iteration);
    void computePartialAverage(unsigned int iteration);

    void prepareSubtomogram(Particle* originalSubtomo, vector<float>& rotatedMask);
    void computeSubtomoRefCC(Particle* subtomo, vector<float>& ccf);

    void  divideByMeanAmplitude(vector<complex<double>>& data);

    void rotateReference(vector<float>& volume, vector<float>& output, Quaternion& rotation);

    void prepareReferenceTexture();

    cudaExtent volumeSize;
    cudaMemcpy3DParms copyParams;

    cudaArray* deviceRef = 0;
    float* rotatedRef;

    size_t volSizeBytes;

	GPURoutines* gpu;

};
