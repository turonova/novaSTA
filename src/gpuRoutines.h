#pragma once

#include <vector>
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include "quaternion.h"

using namespace std;

class GPURoutines
{
public:

    GPURoutines(size_t dimx, size_t dimy, size_t dim_z, bool useRoseman,int processID);
    ~GPURoutines();

    static void fft2DR2C(float* input, float* output_real, float* output_img, size_t x, size_t y);
    static void fft2DC2R(float* input_real, float* input_img, float* output, size_t x, size_t y);
    static void fft3DR2C(float* input, float* output_real, float* output_img, size_t x, size_t y, size_t z);
    static void fft3DC2R(float* input_real, float* input_img, float* output, size_t x, size_t y, size_t z);
    static void fft3DC2C(float* input_real, float* input_img, float* output, size_t x, size_t y, size_t z);

    void prepareReferenceTexture(vector<float>& reference, vector<float>& mask, vector<float>& ccMask, string filterModeName);
    void prepareSubtomoTexture(string filterModeName);
    void copyAndBindWedgeTexture(vector<float>& wedgeMask);
    void allocateDeviceMemory();
   
    void shiftSubtomogram(vector<float>& subtomo, float shiftX, float shiftY, float shiftZ);
    void rotateSubtomogramAndWedge(vector<float>& subtomo, vector<float>& wedge, Quaternion& rotation, bool binarizeMask);
    void rotateTexture(Quaternion& rotation);

    void computeCC(vector<float>& ccVolume, vector<float>& shiftedWedgeMask);

    void getFFTResult(vector<float>& result);
    void getMaskedSubtomogram(vector<float>& outReal, vector<float>& outImag);

    void getShiftFilter(vector<float>& filterReal, vector<float>& filterImag);

    void getRotatedSubtomoAndWedge(vector<float>& subtomo, vector<float>& wedge);
    void getRotatedReferenceAndCCMask(vector<float>& ref, vector<float>& ccMask);

    void maskSubtomogram(vector<float>& subtomogram, vector<float>& mask);

    void getReferenceAndCCMask(vector<float>& reference, vector<float>& ccmask);

   
    double sumArray(vector<float>& data);
    double computeNormFactor(vector<float2>& data);

private:

    double sumArray(float* deviceArray, size_t arraySize,bool squared = false);
    double computeNormFactor(cufftComplex* deviceArray, size_t arraySize);

    void rotateVectorIndices(vector<int>& ids, size_t vectorSize, int direction);
    void getShiftFilter(vector<float2>& filter);
    
    void computeConjugateOfSquare();
	void distributeDevices(int processID);	

    dim3 threadsPerBlock1D; 
    dim3 numBlocks1D;
    dim3 threadsPerBlock;  
    dim3 numBlocks;
    dim3 threadsPerBlockFFT;
    dim3 numBlocksFFT;
    dim3 volumeDim;

    cudaExtent volumeExtent;
    cudaChannelFormatDesc channelDesc;

    cudaMemcpy3DParms copyParamsSubtomo;
    cudaMemcpy3DParms copyParamsWedge;

    cudaArray* deviceMask = 0;
    cudaArray* deviceCCMask = 0;
    cudaArray* deviceRef = 0;
    cudaArray* deviceSubtomo = 0;
    cudaArray* deviceWedge = 0;

    cufftComplex* shiftFilter;
    float* fftSubtomoInNew;
    cufftComplex* fftSubtomoOut;
    cufftComplex* conjugateSquare;

    cufftHandle planR2C;
    cufftHandle planC2R;

    float* deviceSubtomoRot;
    float* deviceWedgeRot;
    float* deviceRefRot;
    float* deviceCCMaskRot;
    float* deviceMaskRot;

    size_t volumeSize;
    size_t fftSize;
    size_t fftz;

    bool subtomoAllocated;

    float* numerator;
    cufftComplex* meanMasked;
    float* meanSubtomo;
    float* intensitySubtomo;

    bool useRosemanCC;
    bool textureAllocated;
    unsigned int rotationType;

  //  float* rrrRef;
  //  float* rrrMask;
  //  float* rrrCCMask;

};
