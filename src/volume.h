#pragma once

#include <vector>
#include <complex>
#include "quaternion.h"

using namespace std;

class Volume
{
public:

    enum DataToProcess{
        originalD,
        rotatedD,
        maskedD,
        shiftedD,
        rotatedMaskD,
        rotatedCCMaskD
    };

    Volume(){};
    Volume(vector<float>& volumeData,vector<size_t>& aDim);
    Volume(vector<size_t>& aDim);
    Volume(string filename, vector<size_t> aDim);
    ~Volume(){};

    vector<float>& rotateQuat(Quaternion& q);

    vector<float>& shiftAndRotate(vector<float>& shift, Quaternion& quat);
    vector<float>& rotateAndShift(Quaternion& quat, vector<float>& shift);
    vector<float>& fftShift(vector<float>& shift);
    vector<float>& Data(DataToProcess dataType = DataToProcess::originalD);
    void shiftAndRotate(vector<float>& shift, Quaternion& quat, vector<float>& mask, bool binarizeMask);
    
    void rotate(Quaternion& q);

    void rotate(Quaternion& q, vector<float>& mask, unsigned int rotationType);

    void rotate(Quaternion& q, vector<float>& mask, vector<float>& ccMask, unsigned int rotationType);

    void update(vector<float>& newData);

    void initWithValue(float value);
	void initWithCCMask(vector<float>& ccMask);

protected:
    
    float interpolation(vector<float>& inputVolume, size_t inputIndex, float vx1, float vx2, float vy1, float vy2, float vz1, float vz2);
    float binarizeMask(float value);

    vector<float> data;
    vector<float> rotated;
    vector<float> shifted;
    vector<float> fftShifted;
    vector<float> rotatedMask;
    vector<float> rotatedCCMask;

    vector<size_t> dim;
    vector<float> center;
    float sxy;
    size_t fftz;

    unsigned int interpolationType;
    float intepolationGaussCoeff;
};


class Particle : public Volume
{
public:
    
    Particle(string filename, vector<size_t> aDim);
    Particle(vector<size_t>& aDim, unsigned int aInterpolationType = 0, float aInterpolationSigma=1.0f);
    Particle(vector<float>& aData, vector<size_t>& aDim);
    ~Particle(){};
    
    vector<float>& rotatedParticle();
    void maskAndNormalize(vector<float>& aMask, float maskSignal, bool keepOriginal = true);
    void maskFFT(vector<float>& aMask);
    void maskFFT2(vector<float>& aMask, bool computeConjugate = true);
    void complexForwardFFT(Volume::DataToProcess dataType);
    void realForwardFFT(Volume::DataToProcess dataType);
    void getData(string filename);
    vector<complex<float>>& DataFFT();
    vector<complex<float>>& DataFFTSquaredConj();
    vector<float>& DataMasked();

private:

    void divideByMeanAmplitude();

    vector<complex<float>> complexFFT;
    vector<complex<float>> complexFFTSquaredConj;
    vector<float> masked;
};

class Mask : public Volume
{
public:

    Mask(string filename, vector<size_t> aDim);
    ~Mask(){};

    vector<float>& rotatedMask();
    void binarize();
};
