#pragma once

#include "parameterSetup.h"

using namespace std;

class FSC
{
public:

    FSC(ParameterSetup& iParams);
    ~FSC(){};

    void computeFSC(vector<float>& vol1, vector<float>& vol2);
    float getFSCValue(float threshold);
    void printFSCCurve(string outputName);

    static float getFSCValue(vector<float>& fsc, float threshold,unsigned int iBoxSize, float iPixelSize);
    static vector<float> computeFSC(vector<float>& vol1, vector<float>& vol2, vector<float>& fscMask, unsigned int iBoxSize, unsigned int iSymmetry, bool iRandomized);

private:

    vector<float> fftRefVolume;
    vector<float> mask;

    unsigned int symmetry;
    unsigned int boxSize;
    bool randomized;

    vector<float> fscCurve;

    size_t fftz;
    size_t fftSize;

    unsigned int numberOfShells;
    float centerOffset;
    float pixelSize;
};
