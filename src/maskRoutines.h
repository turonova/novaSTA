#pragma once
#include <vector>
#include <algorithm>
#include "defines.h"
#include "vectorRoutines.h"
#include "fftRoutines.h"
#include "volume.h"
#include <complex>

using namespace std;

namespace MaskRoutines{

    void createWedgeMask(vector<float>& wedgeMask, int dim, float minAngle, float maxAngle);
    void createWedgeMask(vector<float>& wedgeMask, vector<float>& shifted, vector<float>& shiftedReduced, int dim, float minAngle, float maxAngle, bool maskEdges, bool shift);

    vector<float> createSphericalMask(size_t x, size_t y, size_t z, float radius, float sigma);

    vector<float> createSphericalMask(size_t x, size_t y, size_t z);

    vector<complex<double>> shiftVolume(vector<float>& volume, float x, float y, float z, float sx, float sy, float sz);

    void symmetrizeVolume(Volume& volume, size_t x, size_t y, size_t z, int symmetry);
    void symmetrizeVolume(vector<float>& volume, size_t x, size_t y, size_t z, int symmetry);

    void rotatePoint(float& x, float& y, float& z, float phi, float psi, float theta);
    vector<float> rotatePoint(vector<float>& point, vector<float>& angles);

}
