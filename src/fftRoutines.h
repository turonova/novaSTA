#pragma once
#include <fftw3.h>
#include <vector>
#include <complex>

class FFTRoutines
{
public:
    static void real1DTransform(size_t size, std::vector<float>& data, std::vector<float>& fft);
    static void real2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<double>& mask);
    static void real2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter);
    static void computePowerSpectrum(std::vector<float>& powerSpectrum, fftw_complex* fftOut, size_t sizeX, size_t nyh, bool logarithmizeData);
    static void complex2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter);
    static void many1DTransform(float* data, int sizeX, int sizeY, int direction);
    static void real3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& data, std::vector<std::complex<double>>& filter);
    static void complex3DTransform(size_t sizeX, size_t sizeY, size_t  sizeZ, std::vector<float>& data, std::vector<std::complex<double>>& filter);

    static void forwardComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<std::complex<float>>& out);
    static void forwardComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<std::complex<double>>& out);
    static void forwardReal3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<std::complex<float>>& out);
    static void forwardComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<float>& out);

    static void forwardComplex3DTransformPS(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<float>& out);

    static void inverseReal3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<std::complex<float>>& in, std::vector<float>& out);
    static void inverseReal3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<std::complex<double>>& in, std::vector<float>& out);
    static void inverseComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<std::complex<float>>& in, std::vector<float>& out);
    static void inverseComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<std::complex<double>>& in, std::vector<float>& out);

    static void computeConjugateOfSquare(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<std::complex<double>>& in, std::vector<std::complex<double>>& out);
    static void computeConjugateOfSquare(size_t sizeX, size_t sizeY, size_t sizeZ,  std::vector<std::complex<float>>& in, std::vector<std::complex<float>>& out);
    static void fftshift(std::vector<float>& in, size_t xdim, int direction);
    static void fftshift(std::vector<float>& in, size_t xdim, size_t ydim, int direction);
    static void fftshift(std::vector<float>& in, std::vector<float>& out, size_t xdim, size_t ydim, int direction);

    template <typename T>
    static void fftshift(std::vector<T>& in, std::vector<T>& out, size_t xdim, size_t ydim, size_t zdim, int direction, size_t fftdim);

    template <typename T>
    static void fftshift(std::vector<T>& in, std::vector<T>& out, size_t xdim, size_t ydim, size_t zdim, int direction);
	//static void fftshift(std::vector<std::complex<float>>& in, std::vector<std::complex<float>>& out, size_t xdim, size_t ydim, size_t zdim, int direction);

    template <typename T>
    static void fftshift(std::vector<T>& in, size_t xdim, size_t ydim, size_t zdim, int direction, size_t fftdim);

    static void fftshiftEven(std::vector<float>& out, std::vector<float>& in, size_t xdim, size_t ydim);

private:
    static double logarithmizeValue(double value, bool logarithmizeData);
    static void maskFFT(fftw_complex* fftOut, std::vector<double>& mask, size_t sizeX, size_t nyh);
    static void filterFFT(fftw_complex* fftOut, std::vector<float>& filter, size_t sizeX, size_t nyh);
    static void filterFFT(fftw_complex* fftOut, std::vector<std::complex<double>>& filter, size_t sizeX, size_t sizeY, size_t sizeZ, size_t nzh);

   // static void normalizeValues(std::vector<float>& normalizedData, std::vector<double>& originalData, size_t dataSize, novaCTF::DataStats& dataStats);
    static void normalize(float *array, float scale, size_t dataSize);

};
