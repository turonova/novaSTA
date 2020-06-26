#include "fftRoutines.h"
#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <algorithm>

using namespace std;

template void FFTRoutines::fftshift<float>(std::vector<float>& in, std::vector<float>& out, size_t xdim, size_t ydim, size_t zdim, int direction);
template void FFTRoutines::fftshift<std::complex<double>>(std::vector<std::complex<double>>& in, std::vector<std::complex<double>>& out, size_t xdim, size_t ydim, size_t zdim, int direction);
template void FFTRoutines::fftshift<std::complex<float>>(std::vector<std::complex<float>>& in, std::vector<std::complex<float>>& out, size_t xdim, size_t ydim, size_t zdim, int direction);

template void FFTRoutines::fftshift<std::complex<float>>(std::vector<std::complex<float>>& in, size_t xdim, size_t ydim, size_t zdim, int direction,size_t fftdim);

template void FFTRoutines::fftshift<float>(std::vector<float>& in, std::vector<float>& out, size_t xdim, size_t ydim, size_t zdim, int direction, size_t fftdim);
template void FFTRoutines::fftshift<std::complex<double>>(std::vector<std::complex<double>>& in, std::vector<std::complex<double>>& out, size_t xdim, size_t ydim, size_t zdim, int direction, size_t fftdim);


template void FFTRoutines::fftshift<float>(std::vector<float>& in, size_t xdim, size_t ydim, size_t zdim, int direction, size_t fftdim);
template void FFTRoutines::fftshift<std::complex<double>>(std::vector<std::complex<double>>& in, size_t xdim, size_t ydim, size_t zdim, int direction, size_t fftdim);



void FFTRoutines::real1DTransform(size_t size, std::vector<float>& data, std::vector<float>& fft)
{

    double* fftIn;
    fftw_complex* fftOut;
    fftw_plan plan_forward;

    fftIn = (double*)fftw_malloc(sizeof (double)* size);

    for (unsigned int i = 0; i<size; i++)
    {
        fftIn[i] = (double)data[i];
    }

    size_t fftSize = (size / 2) + 1;

    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* fftSize);

    plan_forward = fftw_plan_dft_r2c_1d(size, fftIn, fftOut, FFTW_ESTIMATE);

    fftw_execute(plan_forward);

    if (size % 2 == 0) /* N is even */
        fft.push_back(fftOut[size / 2][0]);

    for (int k = (size + 1) / 2 - 1; k >= 1; k--)  // (k < N/2 rounded up)
        fft.push_back(fftOut[k][0]);

    for (int k = 1; k < (int)(size + 1) / 2; k++)  // (k < N/2 rounded up)
        fft.push_back(fftOut[k][0]);
    if (size % 2 == 0) // N is even
        fft.push_back(fftOut[size / 2][0]);  // Nyquist freq.

    fftw_destroy_plan(plan_forward);

    fftw_free(fftOut);
    fftw_free(fftIn);

    return;

}

double FFTRoutines::logarithmizeValue(double value, bool logarithmizeData)
{
    if (logarithmizeData && value >= 1.0)
        return log(value);
    else
        return value;
}

void inversefftshift(std::vector<float>& out, std::vector<float>& in, size_t xdim, size_t ydim)
{
    for (size_t i = 0; i<ydim / 2; i++)
    {
        for (size_t j = 0; j<xdim / 2; j++)
        {
            size_t outIndex1 = (j + xdim / 2) + (i + ydim / 2)*xdim;
            size_t inIndex1 = j + i*xdim;
            size_t outIndex2 = j + i*xdim;
            size_t inIndex2 = (j + xdim / 2) + (i + ydim / 2)*xdim;
            out[inIndex1] = in[outIndex1];	//1->4
            out[inIndex2] = in[outIndex2];	//4->1
        }
    }

    for (size_t i = 0; i<ydim / 2; i++)
    {
        for (size_t j = xdim / 2; j<xdim; j++)
        {
            size_t outIndex1 = (j - xdim / 2) + (i + ydim / 2)*xdim;
            size_t inIndex1 = j + i*xdim;
            size_t outIndex2 = j + i*xdim;
            size_t inIndex2 = (j - xdim / 2) + (i + ydim / 2)*xdim;
            out[inIndex1] = in[outIndex1];	//2->3
            out[inIndex2] = in[outIndex2];	//3->2
        }
    }
}


void FFTRoutines::fftshiftEven(std::vector<float>& out, std::vector<float>& in, size_t xdim, size_t ydim)
{
    size_t mx = (size_t)floor((float)xdim / 2.0f);
    size_t my = (size_t)floor((float)ydim / 2.0f);

    for (size_t i = 0; i< my; i++)
    {
        for (size_t j = 0; j< mx; j++)
        {
            size_t outIndex1 = (j + mx) + (i + my)*xdim;
            size_t inIndex1 = j + i*xdim;
            size_t outIndex2 = j + i*xdim;
            size_t inIndex2 = (j + mx) + (i + my)*xdim;
            out[outIndex1] = in[inIndex1];	//1->4
            out[outIndex2] = in[inIndex2];	//4->1
        }
    }

    for (size_t i = 0; i<my; i++)
    {
        for (size_t j = mx; j<xdim; j++)
        {
            size_t outIndex1 = (j - mx) + (i + my)*xdim;
            size_t inIndex1 = j + i*xdim;
            size_t outIndex2 = j + i*xdim;
            size_t inIndex2 = (j - mx) + (i + my)*xdim;
            out[outIndex1] = in[inIndex1];	//2->3
            out[outIndex2] = in[inIndex2];	//3->2
        }
    }
}

void FFTRoutines::fftshift(std::vector<float>& in, size_t xdim, int direction)
{
    size_t mx;

    // Necessary for odd dimensions, for even ones it is the same
    if (direction == 1)   //forward
    {
        mx = (size_t)ceil((float)xdim / 2.0f);
    }
    else if (direction == -1)   //inverse
    {
        mx = (size_t)floor((float)xdim / 2.0f);
    }
    
    rotate(in.begin(), in.begin() + mx, in.end());
}

void FFTRoutines::fftshift(std::vector<float>& in, size_t xdim, size_t ydim, int direction)
{   
    size_t mx;
    size_t my;

    // Necessary for odd dimensions, for even ones it is the same
    // ceil and floor are swapped in comparison to the other 2D fftshift routine here but it is not a bug
    if (direction == 1)   //forward
    {
        mx = (size_t)ceil((float)xdim / 2.0f);
        my = (size_t)ceil((float)ydim / 2.0f);
    }
    else if (direction == -1)   //inverse
    {
        mx = (size_t)floor((float)xdim / 2.0f);
        my = (size_t)floor((float)ydim / 2.0f);
    }

    vector<float>::iterator sx = in.begin();
    
    for (size_t i = 0; i < ydim; i++)
    {
        sx = in.begin()+i*xdim;
        rotate(sx, sx + mx, sx + xdim);
    }

    vector<float> temp;
    temp.resize(ydim);

    for (size_t i = 0; i < xdim; i++)
    {
        for (size_t j = 0; j < ydim; j++)
            temp[j] = in[i + j*xdim];

        rotate(temp.begin(), temp.begin() + my, temp.end());

        for (size_t j = 0; j < ydim; j++)
            in[i + j*xdim] = temp[j];
    }
}

void FFTRoutines::fftshift(std::vector<float>& in, std::vector<float>& out, size_t xdim, size_t ydim, int direction)
{
    size_t mx;
    size_t my;

    // Necessary for odd dimensions, for even ones it is the same
    if (direction == 1)   //forward
    {
        mx = (size_t)floor((float)xdim / 2.0f);
        my = (size_t)floor((float)ydim / 2.0f);
    }
    else if (direction == -1)   //inverse
    {
        mx = (size_t)ceil((float)xdim / 2.0f);
        my = (size_t)ceil((float)ydim / 2.0f);
    }

    vector<size_t> indx;
    indx.resize(xdim);
    for (size_t i = 0; i < xdim; i++)
        indx[i] = i;

    rotate(indx.begin(), indx.begin() + mx, indx.end());

    vector<size_t> indy;
    indy.resize(ydim);
    for (size_t i = 0; i < ydim; i++)
        indy[i] = i;

    rotate(indy.begin(), indy.begin() + my, indy.end());

    for (size_t i = 0; i < ydim; i++)
    {
        for (size_t j = 0; j < xdim; j++)
        {
            out[indx[j] + indy[i]*xdim] = in[j + i*xdim];
        }
    }
}

template <typename T>
void FFTRoutines::fftshift(std::vector<T>& in, std::vector<T>& out, size_t xdim, size_t ydim, size_t zdim, int direction)
{
    fftshift(in,out,xdim,ydim,zdim,direction,0);
}

/*void FFTRoutines::fftshift(std::vector<std::complex<float>>& in, std::vector<std::complex<float>>& out, size_t xdim, size_t ydim, size_t zdim, int direction)
{
	fftshift(in,out,xdim,ydim,zdim,direction,0);
}*/

template <typename T>
void FFTRoutines::fftshift(std::vector<T>& in, size_t xdim, size_t ydim, size_t zdim, int direction, size_t fftdim)
{
    vector<T> out;
    out.resize(in.size());

    fftshift(in,out, xdim, ydim, zdim, direction, fftdim);

    copy(out.begin(), out.end(), in.begin());
}

template <typename T>
void FFTRoutines::fftshift(std::vector<T>& in, std::vector<T>& out, size_t xdim, size_t ydim, size_t zdim, int direction, size_t fftdim)
{
    size_t mx;
    size_t my;
    size_t mz;

    // Necessary for odd dimensions, for even ones it is the same
    if (direction == 1)   //forward
    {
        mx = (size_t)floor((float)xdim / 2.0f);
        my = (size_t)floor((float)ydim / 2.0f);
        mz = (size_t)floor((float)zdim / 2.0f);
    }
    else if (direction == -1)   //inverse
    {
        mx = (size_t)ceil((float)xdim / 2.0f);
        my = (size_t)ceil((float)ydim / 2.0f);
        mz = (size_t)ceil((float)zdim / 2.0f);
    }

    size_t is = mz;

    if (fftdim == 0)
    {
        fftdim = zdim;
        is = 0;
    }

    vector<size_t> indx;
    indx.resize(xdim);
    for (size_t i = 0; i < xdim; i++)
        indx[i] = i;

    rotate(indx.begin(), indx.begin() + mx, indx.end());

    vector<size_t> indy;
    indy.resize(ydim);
    for (size_t i = 0; i < ydim; i++)
        indy[i] = i;

    rotate(indy.begin(), indy.begin() + my, indy.end());

    vector<size_t> indz;
    indz.resize(zdim);
    for (size_t i = 0; i < zdim; i++)
        indz[i] = i;

    rotate(indz.begin(), indz.begin() + mz, indz.end());

    for (size_t i = is; i < fftdim; i++)
    {
        for (size_t j = 0; j < ydim; j++)
        {
            for (size_t k = 0; k < xdim; k++)
            {
                out[indx[k] + indy[j]*fftdim +indz[i]*fftdim*ydim] = in[k + j*xdim + i*xdim*ydim];
            }
        }
    }
}

void FFTRoutines::computePowerSpectrum(std::vector<float>& powerSpectrum, fftw_complex* fftOut, size_t sizeX, size_t nyh, bool logarithmizeData)
{
    size_t k = 0;
    for (size_t j = 0; j<nyh; j++)
    {
        for (size_t i = 0; i<sizeX; i++)
        {
            powerSpectrum[k] = logarithmizeValue(fftOut[j + i*nyh][0] * fftOut[j + i*nyh][0] + fftOut[j + i*nyh][1] * fftOut[j + i*nyh][1], logarithmizeData);
            k++;
        }
    }

    for (int j = nyh - 2; j>0; j--)
    {
        powerSpectrum[k] = fftOut[j][0] * fftOut[j][0] + fftOut[j][1] * fftOut[j][1];
        k++;
        for (int i = sizeX - 1; i>0; i--)
        {
            powerSpectrum[k] = logarithmizeValue(fftOut[j + i*nyh][0] * fftOut[j + i*nyh][0] + fftOut[j + i*nyh][1] * fftOut[j + i*nyh][1], logarithmizeData);
            k++;
        }
    }
}

void FFTRoutines::maskFFT(fftw_complex* fftOut, std::vector<double>& mask, size_t sizeX, size_t nyh)
{
    size_t k = 0;
    for (size_t j = 0; j<nyh; j++)
    {
        for (size_t i = 0; i<sizeX; i++)
        {
            if (mask[k] != 1.0f)
            {
                fftOut[j + i*nyh][0] = fftOut[j + i*nyh][0] * mask[k];
                fftOut[j + i*nyh][1] = fftOut[j + i*nyh][1] * mask[k];
            }
            k++;
        }
    }
}

void FFTRoutines::filterFFT(fftw_complex* fftOut, std::vector<float>& filter, size_t sizeX, size_t nyh)
{
    size_t k = 0;
    for (size_t j = 0; j<nyh; j++)
    {
        for (size_t i = 0; i<sizeX; i++)
        {
            if (filter[k] != 1.0f)
            {
                fftOut[j + i*nyh][0] = fftOut[j + i*nyh][0] * filter[k];
                fftOut[j + i*nyh][1] = fftOut[j + i*nyh][1] * filter[k];
            }
            k++;
        }
    }
}

void FFTRoutines::filterFFT(fftw_complex* fftOut, std::vector<std::complex<double>>& filter, size_t sizeX, size_t sizeY, size_t sizeZ, size_t nzh)
{
    for (size_t k = 0; k < sizeX; k++)
    {
        for (size_t j = 0; j<sizeY; j++)
        {
            for (size_t i = 0; i< nzh; i++)
            {
                double a = fftOut[i + j*nzh + k*nzh*sizeY][0];
                fftOut[i + j*nzh + k*nzh*sizeY][0] = fftOut[i + j*nzh + k*nzh*sizeY][0] * filter[i + j*sizeZ + k*sizeZ*sizeY].real() - fftOut[i + j*nzh + k*nzh*sizeY][1] * filter[i + j*sizeZ + k*sizeZ*sizeY].imag();
                fftOut[i + j*nzh + k*nzh*sizeY][1] = a * filter[i + j*sizeZ + k*sizeZ*sizeY].imag() + fftOut[i + j*nzh + k*nzh*sizeY][1] * filter[i + j*sizeZ + k*sizeZ*sizeY].real();
            }
        }
    }
}

/*
void FFTRoutines::normalizeValues(std::vector<float>& normalizedData, std::vector<double>& originalData, size_t dataSize, novaCTF::DataStats& dataStats)
{
    for (size_t i = 0; i < dataSize; i++)
    {
        normalizedData[i] = originalData[i] / (dataSize);

        dataStats.mean += normalizedData[i];
        dataStats.max = max(dataStats.max, normalizedData[i]);
        dataStats.min = min(dataStats.min, normalizedData[i]);
    }
}*/

void FFTRoutines::complex2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter)
{
    fftw_complex* fftOut;
    fftw_complex* fftIn;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY);
    fftIn = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY);

    for (size_t i = 0; i < sizeX*sizeY; i++)
    {
        fftIn[i][0] = (double)data[i];
        fftIn[i][1] = 0.0;
    }

    fftw_plan plan_forward = fftw_plan_dft_2d(sizeX, sizeY, fftIn, fftOut, -1, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    std::vector<float> isFilter;
    isFilter.resize(filter.size());

    inversefftshift(isFilter, filter, sizeX, sizeY);
    filterFFT(fftOut, isFilter, sizeX, sizeY);

    fftw_plan plan_backward = fftw_plan_dft_2d(sizeX, sizeY, fftOut, fftIn, 1, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    fftw_destroy_plan(plan_backward);
    fftw_destroy_plan(plan_forward);
    fftw_free(fftOut);

    for (size_t i = 0; i < data.size(); i++)
    {
        data[i] = fftIn[i][0] / data.size();
    }

    fftw_free(fftIn);

}

void FFTRoutines::real2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter)
{
    std::vector<double> fftIn;
    size_t nyh;
    fftw_complex* fftOut;
    fftw_plan plan_backward;
    fftw_plan plan_forward;

    size_t sliceSize = sizeX*sizeY;
    fftIn.resize(sizeX*sizeY);

    for (size_t i = 0; i < sliceSize; i++)
    {
        fftIn[i] = (double)data[i];
    }

    nyh = (sizeY / 2) + 1;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * nyh);

    plan_forward = fftw_plan_dft_r2c_2d(sizeX, sizeY, &fftIn[0], fftOut, FFTW_ESTIMATE);
    fftw_execute(plan_forward);


    std::vector<float> isFilter;
    isFilter.resize(filter.size());


    inversefftshift(isFilter, filter, sizeX, sizeY);
    filterFFT(fftOut, isFilter, sizeX, nyh);

    plan_backward = fftw_plan_dft_c2r_2d(sizeX, sizeY, fftOut, &fftIn[0], FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    fftw_destroy_plan(plan_backward);
    fftw_destroy_plan(plan_forward);
    fftw_free(fftOut);

    for (size_t i = 0; i < data.size(); i++)
    {
        data[i] = fftIn[i] / data.size();
    }
}

void FFTRoutines::real2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<double>& mask)
{
    std::vector<double> fftIn;
    size_t nyh;
    fftw_complex* fftOut;
    fftw_plan plan_backward;
    fftw_plan plan_forward;

    size_t sliceSize = sizeX*sizeY;
    fftIn.resize(sizeX*sizeY);

    for (size_t i = 0; i < sliceSize; i++)
    {
        fftIn[i] = (double)data[i];
    }

    nyh = (sizeY / 2) + 1;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * nyh);

    plan_forward = fftw_plan_dft_r2c_2d(sizeX, sizeY, &fftIn[0], fftOut, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    maskFFT(fftOut, mask, sizeX, nyh);

    plan_backward = fftw_plan_dft_c2r_2d(sizeX, sizeY, fftOut, &fftIn[0], FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    //normalizeValues(data, fftIn, sliceSize, dataStats);

    fftw_destroy_plan(plan_backward);
    fftw_destroy_plan(plan_forward);
    fftw_free(fftOut);

}

void FFTRoutines::many1DTransform(float* data, int sizeX, int sizeY, int direction)
{
    size_t nx = sizeX;
    size_t ny = sizeY;

    int nxpad = 2 * (nx / 2 + 1);
    float invScale = (float)(1.0f / sqrt((double)nx));

    fftwf_plan plan;

    if (direction == 0)
    {
        plan = fftwf_plan_many_dft_r2c(1, &sizeX, ny, data, NULL, 1, nxpad, (fftwf_complex *)data, NULL, 1, nxpad / 2, FFTW_ESTIMATE);
    }
    else if (direction == 1)
    {
        plan = fftwf_plan_many_dft_c2r(1, &sizeX, ny, (fftwf_complex *)data, NULL, 1, nxpad / 2, data, NULL, 1, nxpad, FFTW_ESTIMATE);
    }

    fftwf_execute(plan);

    normalize(data, invScale, (size_t)nxpad * ny);

    fftwf_destroy_plan(plan);
}


/* Normalize the given number of real elements by the scaling factor */
void FFTRoutines::normalize(float *array, float scale, size_t dataSize)
{
    size_t i;
    for (i = 0; i < dataSize; i++)
        array[i] *= scale;
}

void FFTRoutines::real3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& data, std::vector<std::complex<double>>& filter)
{
    std::vector<double> fftIn;
    fftw_complex* fftOut;
    fftw_plan plan_backward;
    fftw_plan plan_forward;

    size_t volumeSize = sizeX*sizeY*sizeZ;
    fftIn.resize(volumeSize);

    for (size_t i = 0; i < volumeSize; i++)
    {
        fftIn[i] = (double)data[i];
    }

    size_t nzh = (sizeZ / 2) + 1;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY * nzh);

    plan_forward = fftw_plan_dft_r2c_3d(sizeX, sizeY, sizeZ, &fftIn[0], fftOut, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    std::vector<std::complex<double>> isFilter;
    isFilter.resize(sizeX * sizeY * sizeZ);

    fftshift(filter,isFilter, sizeX, sizeY,sizeZ, -1);
    
    /*for (size_t i = 0; i < nzh * sizeY * sizeZ; i++)
    {
        double a = fftOut[i][0];

        fftOut[i][0] = fftOut[i][0] * isFilter[i].real() - fftOut[i][1] * isFilter[i].imag();
        fftOut[i][1] = a * isFilter[i].imag() + fftOut[i][1] * isFilter[i].real();
    }*/

    filterFFT(fftOut, isFilter, sizeX, sizeY, sizeZ, nzh);

    plan_backward = fftw_plan_dft_c2r_3d(sizeX, sizeY, sizeZ, fftOut, &fftIn[0], FFTW_ESTIMATE);
    fftw_execute(plan_backward);
    
    fftw_destroy_plan(plan_backward);
    fftw_destroy_plan(plan_forward);
    fftw_free(fftOut);
    
    for (size_t i = 0; i < data.size(); i++)
    {
        data[i] = fftIn[i] / data.size();
    }
}

void FFTRoutines::complex3DTransform(size_t sizeX, size_t sizeY, size_t  sizeZ, std::vector<float>& data, std::vector<std::complex<double>>& filter)
{
    fftw_complex* fftOut;
    fftw_complex* fftIn;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);
    fftIn = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);

    for (size_t i = 0; i < sizeX*sizeY*sizeZ; i++)
    {
        fftIn[i][0] = (double)data[i];
        fftIn[i][1] = 0.0;
    }

    fftw_plan plan_forward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, fftIn, fftOut, -1, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    std::vector<std::complex<double>> isFilter;
    isFilter.resize(filter.size());

    fftshift(filter, isFilter, sizeX, sizeY, sizeZ, -1);

    for (size_t i = 0; i < sizeX * sizeY * sizeZ; i++)
    {
        double a = fftOut[i][0];

        fftOut[i][0] = fftOut[i][0] * isFilter[i].real() - fftOut[i][1] * isFilter[i].imag();
        fftOut[i][1] = a * isFilter[i].imag() + fftOut[i][1] * isFilter[i].real();
    }

    fftw_plan plan_backward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, fftOut, fftIn, 1, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    fftw_destroy_plan(plan_backward);
    fftw_destroy_plan(plan_forward);
    fftw_free(fftOut);

    for (size_t i = 0; i < data.size(); i++)
    {
        data[i] = fftIn[i][0] / data.size();
    }

    fftw_free(fftIn);

}

void FFTRoutines::forwardComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<std::complex<float>>& out)
{
    fftw_complex* fftOut;
    fftw_complex* fftIn;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);
    fftIn = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);

    for (size_t i = 0; i < sizeX*sizeY*sizeZ; i++)
    {
        fftIn[i][0] = (double)in[i];
        fftIn[i][1] = 0.0;
    }

    fftw_plan plan_forward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, fftIn, fftOut, -1, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    for (size_t i = 0; i < in.size(); i++)
    {
        out[i] = complex<float>(fftOut[i][0],fftOut[i][1]);
    }

    fftw_destroy_plan(plan_forward);
    fftw_free(fftOut);
    fftw_free(fftIn);

}

void FFTRoutines::forwardComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<std::complex<double>>& out)
{
    fftw_complex* fftIn;
    fftIn = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);

    for (size_t i = 0; i < sizeX*sizeY*sizeZ; i++)
    {
        fftIn[i][0] = (double)in[i];
        fftIn[i][1] = 0.0;
    }

    fftw_plan plan_forward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, fftIn, reinterpret_cast<fftw_complex*>(&out[0]), -1, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    fftw_destroy_plan(plan_forward);
    fftw_free(fftIn);

}

/*void FFTRoutines::inverseReal3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<complex<float>>& in, std::vector<float>& out)
{
    fftw_complex* fftIn;
    fftIn = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);

    for (size_t i = 0; i < sizeX*sizeY*sizeZ; i++)
    {
        fftIn[i][0] = (double)in[i].real();
        fftIn[i][1] = (double)in[i].imag();
    }

    vector<double> fftOut;
    fftOut.resize(out.size());

    fftw_plan plan_backward = fftw_plan_dft_c2r_3d(sizeX, sizeY, sizeZ, fftIn, &fftOut[0], FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    for (size_t i = 0; i < out.size(); i++)
    {
        out[i] = fftOut[i] / out.size();
    }

    fftw_destroy_plan(plan_backward);
    fftw_free(fftIn);
}*/

void FFTRoutines::inverseComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<complex<float>>& in, std::vector<float>& out)
{
    fftw_complex* fftOut;
    fftw_complex* fftIn;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);
    fftIn = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);

    for (size_t i = 0; i < sizeX*sizeY*sizeZ; i++)
    {
        fftIn[i][0] = (double)in[i].real();
        fftIn[i][1] = (double)in[i].imag();
    }

    fftw_plan plan_backward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, fftIn, fftOut, 1, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    fftw_destroy_plan(plan_backward);

    for (size_t i = 0; i < out.size(); i++)
    {
        out[i] = fftOut[i][0] / out.size();
    }

    fftw_free(fftOut);
    fftw_free(fftIn);

}

void FFTRoutines::inverseComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<complex<double>>& in, std::vector<float>& out)
{
    fftw_complex* fftOut;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);
 
    fftw_plan plan_backward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, reinterpret_cast<fftw_complex*>(&in[0]), fftOut, 1, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    fftw_destroy_plan(plan_backward);

    for (size_t i = 0; i < out.size(); i++)
    {
        out[i] = fftOut[i][0] / out.size();
    }

    fftw_free(fftOut);

}

void FFTRoutines::inverseReal3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<complex<float>>& in, std::vector<float>& out)
{
    std::vector<float> fftOut(sizeX*sizeY*sizeZ);

    size_t fftZ = floor(sizeZ / 2.0f) + 1;
    std::vector<std::complex<float>> fftIn(sizeX * sizeY*fftZ);
    
    fftwf_plan plan_backward = fftwf_plan_dft_c2r_3d(sizeX, sizeY, sizeZ, reinterpret_cast<fftwf_complex*>(&fftIn[0]), &fftOut[0], FFTW_ESTIMATE);
    
    std::copy(in.begin(), in.end(), fftIn.begin());

    fftwf_execute(plan_backward);

    fftwf_destroy_plan(plan_backward);

    for (size_t i = 0; i < out.size(); i++)
    {
        out[i] = (float)(fftOut[i] / out.size());
    }
}

void FFTRoutines::forwardReal3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<std::complex<float>>& out)
{
    vector<float> fftIn(sizeX * sizeY*sizeZ);

    for (size_t i = 0; i < sizeX*sizeY*sizeZ; i++)
    {
        fftIn[i] = in[i];
    }

    fftwf_plan plan_forward = fftwf_plan_dft_r2c_3d(sizeX, sizeY, sizeZ, &fftIn[0], reinterpret_cast<fftwf_complex*>(&out[0]), FFTW_ESTIMATE);

    fftwf_execute(plan_forward);
    fftwf_destroy_plan(plan_forward);
}

void FFTRoutines::forwardComplex3DTransform(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<float>& out)
{
    size_t volumeSize = sizeX*sizeY*sizeZ;

    fftw_complex* fftIn;
    fftIn = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)*volumeSize);

    fftw_complex* fftOut;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* volumeSize);

    for (size_t i = 0; i < volumeSize; i++)
    {
        fftIn[i][0] = (double)in[i];
        fftIn[i][1] = 0.0;
    }

    fftw_plan plan_forward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, fftIn,fftOut, -1, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    
    for (size_t i = 0; i < volumeSize; i++)
    {
        out[i] = fftOut[i][0] / volumeSize;
    }

    
    fftw_destroy_plan(plan_forward);
    fftw_free(fftOut);
    fftw_free(fftIn);
}

void FFTRoutines::forwardComplex3DTransformPS(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<float>& in, std::vector<float>& out)
{
    size_t volumeSize = sizeX*sizeY*sizeZ;

    fftw_complex* fftIn;
    fftIn = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)*volumeSize);

    fftw_complex* fftOut;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* volumeSize);

    for (size_t i = 0; i < volumeSize; i++)
    {
        fftIn[i][0] = (double)in[i];
        fftIn[i][1] = 0.0;
    }

    fftw_plan plan_forward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, fftIn, fftOut, -1, FFTW_ESTIMATE);
    fftw_execute(plan_forward);


    for (size_t i = 0; i < volumeSize; i++)
    {
        out[i] = sqrt(fftOut[i][0] * fftOut[i][0] + fftOut[i][1] * fftOut[i][1]);// / (float)volumeSize;
    }


    fftw_destroy_plan(plan_forward);
    fftw_free(fftOut);
    fftw_free(fftIn);
}

void FFTRoutines::computeConjugateOfSquare(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<std::complex<double>>& in, std::vector<std::complex<double>>& out)
{
    fftw_complex* fftOut;
    fftOut = (fftw_complex*)fftw_malloc(sizeof (fftw_complex)* sizeX * sizeY*sizeZ);

    fftw_plan plan_backward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, reinterpret_cast<fftw_complex*>(&in[0]), fftOut, 1, FFTW_ESTIMATE);
    fftw_execute(plan_backward);

    fftw_destroy_plan(plan_backward);

    float volumeSize = 1.0f/(sizeX*sizeY*sizeZ);
    for (size_t i = 0; i < out.size(); i++)
    {
        fftOut[i][0] = fftOut[i][0] * fftOut[i][0] * volumeSize;
        fftOut[i][1] = 0.0 ;
    }

    fftw_plan plan_forward = fftw_plan_dft_3d(sizeX, sizeY, sizeZ, fftOut, reinterpret_cast<fftw_complex*>(&out[0]), -1, FFTW_ESTIMATE);
    fftw_execute(plan_forward);

    fftw_destroy_plan(plan_forward);

    fftw_free(fftOut);
}

void FFTRoutines::computeConjugateOfSquare(size_t sizeX, size_t sizeY, size_t sizeZ, std::vector<std::complex<float>>& in, std::vector<std::complex<float>>& out)
{
    size_t fftZ = floor(sizeZ / 2.0f) + 1;
    vector<float> fftOut(sizeX * sizeY*sizeZ);
    std::vector<std::complex<float>> fftIn(sizeX * sizeY*fftZ);
    
    fftwf_plan plan_backward = fftwf_plan_dft_c2r_3d(sizeX, sizeY, sizeZ, reinterpret_cast<fftwf_complex*>(&fftIn[0]), &fftOut[0], FFTW_ESTIMATE);
    
    std::copy(in.begin(), in.end(), fftIn.begin());
    
    fftwf_execute(plan_backward);
    fftwf_destroy_plan(plan_backward);

    float volumeSize = (sizeX*sizeY*sizeZ*sizeX*sizeY*sizeZ);
    for (size_t i = 0; i < sizeX*sizeY*sizeZ; i++)
    {
        fftOut[i] = fftOut[i] * fftOut[i] / volumeSize;
        //fftOut[i] = fftOut[i] / volumeSize;
    }

    fftwf_plan plan_forward = fftwf_plan_dft_r2c_3d(sizeX, sizeY, sizeZ, &fftOut[0], reinterpret_cast<fftwf_complex*>(&out[0]),  FFTW_ESTIMATE);
    fftwf_execute(plan_forward);

    fftwf_destroy_plan(plan_forward);
}
