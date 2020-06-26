#include "volume.h"
#include "maskRoutines.h"
#include "emFile.h"

using namespace vr;

Volume::Volume(vector<float>& volumeData,vector<size_t>& aDim):Volume(aDim)
{
    for(size_t i=0; i < data.size(); i++)
    {
        data[i]=volumeData[i];
    }
}

Volume::Volume(vector<size_t>& aDim)
{
    interpolationType = 0;
    intepolationGaussCoeff = 1.0f;
    size_t volSize = 1;
    dim.resize(aDim.size());
    center.resize(aDim.size());

    for (size_t i = 0; i < aDim.size(); i++)
    {
        dim[i] = aDim[i];
        volSize *= dim[i];
        center[i] = floor((float)aDim[i] / 2.0f);
    }

    sxy = dim[0] * dim[1];

    fftz = floor(dim[2] / 2.0f) + 1;

    data.resize(volSize, 0.0f);
    rotated.resize(volSize, 0.0f);
    //rotatedMask.resize(volSize, 0.0f);
}

//vector<float>& Volume::rotate(float phi, float psi, float theta)
//{
 //   MaskRoutines::rotateVolume(&data[0], rotated, dim[0], dim[1], dim[2], phi, psi, theta);
 //   return rotated;
//}

vector<float>& Volume::rotateQuat(Quaternion& q)
{
    rotate(q);
    return rotated;
}

void Volume::update(vector<float>& newData)
{
    if(newData.size()!=data.size())
    {
        cout << endl << "ERROR: The volumes do not have the same size!!!" << endl;
        exit(EXIT_FAILURE);
    }

    for(size_t i=0; i< data.size(); i++)
        data[i]=newData[i];
}

//vector<float>& Volume::rotateAndShift(float phi, float psi, float theta, vector<float>& shift)
//{
 //   MaskRoutines::rotateVolume(&data[0], rotated, dim[0], dim[1], dim[2], phi, psi, theta);
  //  MaskRoutines::shiftVolume(rotated, dim[0], dim[1], dim[2], shift[0], shift[1], shift[2]);
   // return rotated;
//}

//vector<float>& Volume::shiftAndRotate(vector<float>& shift, float phi, float psi, float theta)
//{
//    MaskRoutines::shiftVolume(data, dim[0], dim[1], dim[2], shift[0], shift[1], shift[2]);
 //   MaskRoutines::rotateVolume(&data[0], rotated, dim[0], dim[1], dim[2], phi, psi, theta);
//
//
  //  return rotated;
//}

vector<float>& Volume::rotateAndShift(Quaternion& quat, vector<float>& shift)
{
    rotate(quat);
    MaskRoutines::shiftVolume(rotated, dim[0], dim[1], dim[2], shift[0], shift[1], shift[2]);
    return rotated;
}

vector<float>& Volume::shiftAndRotate(vector<float>& shift, Quaternion& quat)
{
    MaskRoutines::shiftVolume(data, dim[0], dim[1], dim[2], shift[0], shift[1], shift[2]);
    rotate(quat);
    return rotated;
}

void Volume::shiftAndRotate(vector<float>& shift, Quaternion& quat, vector<float>& mask, bool binarizeMask)
{
    MaskRoutines::shiftVolume(data, dim[0], dim[1], dim[2], shift[0], shift[1], shift[2]);

    if (binarizeMask)
        rotate(quat, mask, 1);
    else
        rotate(quat, mask, 4);
}

void Particle::complexForwardFFT(Volume::DataToProcess dataType)
{
    if (dataType == Volume::DataToProcess::maskedD)
        FFTRoutines::forwardComplex3DTransform(dim[0], dim[1], dim[2], masked, complexFFT);
    else if (dataType == Volume::DataToProcess::rotatedD)
        FFTRoutines::forwardComplex3DTransform(dim[0], dim[1], dim[2], rotated, complexFFT);
    else if (dataType == Volume::DataToProcess::originalD)
        FFTRoutines::forwardComplex3DTransform(dim[0], dim[1], dim[2], data, complexFFT);
}

void Particle::realForwardFFT(Volume::DataToProcess dataType)
{
    if (dataType == Volume::DataToProcess::maskedD)
        FFTRoutines::forwardReal3DTransform(dim[0], dim[1], dim[2], masked, complexFFT);
    else if (dataType == Volume::DataToProcess::rotatedD)
        FFTRoutines::forwardReal3DTransform(dim[0], dim[1], dim[2], rotated, complexFFT);
    else if (dataType == Volume::DataToProcess::originalD)
        FFTRoutines::forwardReal3DTransform(dim[0], dim[1], dim[2], data, complexFFT);
}

void Particle::maskAndNormalize(vector<float>& aMask, float maskSignal, bool keepOriginal)
{
   // float ds = vec_sum(data);
    masked = data * aMask;

   // float ms = vec_sum(masked);

    // Subtract the mask - weighted mean value
    float normFactor = vec_sum(masked) / maskSignal;

    for (size_t i = 0; i < masked.size(); i++)
    {
        masked[i] = masked[i] - aMask[i] * normFactor;
        
        if (!keepOriginal)
            data[i] = masked[i];
    }

   // double ss = 0;
   // for (size_t i = 0; i < masked.size(); i++)
   //     ss += (double)data[i];
}

void Particle::maskFFT(vector<float>& aMask)
{

    double normFactor = 0.0;

    for (size_t i = 0; i < fftz; i++)
        for (size_t j = 0; j < dim[1]; j++)
            for (size_t k = 0; k < dim[2]; k++)
            {
                size_t index = i + j*fftz + k*fftz*dim[1];
                if (index == 0)
                    continue;

                complexFFT[index] = complex<double>(complexFFT[index].real() * aMask[index], complexFFT[index].imag() * aMask[index]);
                if (i==0 )
                    normFactor = normFactor + (double)(complexFFT[index].real()*complexFFT[index].real()) + (double)(complexFFT[index].imag()*complexFFT[index].imag());
                else
                    normFactor = normFactor + 2.0*(double)(complexFFT[index].real()*complexFFT[index].real()) + 2.0*(double)(complexFFT[index].imag()*complexFFT[index].imag());
             }

    // Set 0 - frequency peak to zero
    complexFFT[0] = complex<float>(0.0, 0.0);


    normFactor = (dim[0]*dim[1]*dim[2]) / sqrt(normFactor);
    vec_mult(complexFFT, normFactor);
}

void Particle::maskFFT2(vector<float>& aMask,bool computeConjugate)
{
    // Set 0 - frequency peak to zero
    complexFFT[0] = complex<float>(0.0, 0.0);

    // mask volume in Fourier space
    for (size_t i = 1; i < complexFFT.size(); i++)
    {
        complexFFT[i] = complex<float>(complexFFT[i].real() * aMask[i], complexFFT[i].imag() * aMask[i]);
    }

   // float elements = dim[0]* dim[1]* dim[2];
    if(computeConjugate)
    {
      //  float volumeSize = 1.0f / (dim[0] * dim[1] * dim[2]);
        FFTRoutines::computeConjugateOfSquare(dim[0],dim[1],dim[2],complexFFT,complexFFTSquaredConj);

        for (size_t i = 0; i < complexFFT.size(); i++)
        {
            // Store complex conjugate
            complexFFT[i] = complex<float>(complexFFT[i].real(), -complexFFT[i].imag());

            // Store complex conjugate of square
            complexFFTSquaredConj[i] = complex<float>(complexFFTSquaredConj[i].real(), -complexFFTSquaredConj[i].imag());
        }
    }
}


vector<float>& Volume::Data(DataToProcess dataType)
{
    if (dataType == Volume::DataToProcess::rotatedD)
        return rotated;
    else if (dataType == Volume::DataToProcess::rotatedMaskD)
        return rotatedMask;
    else if (dataType == Volume::DataToProcess::rotatedCCMaskD)
        return rotatedCCMask;
    else
        return data;
}



Particle::Particle(vector<size_t>& aDim, unsigned int aInterpolationType, float aInterpolationSigma):Volume(aDim)
{
    interpolationType = aInterpolationType;
    intepolationGaussCoeff = -1.0f/(2.0f*aInterpolationSigma*aInterpolationSigma);

    complexFFT.resize(dim[0] * dim[1] * fftz);
    complexFFTSquaredConj.resize(dim[0] * dim[1] * fftz);
}

vector<complex<float>>& Particle::DataFFT()
{
    return complexFFT;
}


vector<complex<float>>& Particle::DataFFTSquaredConj()
{
    return complexFFTSquaredConj;
}


vector<float>& Particle::DataMasked()
{
    return masked;
}

Particle::Particle(vector<float>& aData, vector<size_t>& aDim) :Volume(aDim)
{
    for (size_t i = 0; i < aData.size(); i++)
        data[i] = aData[i];

    complexFFT.resize(dim[0] * dim[1] * fftz);
}

vector<float>& Particle::rotatedParticle()
{
    return rotated;
}

void Particle::getData(string filename)
{
    emFile::read(filename, data);
}

void Particle::divideByMeanAmplitude()
{
    double normFactor = 0.0;
    for (size_t i = 0; i< complexFFT.size(); i++)
    {
        normFactor = normFactor + (complexFFT[i].real()*complexFFT[i].real() + complexFFT[i].imag()*complexFFT[i].imag());
    }
    normFactor = complexFFT.size() / sqrt(normFactor);
    vec_mult(complexFFT, normFactor);
}

Mask::Mask(string filename, vector<size_t> aDim) :Volume(filename,aDim)
{
}

Particle::Particle(string filename, vector<size_t> aDim) :Volume(filename,aDim)
{
}

Volume::Volume(string filename, vector<size_t> aDim):Volume(aDim)
{
    emFile::read(filename, data);
}

void Volume::initWithValue(float value)
{
    rotatedMask.resize(dim[0]*dim[1]*dim[2], value);
    rotatedCCMask.resize(dim[0]*dim[1]*dim[2], value);
}

void Volume::rotate(Quaternion& q)
{
    vector<float> rm;
    q.getRotationMatrix(rm);

    //size_t outputIndex = 0;
    for (size_t k = 0; k < dim[2]; k++)
    {
        for (size_t j = 0; j < dim[1]; j++)
        {
            for (size_t i = 0; i < dim[0]; i++)
            {
                long pi = i - center[0];
                long pj = j - center[1];
                long pk = k - center[2];

                /* transformation of coordinates */
                float rx = center[0] + rm[0] * pi + rm[3] * pj + rm[6] * pk;
                float floorx = floor(rx);

                size_t outputIndex = i+j*dim[0]+k*dim[0]*dim[1];

                if (rx < 0 || floorx >= (dim[0]-1))
                {
                    rotated[outputIndex] = 0.0f;
                   // outputIndex++;
                    continue;
                }
                float ry = center[1] + rm[1] * pi + rm[4] * pj + rm[7] * pk;
                float floory = floor(ry);
                if (ry < 0 || floory >= (dim[1] - 1))
                {
                    rotated[outputIndex] = 0.0f;
                    //outputIndex++;
                    continue;
                }
                float rz = center[2] + rm[2] * pi + rm[5] * pj + rm[8] * pk;
                float floorz = floor(rz);
                if (rz < 0 || floorz >= (dim[2] - 1))
                {
                    rotated[outputIndex] = 0.0f;
                    //outputIndex++;
                    continue;
                }

                /* Interpolation */

                float vx2 = rx - floorx;
                float vx1 = 1 - vx2;

                float vy2 = ry - floory;
                float vy1 = 1 - vy2;

                float vz2 = rz - floorz;
                float vz1 = 1 - vz2;

                /* the following section detects border pixels to avoid exceeding dimensions */
                size_t inputIndex = floorx + floory * dim[0] + floorz * sxy;

                /* interpolation */
                float fb = data[inputIndex] + (data[inputIndex+1] - data[inputIndex]) * vx2;
                float ft = data[inputIndex + dim[0]] * vx1 + data[inputIndex + dim[0] + 1] * vx2;
                float rb = data[inputIndex + sxy] * vx1 + data[inputIndex + sxy + 1] * vx2;
                float rt = data[inputIndex + dim[0] + sxy] * vx1 + data[inputIndex + dim[0] + sxy +1] * vx2;

                rotated[outputIndex] = (fb * vy1 + ft * vy2) * vz1 + (rb * vy1 + rt * vy2) * vz2;
               // outputIndex++;
            }
        }
    }
}

void Mask::binarize()
{
    for (size_t i = 0; i < dim[0] * dim[1] * dim[2]; i++)
    {
        if (data[i]>0.0f)
            data[i] = 1.0f;
        else
            data[i] = 0.0f;
    }
}

void Volume::rotate(Quaternion& q, vector<float>& mask, unsigned int rotationType)
{
    vector<float> dummyMask;
    rotate(q, mask, dummyMask, rotationType);
}

float Volume::interpolation(vector<float>& inputVolume, size_t inputIndex, float vx1, float vx2, float vy1, float vy2, float vz1, float vz2)
{
    float fb = inputVolume[inputIndex] + (inputVolume[inputIndex + 1] - inputVolume[inputIndex]) * vx2;
    float ft = inputVolume[inputIndex + dim[0]] * vx1 + inputVolume[inputIndex + dim[0] + 1] * vx2;
    float rb = inputVolume[inputIndex + sxy] * vx1 + inputVolume[inputIndex + sxy + 1] * vx2;
    float rt = inputVolume[inputIndex + dim[0] + sxy] * vx1 + inputVolume[inputIndex + dim[0] + sxy + 1] * vx2;

    float finalValue = (fb * vy1 + ft * vy2) * vz1 + (rb * vy1 + rt * vy2) * vz2;

    return finalValue;
}

void Volume::rotate(Quaternion& q, vector<float>& mask, vector<float>& ccMask, unsigned int rotationType)
{
    vector<float> rm;
    q.getRotationMatrix(rm);

    //size_t outputIndex = 0;
    for (size_t k = 0; k < dim[2]; k++)
    {
        for (size_t j = 0; j < dim[1]; j++)
        {
            for (size_t i = 0; i < dim[0]; i++)
            {
                long pi = i - center[0];
                long pj = j - center[1];
                long pk = k - center[2];

                /* transformation of coordinates */
                float rx = center[0] + rm[0] * pi + rm[3] * pj + rm[6] * pk;
                float floorx = floor(rx);

                size_t outputIndex = i+j*dim[0]+k*dim[0]*dim[1];

                if (rx < 0 || floorx >= (dim[0]-1))
                {
                    if (rotationType != 1 && rotationType != 4)
                    {
                        rotated[outputIndex] = 0.0f;
                        rotatedMask[outputIndex] = 0.0f;
                        rotatedCCMask[outputIndex] = 0.0f;
                    }

                    continue;
                }
                float ry = center[1] + rm[1] * pi + rm[4] * pj + rm[7] * pk;
                float floory = floor(ry);
                if (ry < 0 || floory >= (dim[1] - 1))
                {
                    if (rotationType != 1 && rotationType != 4)
                    {
                        rotated[outputIndex] = 0.0f;
                        rotatedMask[outputIndex] = 0.0f;
                        rotatedCCMask[outputIndex] = 0.0f;
                    }

                    continue;
                }
                float rz = center[2] + rm[2] * pi + rm[5] * pj + rm[8] * pk;
                float floorz = floor(rz);
                if (rz < 0 || floorz >= (dim[2] - 1))
                {
                    if (rotationType != 1 && rotationType != 4)
                    {
                        rotated[outputIndex] = 0.0f;
                        rotatedMask[outputIndex] = 0.0f;
                        rotatedCCMask[outputIndex] = 0.0f;
                    }

                    continue;
                }

                /* Interpolation */

                float vx2 = rx - floorx;
                float vx1 = 1 - vx2;

                float vy2 = ry - floory;
                float vy1 = 1 - vy2;

                float vz2 = rz - floorz;
                float vz1 = 1 - vz2;

                // gauss
                if (interpolationType==1)
                {
                    vx1 = exp((vx2*vx2) * intepolationGaussCoeff);
                    vx2 = exp((vx1*vx1) * intepolationGaussCoeff);
                    vy1 = exp((vy2*vy2) * intepolationGaussCoeff);
                    vy2 = exp((vy1*vy1) * intepolationGaussCoeff);
                    vz1 = exp((vz2*vz2) * intepolationGaussCoeff);
                    vz2 = exp((vz1*vz1) * intepolationGaussCoeff);
                }

                /* the following section detects border pixels to avoid exceeding dimensions */
                size_t inputIndex = floorx + floory * dim[0] + floorz * sxy;
                           
                if (rotationType == 1)
                {
                    rotated[outputIndex] += interpolation(data, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2);
                    rotatedMask[outputIndex] += binarizeMask(interpolation(mask, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2));
                }
                else if (rotationType == 4)
                {
                    rotated[outputIndex] += interpolation(data, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2);
                    rotatedMask[outputIndex] += interpolation(mask, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2);
                }
                else
                {
                    rotated[outputIndex] = interpolation(data, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2);
                    rotatedCCMask[outputIndex] = interpolation(ccMask, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2);
                }

                if (rotationType == 3)
                {
                    rotatedMask[outputIndex] = binarizeMask(interpolation(mask, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2));
                }
            }
        }
    }
}

float Volume::binarizeMask(float value)
{
    if (value >= 0.5f)
        value = 1.0f;
    else
        value = 0.0f;

    return value;
}
