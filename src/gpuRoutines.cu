#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <algorithm>
#include "gpuRoutines.h"
#include "math_constants.h"
#include "quaternion.h"
#include "emFile.h"

texture<float, 3, cudaReadModeElementType> texRef;      // 3D texture for reference
texture<float, 3, cudaReadModeElementType> texSubtomo;  // 3D texture for subtomo
texture<float, 3, cudaReadModeElementType> texWedge;    // 3D texture for wedge
texture<float, 3, cudaReadModeElementType> texCCMask;    // 3D texture for ccmask
texture<float, 3, cudaReadModeElementType> texMask;    // 3D texture for mask

using namespace std;

/********************/
/* CUDA ERROR CHECK */
/********************/
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) { getchar(); exit(code); }
    }
}

__global__ void gpu_sumArray(float *input, double *output, size_t elements)
{
    double localSum = 0;

    for (size_t i = threadIdx.x; i< elements; i += gridDim.x*blockDim.x)
    {
        if ((blockIdx.x*blockDim.x + i) >= elements)
            break;

        localSum += input[blockIdx.x*blockDim.x + i];
    }

    output[blockIdx.x*blockDim.x + threadIdx.x] = localSum;
}

__global__ void gpu_sumArraySquared(float* input, double *output, size_t elements)
{
    double localSum = 0;

    for (size_t i = threadIdx.x; i< elements; i += gridDim.x*blockDim.x)
    {
        if ((blockIdx.x*blockDim.x + i) >= elements)
            break;

        localSum = localSum + input[blockIdx.x*blockDim.x + i] * input[blockIdx.x*blockDim.x + i];
    }

    output[blockIdx.x*blockDim.x + threadIdx.x] = localSum;
}


__global__ void gpu_prepareNormValue(cufftComplex* input, float* output, dim3 dim)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;
    size_t k = blockIdx.z * blockDim.z + threadIdx.z;
    
    if( i >= dim.x || j >= dim.y || k>= dim.z)
        return;
         
    size_t index = i + j*dim.x + k*dim.x*dim.y;
    float value = input[index].x*input[index].x + input[index].y*input[index].y;   
    
    if ( i==0 )
        output[index] = value;
    else
         output[index] = 2.0*value;
}

__global__ void gpu_computeNormFactor(cufftComplex *input, double *output, size_t elements)
{
    double localSum = 0;

    for (size_t i = threadIdx.x; i< elements; i += gridDim.x*blockDim.x)
    {
        if ((blockIdx.x*blockDim.x + i) >= elements)
            break;

        localSum = localSum + (input[blockIdx.x*blockDim.x + i].x*input[blockIdx.x*blockDim.x + i].x) + (input[blockIdx.x*blockDim.x + i].y*input[blockIdx.x*blockDim.x + i].y);
    }

    output[blockIdx.x*blockDim.x + threadIdx.x] = localSum;
}

__device__  cuComplex my_expf(float scale, cuComplex value)
{
    cuComplex res;

    float t = expf(scale*value.x);
    sincosf(scale*value.y, &res.y, &res.x);

    res.x *= t;
    res.y *= t;
    
    return res;

}

__global__ void gpu_computeFLCC(float* ccf, float* numerator, float* intensitySubtomo, float* meanSubtomo, float maskSum, float sigmaRef, dim3 dim)
{
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x >= dim.x || y >= dim.y || z >= dim.z)
        return;

    size_t shift = floor(dim.x / 2.0f) + 1;
    size_t volSize = dim.x*dim.y*dim.z;

    size_t inputIndex = volSize - 1 - (x + y*dim.x + z*dim.x*dim.y);
    size_t outputIndex = (x + shift) % dim.x + ((y + shift) % dim.y)*dim.x + ((z + shift) % dim.z)*dim.x*dim.y;
    double denominator = sqrt(intensitySubtomo[inputIndex] - meanSubtomo[inputIndex] * meanSubtomo[inputIndex] / maskSum)*sigmaRef;
    if (denominator != 0.0)
        ccf[outputIndex] = numerator[inputIndex] / denominator;
    else
        ccf[outputIndex] = -1.0f;
}

__global__ void gpu_generateShiftFilter(cufftComplex* filter, int x, int y, int z, float3 shift)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;

    float fftz = floor(z / 2.0f) + 1;

    float ox = floor(x / 2.0f);
    float oy = floor(y / 2.0f);
    float oz = floor(z / 2.0f);

    if ((i + (int)ox) % x >= fftz || j >= y || k >=z)
        return;

    shift.x /= (float)x;
    shift.y /= (float)y;
    shift.z /= (float)z;

    cuComplex c = { 0, 1 };

    double expValue = (i - ox)*shift.x + (j - oy)*shift.y + (k - oz)*shift.z;
  //  int voxelIndex = (k + (int)oz) % z + ((i + (int)ox) % x) * fftz + ((j + (int)oy) % y) * fftz*x;
    size_t voxelIndex = (i + (int)ox) % x + (j + (int)oy) % y * fftz + (k + (int)oz) % z* fftz*x;
    filter[voxelIndex] = my_expf(-2.0* CUDART_PI_F *expValue, c);
}

__global__ void createWedgeMask(float* wedgeMask, int dim, float minAngle, float maxAngle)
{
    double minAngleRad = minAngle*CUDART_PI_F / 180.0;
    double maxAngleRad = maxAngle*CUDART_PI_F / 180.0;

    double tan_min = tan((-CUDART_PI_F / 2.0) - minAngleRad);
    double tan_max = tan(CUDART_PI_F / 2.0 - maxAngleRad);

    //wedgeMask.resize((size_t)dim*dim*dim);
//    fill(wedgeMask.begin(), wedgeMask.end(), 1.0f);

    float halfDim = floor(dim / 2.0);
    for (float z = -halfDim; z < -halfDim + dim; z++)
    {
        if (z == 0)
            continue;

        for (float x = -halfDim; x < -halfDim + dim; x++)
        {
            if ((tan_max > x / z) && (tan_min < x / z))
            {
                for (int y = 0; y < dim; y++)
                {
                    wedgeMask[(size_t)((x + halfDim) + y*dim + (z + halfDim)*dim*dim)] = 0.0f;
                }
            }
        }
    }
}

// Complex pointwise multiplication
__global__ void complexMultiplication(cufftComplex* a, cufftComplex* b, size_t elements)
{

    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= elements)
        return;

   // size_t i = ix + j*dim.x + k*dim.x*dim.y;
  
    float orValue = a[i].x;
    a[i].x = orValue * b[i].x - a[i].y * b[i].y;
    a[i].y = orValue * b[i].y + a[i].y * b[i].x;

}

// Complex pointwise multiplication
__global__ void gpu_complexMultiplication(cufftComplex* a, cufftComplex* b, cufftComplex* c, size_t elements)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= elements)
        return;

    c[i].x = a[i].x * b[i].x - a[i].y * b[i].y;
    c[i].y = a[i].x * b[i].y + a[i].y * b[i].x;
}

// Complex pointwise conjugate multiplication with normalization
__global__ void gpu_conjugateMultiplication(cufftComplex *a, cufftComplex *b, float normFactor, size_t elements)
{
    //blockDim.x * gridDim.x;
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= elements)
        return;

    cufftComplex cn;
    cn.x =  a[i].x * normFactor;
    cn.y = -a[i].y * normFactor;

    a[i].x = cn.x * b[i].x - cn.y * b[i].y;
    a[i].y = cn.x * b[i].y + cn.y * b[i].x;

}

// Complex pointwise multiplication
__global__ void gpu_applyCCMask(float* ref, float* mask, float* output, float volumeSize, dim3 dim, dim3 halfDim)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;


    if (i >= dim.x || j >= dim.y || k >= dim.z)
        return;

    size_t input_index = i + j*dim.x + k*dim.x*dim.y;
    //size_t output_index = (k + halfDim.z) % dim.z + ((i + halfDim.x) % dim.x)*dim.z + ((j + halfDim.y) % dim.y)*dim.x*dim.z;
	size_t output_index = (i + halfDim.x) % dim.x + ((j + halfDim.y) % dim.y)*dim.z + ((k + halfDim.z) % dim.z)*dim.x*dim.z;
	
    output[output_index] = ref[input_index] * mask[output_index] / volumeSize;
}

// Complex pointwise masking
__global__ void complexMultiplicationDDD(cufftComplex *a, float *mask, dim3 dim)
{
    //blockDim.x * gridDim.x;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= (dim.z * dim.y * dim.x))
        return;

    if (tid == 0)
    {
        a[tid].x = 0.0f;
        a[tid].y = 0.0f;
    }
    else
    {
        a[tid].x = a[tid].x * mask[tid];
        a[tid].y = a[tid].y * mask[tid];
    }
}

// Complex pointwise multiplication with a constant value
__global__ void gpu_complexMultiplication(cufftComplex *a, float factor, size_t elements)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= elements)
        return;

    a[tid].x = a[tid].x * factor;
    a[tid].y = a[tid].y * factor;
}

// Complex pointwise multiplication with a constant value
__global__ void gpu_multiplication(float *a, float* b, size_t elements)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= elements)
        return;

    a[tid] = a[tid] * b[tid];
}

// Pointwise division
__global__ void division(float* a, float div, size_t elements)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= elements)
        return;

    a[tid] = a[tid] / div;
}

// Pointwise division
__global__ void gpu_subtract(float* a, float value, size_t elements)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= elements)
        return;
    
    if (a[tid]!=0.0f)
        a[tid] = a[tid] - value;
}

// Pointwise multiplication
__global__ void gpu_pointwiseMult(float* a, float* b, float* c, size_t elements)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= elements)
        return;

    c[tid] = a[tid] * b[tid];
}

// compute normalized square
__global__ void gpu_normalizedSquare(float *a, float factor, size_t elements)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= elements)
        return;

    a[tid] = a[tid] * a[tid] / factor;
}

// create complex conjugate
__global__ void gpu_complexConjugate(cufftComplex* a, cufftComplex* b, size_t elements)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid >= elements)
        return;

    a[tid].y = -a[tid].y;
    b[tid].y = -b[tid].y;
}

__device__  float binarizeValue(float value)
{
    float binValue;

    if (value >= 0.5f)
        binValue = 1.0f;
    else
        binValue = 0.0f;

    return binValue;
}


__device__ float interpolation(float* inputVolume, size_t inputIndex, float vx1, float vx2, float vy1, float vy2, float vz1, float vz2, size_t x, size_t y)
{
    float fb = inputVolume[inputIndex] + (inputVolume[inputIndex + 1] - inputVolume[inputIndex]) * vx2;
    float ft = inputVolume[inputIndex + x] * vx1 + inputVolume[inputIndex + x + 1] * vx2;
    float rb = inputVolume[inputIndex + x*y] * vx1 + inputVolume[inputIndex + x*y + 1] * vx2;
    float rt = inputVolume[inputIndex + x + x*y] * vx1 + inputVolume[inputIndex + x + x*y + 1] * vx2;

    float finalValue = (fb * vy1 + ft * vy2) * vz1 + (rb * vy1 + rt * vy2) * vz2;

    return finalValue;
}

__global__
void gpu_rotateVolume(float* volume, float* output, float* mask, float*maskOutput, float* ccMask, float*ccMaskOutput, size_t x, size_t y, size_t z, float4 quat, unsigned int rotationType)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= x || j >= y || k >= z)
        return;

    float cx = floor((float)x / 2.0f);
    float cy = floor((float)y / 2.0f);
    float cz = floor((float)z / 2.0f);

    long sxy = x * y;

    float3 rm1, rm2, rm3;
    float qi = quat.x;
    float qj = quat.y;
    float qk = quat.z;
    float w = quat.w;
    rm1.x = 1.0f - 2.0f*(qj*qj + qk*qk);
    rm1.y = 2.0f*(qi*qj + w*qk);
    rm1.z = 2.0f*(qi*qk - w*qj);

    rm2.x = 2.0f*(qi*qj - w*qk);
    rm2.y = 1.0f - 2.0f*(qi*qi + qk*qk);
    rm2.z = 2.0f*(qj*qk + w*qi);

    rm3.x = 2.0f*(qi*qk + w*qj);
    rm3.y = 2.0f*(qj*qk - w*qi);
    rm3.z = 1.0f - 2.0f*(qi*qi + qj*qj);

    long pi = i - cx;
    long pj = j - cy;
    long pk = k - cz;

    size_t outputIndex = i + j * x + k * sxy;

    /* transformation of coordinates */
    float rx = cx + rm1.x * pi + rm1.y * pj + rm1.z * pk;
    float ry = cy + rm2.x * pi + rm2.y * pj + rm2.z * pk;
    float rz = cz + rm3.x * pi + rm3.y * pj + rm3.z * pk;

    
    float floorx = floor(rx);

    if (rx < 0 || floorx >= (x - 1))
    {
        if (rotationType != 1 && rotationType != 4)
        {
            output[outputIndex] = 0.0f;
            maskOutput[outputIndex] = 0.0f;
            ccMaskOutput[outputIndex] = 0.0f;
        }
        return;
    }
    
    float floory = floor(ry);
    if (ry < 0 || floory >= (y - 1))
    {
        if (rotationType != 1 && rotationType != 4)
        {
            output[outputIndex] = 0.0f;
            maskOutput[outputIndex] = 0.0f;
            ccMaskOutput[outputIndex] = 0.0f;
        }
        return;
    }
    
    float floorz = floor(rz);
    if (rz < 0 || floorz >= (z - 1))
    {
        if (rotationType != 1 && rotationType != 4)
        {
            output[outputIndex] = 0.0f;
            maskOutput[outputIndex] = 0.0f;
            ccMaskOutput[outputIndex] = 0.0f;
        }
        return;
    }

    /* Interpolation */

    float vx2 = rx - floorx;
    float vx1 = 1 - vx2;

    float vy2 = ry - floory;
    float vy1 = 1 - vy2;

    float vz2 = rz - floorz;
    float vz1 = 1 - vz2;

    /* the following section detects border pixels to avoid exceeding dimensions */
    size_t inputIndex = floorx + floory * x + floorz * sxy;

    if (rotationType == 1)
    {
        output[outputIndex] += interpolation(volume, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2,x,y);
        maskOutput[outputIndex] += binarizeValue(interpolation(mask, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2, x, y));
    }
	else if (rotationType == 4)
    {
        output[outputIndex] += interpolation(volume, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2,x,y);
        maskOutput[outputIndex] += interpolation(mask, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2, x, y);
    }
    else
    {
        output[outputIndex] = interpolation(volume, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2, x, y);
        ccMaskOutput[outputIndex] = interpolation(ccMask, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2, x, y);
    }

    if (rotationType == 3)
    {
        maskOutput[outputIndex] = binarizeValue(interpolation(mask, inputIndex, vx1, vx2, vy1, vy2, vz1, vz2, x, y));
    }
}


__global__
void gpu_rotateTexture(float* output, float* mask, float* ccmask, float4 quat, size_t x, size_t y, size_t z, unsigned int type)
{

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;

    if(i>=x || j>=y || k>=z)
        return;

    float cx = floor((float)x / 2.0f);
    float cy = floor((float)y / 2.0f);
    float cz = floor((float)z / 2.0f);

    // Rotation matrix computation

    float3 rm1, rm2, rm3;
    float qi = quat.x;
    float qj = quat.y;
    float qk = quat.z;
    float w = quat.w;
    rm1.x = 1.0f - 2.0f*(qj*qj + qk*qk);
    rm1.y = 2.0f*(qi*qj + w*qk);
    rm1.z = 2.0f*(qi*qk - w*qj);

    rm2.x = 2.0f*(qi*qj - w*qk);
    rm2.y = 1.0f - 2.0f*(qi*qi + qk*qk);
    rm2.z = 2.0f*(qj*qk + w*qi);

    rm3.x = 2.0f*(qi*qk + w*qj);
    rm3.y = 2.0f*(qj*qk - w*qi);
    rm3.z = 1.0f - 2.0f*(qi*qi + qj*qj);


    long pi = i - cx;
    long pj = j - cy;
    long pk = k - cz;

    // transformation of coordinates 
    float rx = cx + rm1.x * pi + rm1.y * pj + rm1.z * pk + 0.5f;
    float ry = cy + rm2.x * pi + rm2.y * pj + rm2.z * pk + 0.5f;
    float rz = cz + rm3.x * pi + rm3.y * pj + rm3.z * pk + 0.5f;

    // read from 3D texture
    if (type == 1)
    {
        output[i + j*x + k*x*y] += tex3D(texSubtomo, rx, ry, rz);
        mask[i + j*x + k*x*y] += binarizeValue(tex3D(texWedge, rx, ry, rz));
    }
    else if (type == 4)
    {
        output[i + j*x + k*x*y] += tex3D(texSubtomo, rx, ry, rz);
        mask[i + j*x + k*x*y] += tex3D(texWedge, rx, ry, rz);
    }
    else
    {
        output[i + j*x + k*x*y] = tex3D(texRef, rx, ry, rz);
        ccmask[i + j*x + k*x*y] = tex3D(texCCMask, rx, ry, rz);
    }

    if (type == 3)
    {
        mask[i + j*x + k*x*y] = binarizeValue(tex3D(texMask, rx, ry, rz));
    }
}

void GPURoutines::fft2DR2C(float* input, float* output_real, float* output_img, size_t x, size_t y)
{
    size_t ty = floor(y / 2) + 1;

    size_t numbytes = x*y*sizeof(float);
    size_t numel = x*y;

    float* fft_in;
    cudaMalloc((void**)&fft_in, sizeof(float)*numel);

    //cudaMalloc((void**)&dev_in, numbytes);
    cudaMemcpy(fft_in, input, numbytes, cudaMemcpyHostToDevice);

    cufftComplex* fft_out;
    cudaMalloc((void**)&fft_out, sizeof(cufftComplex)*x*ty);

    // CUFFT plan simple API
    cufftHandle plan;
    cufftPlan2d(&plan, x, y, CUFFT_R2C);

    // Transform signal and kernel
    cufftExecR2C(plan, fft_in, fft_out);

    // Copy device memory to host
    vector<float2> outcome(x*ty);

    cudaMemcpy(&outcome[0], fft_out, x*ty*sizeof(float2), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < x*ty; i++)
    {
        output_real[i] = outcome[i].x;
        output_img[i] = outcome[i].y;
    }

}

void GPURoutines::fft2DC2R(float* input_real, float* input_img, float* output, size_t x, size_t y)
{
    size_t ty = floor(y / 2) + 1;
    vector<float2> complex_in(x*ty);

    for (size_t i = 0; i < x*ty; i++)
    {
        complex_in[i].x = input_real[i];
        complex_in[i].y = input_img[i];
    }

    cufftComplex* fft_in;
    cudaMalloc((void**)&fft_in, sizeof(cufftComplex)*x*ty);

    cudaMemcpy(fft_in, &complex_in[0], sizeof(cufftComplex)*x*ty, cudaMemcpyHostToDevice);

    cufftHandle planBack;
    cufftPlan2d(&planBack, x, y, CUFFT_C2R);

    float* fft_out;
    cudaMalloc((void**)&fft_out, sizeof(float)*x*y);

    // Transform signal and kernel
    cufftExecC2R(planBack, fft_in, fft_out);

    cudaMemcpy(&output[0], fft_out, x*y*sizeof(float), cudaMemcpyDeviceToHost);
}



void GPURoutines::fft3DR2C(float* input, float* output_real, float* output_img, size_t x, size_t y, size_t z)
{
    size_t tz = floor(z / 2) + 1;

    size_t numbytes = x*y*z*sizeof(float);
    size_t numel = x*y*z;

    float* fft_in;
    cudaMalloc((void**)&fft_in, sizeof(float)*numel);

    //cudaMalloc((void**)&dev_in, numbytes);
    cudaMemcpy(fft_in, input, numbytes, cudaMemcpyHostToDevice);

    cufftComplex* fft_out;
    cudaMalloc((void**)&fft_out, sizeof(cufftComplex)*x*y*tz);

    // CUFFT plan simple API
    cufftHandle plan;
    cufftPlan3d(&plan, x, y, z, CUFFT_R2C);

    // Transform signal and kernel
    cufftExecR2C(plan, fft_in, fft_out);

    // Copy device memory to host
    vector<float2> outcome(x*y*tz);

    cudaMemcpy(&outcome[0], fft_out, x*y*tz*sizeof(float2), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < x*y*tz; i++)
    {
        output_real[i] = outcome[i].x;
        output_img[i] = outcome[i].y;
    }
}


void GPURoutines::fft3DC2R(float* input_real, float* input_img, float* output, size_t x, size_t y, size_t z)
{
    size_t tz = floor(z / 2) + 1;
    vector<float2> complex_in(x*y*tz);

    for (size_t i = 0; i < x*y*tz; i++)
    {
        complex_in[i].x = input_real[i];
        complex_in[i].y = input_img[i];
    }

    cufftComplex* fft_in;
    cudaMalloc((void**)&fft_in, sizeof(cufftComplex)*x*y*tz);

    cudaMemcpy(fft_in, &complex_in[0], sizeof(cufftComplex)*x*y*tz, cudaMemcpyHostToDevice);

    cufftHandle planBack;
    cufftPlan3d(&planBack, x, y, z, CUFFT_C2R);

    float* fft_out;
    cudaMalloc((void**)&fft_out, sizeof(float)*x*y*z);

    // Transform signal and kernel
    cufftExecC2R(planBack, fft_in, fft_out);

    cudaMemcpy(&output[0], fft_out, x*y*z*sizeof(float), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < x*y*z; i++)
    {
        output[i] = output[i] / (float)(x*y*z);
    }
}

void GPURoutines::fft3DC2C(float* input_real, float* input_img, float* output, size_t x, size_t y, size_t z)
{
    vector<float2> complex_in(x*y*z);

    for (size_t i = 0; i < x*y*z; i++)
    {
        complex_in[i].x = input_real[i];
        complex_in[i].y = input_img[i];
    }

    cufftComplex* fft_in;
    cudaMalloc((void**)&fft_in, sizeof(cufftComplex)*x*y*z);

    cudaMemcpy(fft_in, &complex_in[0], sizeof(cufftComplex)*x*y*z, cudaMemcpyHostToDevice);

    cufftHandle planBack;
    cufftPlan3d(&planBack, x, y, z, CUFFT_C2C);

    cufftComplex* fft_out;
    cudaMalloc((void**)&fft_out, sizeof(cufftComplex)*x*y*z);

    // Transform signal and kernel
    cufftExecC2C(planBack, fft_in, fft_out, CUFFT_INVERSE);

    vector<float2> tempOut(x*y*z);
    cudaMemcpy(&tempOut[0], fft_out, x*y*z*sizeof(cufftComplex), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < x*y*z; i++)
    {
        output[i] = tempOut[i].x/(float)(x*y*z);
    }

}

GPURoutines::GPURoutines(size_t dimx, size_t dimy, size_t dimz, bool useRoseman,int processID)
{
    volumeDim = dim3( dimx, dimy, dimz );
    volumeSize = dimx*dimy*dimz;

    fftz = floor((float)volumeDim.z / 2) + 1;
    fftSize = dimx*dimy*fftz;

    // Create 3D array
    volumeExtent = make_cudaExtent(dimx, dimy, dimz);

    // Create a channel description
    channelDesc = cudaCreateChannelDesc<float>();

    threadsPerBlock1D = (8*8*8);  // 576 threads
    numBlocks1D = ceil(volumeSize / (float)threadsPerBlock1D.x);

    threadsPerBlock = { 8, 8, 8 };  // 512 threads
    numBlocks = dim3(ceil(dimx / (float)threadsPerBlock.x), ceil(dimy / (float)threadsPerBlock.y), ceil(dimz / (float)threadsPerBlock.z));

    threadsPerBlockFFT = { 8, 8, 8 };  // 576 threads
    numBlocksFFT = dim3(ceil(dimx / (float)threadsPerBlockFFT.x), ceil(dimy / (float)threadsPerBlockFFT.y), ceil(dimz / (float)threadsPerBlockFFT.z));

    subtomoAllocated = false;

    useRosemanCC = useRoseman;

    if (useRosemanCC)
        rotationType = 3;
    else
        rotationType = 2;
        
    textureAllocated = false;
  	
	distributeDevices(processID);
}

void GPURoutines::distributeDevices(int processID)
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    
    int deviceID=processID%deviceCount;
    cudaSetDevice(deviceID);
  
    cout << "Process #" << processID << " got card #" << deviceID << endl;   
}

GPURoutines::~GPURoutines()
{
    cudaFree(shiftFilter);
    
    if (subtomoAllocated)
        cudaFree(fftSubtomoInNew);

    cudaFree(fftSubtomoOut);

    cudaFree(deviceSubtomoRot);
    cudaFree(deviceWedgeRot);
    cudaFree(deviceCCMaskRot);
    cudaFree(deviceRefRot);
   
    cufftDestroy(planR2C);
    cufftDestroy(planC2R);

 //   cudaFree(rrrRef);
 //   cudaFree(rrrMask);
 //   cudaFree(rrrCCMask);

    if (useRosemanCC)
    {
        cudaFree(deviceMaskRot);
        cudaFree(conjugateSquare);
        cudaFree(intensitySubtomo);
        cudaFree(meanSubtomo);
        cudaFree(meanMasked);
        cudaFree(numerator);
    }

    if(textureAllocated)
    {
        cudaFreeArray(deviceRef);
        cudaFreeArray(deviceCCMask);
        cudaFreeArray(deviceMask);
    }
}


void GPURoutines::prepareReferenceTexture(vector<float>& reference, vector<float>& mask, vector<float>& ccMask, string filterModeName)
{
    cudaTextureFilterMode filterMode;

    if (filterModeName == "linear")
        filterMode = cudaFilterModeLinear;
    else
        filterMode = cudaFilterModePoint;

    // Allocate memory on device
    cudaMalloc3DArray(&deviceRef, &channelDesc, volumeExtent);
    cudaMalloc3DArray(&deviceCCMask, &channelDesc, volumeExtent);
    cudaMalloc3DArray(&deviceMask, &channelDesc, volumeExtent);
    
    textureAllocated = true;
    
    // --- Set texture parameters
    texRef.normalized = false;                          // access with normalized texture coordinates
    texRef.filterMode = filterMode;                     // trilinear or no interpolation
    texRef.addressMode[0] = cudaAddressModeBorder;      // wrap texture coordinates
    texRef.addressMode[1] = cudaAddressModeBorder;
    texRef.addressMode[2] = cudaAddressModeBorder;

    texCCMask.normalized = false;                       // access with normalized texture coordinates
    texCCMask.filterMode = filterMode;                     // trilinear or no interpolation
    texCCMask.addressMode[0] = cudaAddressModeBorder;      // wrap texture coordinates
    texCCMask.addressMode[1] = cudaAddressModeBorder;
    texCCMask.addressMode[2] = cudaAddressModeBorder;

    // Set copy parameters for 3D array
    cudaMemcpy3DParms copyParamsReference = { 0 };
    copyParamsReference.dstArray = deviceRef;
    copyParamsReference.extent = volumeExtent;
    copyParamsReference.kind = cudaMemcpyHostToDevice;
    copyParamsReference.srcPtr = make_cudaPitchedPtr((void *)&reference[0], volumeExtent.width*sizeof(float), volumeExtent.width, volumeExtent.height);

    // Copy 3D array from host to device
    gpuErrchk(cudaMemcpy3D(&copyParamsReference));

    // Bind array to 3D texture
    gpuErrchk(cudaBindTextureToArray(texRef, deviceRef, channelDesc));

    // Allocate memory for rotated reference
    gpuErrchk(cudaMalloc((void **)&deviceRefRot, volumeSize*sizeof(float)));
    
    // Set copy parameters for 3D array
    copyParamsReference.dstArray = deviceCCMask;
    copyParamsReference.srcPtr = make_cudaPitchedPtr((void *)&ccMask[0], volumeExtent.width*sizeof(float), volumeExtent.width, volumeExtent.height);

    // Copy 3D array from host to device
    gpuErrchk(cudaMemcpy3D(&copyParamsReference));

    // Bind array to 3D texture
    gpuErrchk(cudaBindTextureToArray(texCCMask, deviceCCMask, channelDesc));
    gpuErrchk(cudaMalloc((void **)&deviceCCMaskRot, volumeSize*sizeof(float)));

   
    if (useRosemanCC)
    {
        // Set copy parameters for 3D array
        copyParamsReference.dstArray = deviceMask;
        copyParamsReference.srcPtr = make_cudaPitchedPtr((void *)&mask[0], volumeExtent.width*sizeof(float), volumeExtent.width, volumeExtent.height);

        // Copy 3D array from host to device
        gpuErrchk(cudaMemcpy3D(&copyParamsReference));

        // Bind array to 3D texture
        gpuErrchk(cudaBindTextureToArray(texMask, deviceMask, channelDesc));
        gpuErrchk(cudaMalloc((void **)&deviceMaskRot, volumeSize*sizeof(float)));
    }

  //  gpuErrchk(cudaMalloc((void **)&deviceMaskRot, volumeSize*sizeof(float)));

  //  gpuErrchk(cudaMalloc((void **)&rrrRef, volumeSize*sizeof(float)));
  //  gpuErrchk(cudaMemcpy(rrrRef, &reference[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));
  
  //  gpuErrchk(cudaMalloc((void **)&rrrMask, volumeSize*sizeof(float)));
  //  gpuErrchk(cudaMemcpy(rrrMask, &mask[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));

  //  gpuErrchk(cudaMalloc((void **)&rrrCCMask, volumeSize*sizeof(float)));
   // gpuErrchk(cudaMemcpy(rrrCCMask, &ccMask[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));

    // Just for now to prevent zero rotation
    //gpuErrchk(cudaMemcpy(deviceRefRot, &reference[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));
    //gpuErrchk(cudaMemcpy(deviceCCMaskRot, &ccMask[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));
    //gpuErrchk(cudaMemcpy(deviceMaskRot, &mask[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));
}

void GPURoutines::prepareSubtomoTexture(string filterModeName)
{
    cudaTextureFilterMode filterMode;

    if (filterModeName == "linear")
        filterMode = cudaFilterModeLinear;
    else
        filterMode = cudaFilterModePoint;
    
    // Allocate memory on device
    cudaMalloc3DArray(&deviceSubtomo, &channelDesc, volumeExtent);
    cudaMalloc3DArray(&deviceWedge, &channelDesc, volumeExtent);


    // --- Set texture parameters
    texSubtomo.normalized = false;                      // access with normalized texture coordinates
    texSubtomo.filterMode = filterMode;                 // trilinear interpolation
    texSubtomo.addressMode[0] = cudaAddressModeBorder;  // wrap texture coordinates
    texSubtomo.addressMode[1] = cudaAddressModeBorder;
    texSubtomo.addressMode[2] = cudaAddressModeBorder;

    // --- Set texture parameters
    texWedge.normalized = false;                      // access with normalized texture coordinates
    texWedge.filterMode = filterMode;                 // trilinear interpolation
    texWedge.addressMode[0] = cudaAddressModeBorder;  // wrap texture coordinates
    texWedge.addressMode[1] = cudaAddressModeBorder;
    texWedge.addressMode[2] = cudaAddressModeBorder;

    // Set copy parameters for 3D array
    copyParamsSubtomo = { 0 };
    copyParamsSubtomo.dstArray = deviceSubtomo;
    copyParamsSubtomo.extent = volumeExtent;
    copyParamsSubtomo.kind = cudaMemcpyDeviceToDevice;
    
    // Set copy parameters for 3D array
    copyParamsWedge = { 0 };
    copyParamsWedge.dstArray = deviceWedge;
    copyParamsWedge.extent = volumeExtent;
    copyParamsWedge.kind = cudaMemcpyHostToDevice;
    
}

void GPURoutines::copyAndBindWedgeTexture(vector<float>& wedgeMask)
{
    //double minAngleRad = minAngle*CUDART_PI_F / 180.0;
    //double maxAngleRad = maxAngle*CUDART_PI_F / 180.0;

    //double tan_min = tan((-CUDART_PI_F / 2.0) - minAngleRad);
    //double tan_max = tan(CUDART_PI_F / 2.0 - maxAngleRad);

    //wedgeMask.resize((size_t)dim*dim*dim);
    //fill(wedgeMask.begin(), wedgeMask.end(), 1.0f);

    
    
    copyParamsWedge.srcPtr = make_cudaPitchedPtr((void *)&wedgeMask[0], volumeExtent.width*sizeof(float), volumeExtent.width, volumeExtent.height);

    // Copy 3D array from host to device
    gpuErrchk(cudaMemcpy3D(&copyParamsWedge));

    // Bind array to 3D texture
    gpuErrchk(cudaBindTextureToArray(texWedge, deviceWedge, channelDesc));
}

void GPURoutines::rotateVectorIndices(vector<int>& ids, size_t vectorSize, int direction)
{
    size_t shiftStart;
    
    if (direction==-1)
        shiftStart = (size_t)ceil((float)vectorSize / 2.0f);
    else
        shiftStart = (size_t)floor((float)vectorSize / 2.0f);

    for (int i = 0; i < vectorSize; i++)
        ids.push_back(i);

    rotate(ids.begin(), ids.begin() + shiftStart, ids.end());
}

void GPURoutines::allocateDeviceMemory()
{
    // allocate memory for shift filter
    gpuErrchk(cudaMalloc((void **)&shiftFilter, fftSize*sizeof(cufftComplex)));
    gpuErrchk(cudaMemset(shiftFilter, 0.0, fftSize*sizeof(cufftComplex)));
   
  //  gpuErrchk(cudaMalloc((void **)&deviceSubtomoRot, volumeSize*sizeof(float)));
    gpuErrchk(cudaMalloc((void **)&deviceWedgeRot, volumeSize*sizeof(float)));

    gpuErrchk(cudaMemset(deviceWedgeRot, 0.0f, volumeSize*sizeof(float)));
  //  gpuErrchk(cudaMemset(deviceSubtomoRot, 0.0f, volumeSize*sizeof(float)));

    cufftPlan3d(&planR2C, volumeDim.x, volumeDim.y, volumeDim.z, CUFFT_R2C);
    cufftPlan3d(&planC2R, volumeDim.x, volumeDim.y, volumeDim.z, CUFFT_C2R);

 
    // allocate memory for subtomogram for FFT output
    gpuErrchk(cudaMalloc((void **)&fftSubtomoOut, fftSize*sizeof(cufftComplex)));

    gpuErrchk(cudaMalloc((void **)&deviceSubtomoRot, volumeSize*sizeof(float)));
    gpuErrchk(cudaMemset(deviceSubtomoRot, 0.0f, volumeSize*sizeof(float)));

    cudaDeviceSynchronize();

    if (useRosemanCC)
    {
        gpuErrchk(cudaMalloc((void **)&numerator, volumeSize*sizeof(float)));
        gpuErrchk(cudaMalloc((void **)&meanSubtomo, volumeSize*sizeof(float)));
        gpuErrchk(cudaMalloc((void **)&intensitySubtomo, volumeSize*sizeof(float)));
        gpuErrchk(cudaMalloc((void **)&meanMasked, fftSize*sizeof(cufftComplex)));
    }
}

void GPURoutines::shiftSubtomogram(vector<float>& subtomo, float shiftX, float shiftY, float shiftZ)
{
    float3 shift = { shiftX, shiftY, shiftZ };
   
    gpuErrchk(cudaMalloc((void **)&fftSubtomoInNew, volumeSize*sizeof(float)));
    gpuErrchk(cudaMemcpy(fftSubtomoInNew, &subtomo[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));
    subtomoAllocated = true;

    cudaDeviceSynchronize();

    gpu_generateShiftFilter << <numBlocks, threadsPerBlock >> > (shiftFilter, volumeDim.x, volumeDim.y, volumeDim.z, shift);
   
    // Transform signal and kernel
    cufftExecR2C(planR2C, fftSubtomoInNew, fftSubtomoOut);
    cudaDeviceSynchronize();
    
    complexMultiplication << <numBlocks1D, threadsPerBlock1D >> >(fftSubtomoOut, shiftFilter, fftSize);

    cudaDeviceSynchronize();
 
    // Transform signal and kernel
    cufftExecC2R(planC2R, fftSubtomoOut, fftSubtomoInNew);
    cudaDeviceSynchronize();

    division << <numBlocks1D, threadsPerBlock1D >> > (fftSubtomoInNew, volumeSize, volumeSize);

}

double GPURoutines::sumArray(vector<float>& data)
{
    float* deviceArray;
    cudaMalloc((void**)&deviceArray, sizeof(float)*data.size());
    cudaMemcpy(deviceArray, &data[0], sizeof(float)*data.size(), cudaMemcpyHostToDevice);

    double sum = sumArray(deviceArray, data.size());

    cudaFree(deviceArray);

    return sum;
}

double GPURoutines::computeNormFactor(vector<float2>& data)
{
    cufftComplex* deviceArray;
    cudaMalloc((void**)&deviceArray, sizeof(cufftComplex)*data.size());
    cudaMemcpy(deviceArray, &data[0], sizeof(cufftComplex)*data.size(), cudaMemcpyHostToDevice);

    double nf = computeNormFactor(deviceArray, data.size());

    cudaFree(deviceArray);

    return nf;
}

double GPURoutines::sumArray(float* deviceArray, size_t arraySize, bool squared)
{
    unsigned int threads = 512;
    unsigned int blocks = 2;

    double* deviceOut;
    cudaMalloc((void**)&deviceOut, sizeof(double)*threads*blocks);

    if (squared)
        gpu_sumArraySquared << <blocks, threads >> >(deviceArray, deviceOut, arraySize);
    else
        gpu_sumArray << <blocks, threads >> >(deviceArray, deviceOut, arraySize);

    vector<double> out(threads*blocks);
    cudaMemcpy(&out[0], deviceOut, sizeof(double)*threads*blocks, cudaMemcpyDeviceToHost);

    double globalSum = 0;
    for (size_t i = 0; i < threads*blocks; i++)
        globalSum += out[i];

    cudaFree(deviceOut);
    return globalSum;
}

double GPURoutines::computeNormFactor(cufftComplex* deviceArray, size_t arraySize)
{

    float* conjugateValues;
    cudaMalloc((void**)&conjugateValues, sizeof(float)*fftSize);
    
    dim3 fftDim(fftz, volumeDim.y, volumeDim.z);    
    gpu_prepareNormValue<< <numBlocks, threadsPerBlock >> >(deviceArray,conjugateValues,fftDim);
    
    
    unsigned int threads = 512;
    unsigned int blocks = 1;

    double* deviceOut;
    cudaMalloc((void**)&deviceOut, sizeof(double)*threads*blocks);

    gpu_sumArray << <blocks, threads >> >(conjugateValues, deviceOut, arraySize);

    vector<double> out(threads*blocks);
    cudaMemcpy(&out[0], deviceOut, sizeof(double)*threads*blocks, cudaMemcpyDeviceToHost);

    double globalSum = 0;
    for (size_t i = 0; i < threads*blocks; i++)
    {
        globalSum = globalSum + out[i];
    }
    
    cudaFree(conjugateValues);
	cudaFree(deviceOut);
    return globalSum;
}

void GPURoutines::computeConjugateOfSquare()
{
   // cufftComplex* conjugateSquare;
    gpuErrchk(cudaMalloc((void **)&conjugateSquare, fftSize*sizeof(cufftComplex)));

    float* fft_out;
    cudaMalloc((void**)&fft_out, sizeof(float)*volumeSize);

    cufftHandle planBack;
    cufftPlan3d(&planBack, volumeDim.x, volumeDim.y, volumeDim.z, CUFFT_C2R);

   // cudaDeviceSynchronize();

    // Transform signal and kernel
    cufftExecC2R(planBack, fftSubtomoOut, fft_out);

    float normFactor = volumeSize*volumeSize;

    gpu_normalizedSquare << <numBlocks1D, threadsPerBlock1D >> >(fft_out, normFactor, volumeSize);

   // cudaDeviceSynchronize();

    cufftHandle planForward;
    cufftPlan3d(&planForward, volumeDim.x, volumeDim.y, volumeDim.z, CUFFT_R2C);

    cufftExecR2C(planForward, fft_out, conjugateSquare);

    cufftDestroy(planForward);
    cufftDestroy(planBack);
	cudaFree(fft_out);
}

void GPURoutines::maskSubtomogram(vector<float>& subtomogram, vector<float>& mask)
{
    // Create a stream for async copy
    //cudaStream_t stream1;
   // cudaStreamCreate(&stream1);

    //copy mask
    float* dMask;
    gpuErrchk(cudaMalloc((void **)&dMask, fftSize*sizeof(float)));
    // copy subtomogram
    float* fftIn;
    gpuErrchk(cudaMalloc((void **)&fftIn, volumeSize*sizeof(float)));
    
    
    gpuErrchk(cudaMemcpy(dMask, &mask[0], fftSize*sizeof(float), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(fftIn, &subtomogram[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));

    // Calculate shifted fourier transform of rotated reference
    cufftExecR2C(planR2C, fftIn, fftSubtomoOut);

    //cudaStreamSynchronize(stream1);
    dim3 fftDim(volumeDim.x, volumeDim.y, fftz);
    complexMultiplicationDDD << <numBlocks1D, threadsPerBlock1D >> > (fftSubtomoOut, dMask, fftDim);

    if (useRosemanCC)
    {
        computeConjugateOfSquare();
        gpu_complexConjugate << <numBlocks1D, threadsPerBlock1D >> > (fftSubtomoOut, conjugateSquare, fftSize);
    }
    else
    {
        double normFactor = computeNormFactor(fftSubtomoOut, fftSize);

        normFactor = volumeSize / sqrt(normFactor);
        // Calculate cross correlation and apply rotated ccmask
        gpu_complexMultiplication << <numBlocks1D, threadsPerBlock1D >> >(fftSubtomoOut, normFactor, fftSize);
    }

    cudaFree(fftIn);
    cudaFree(dMask);
    
    //cudaStreamDestroy(stream1);
}


void GPURoutines::computeCC(vector<float>& ccVolume, vector<float>& shiftedWedgeMask)
{
    
    cufftComplex* fftRef;
    gpuErrchk(cudaMalloc((void **)&fftRef, fftSize*sizeof(cufftComplex)));

    float* deviceCCVolume;
    gpuErrchk(cudaMalloc((void **)&deviceCCVolume, volumeSize*sizeof(float)));

  //  vector<float> rotRef(volumeSize);
  //  gpuErrchk(cudaMemcpy(&rotRef[0], deviceRefRot, volumeSize*sizeof(float), cudaMemcpyDeviceToHost));
    
    // Calculate shifted fourier transform of rotated reference
    cufftExecR2C(planR2C, deviceRefRot, fftRef);

    // Allocate device memory for the shifted wedge mask
    float* deviceShiftedWedge;
    gpuErrchk(cudaMalloc((void **)&deviceShiftedWedge, fftSize*sizeof(float)));
    gpuErrchk(cudaMemcpy(deviceShiftedWedge, &shiftedWedgeMask[0], fftSize*sizeof(float), cudaMemcpyHostToDevice));

    // Apply bandpass filtered wedge
    dim3 fftDim(volumeDim.x, volumeDim.y, fftz);
    complexMultiplicationDDD << < numBlocks1D, threadsPerBlock1D >> > (fftRef, deviceShiftedWedge, fftDim);
    
    cudaFree(deviceShiftedWedge);

    if (useRosemanCC)
    {

        // Get real space version of masked reference
        cufftExecC2R(planC2R, fftRef, deviceCCVolume);
        division << <numBlocks1D, threadsPerBlock1D >> > (deviceCCVolume, volumeSize, volumeSize);

        // sum rotated Mask
        double maskSum = sumArray(deviceMaskRot, volumeSize);

        // mask and normalize reference
        float* deviceMaskedRef;
        gpuErrchk(cudaMalloc((void **)&deviceMaskedRef, volumeSize*sizeof(float)));
        gpu_pointwiseMult << < numBlocks1D, threadsPerBlock1D >> > (deviceMaskRot, deviceCCVolume, deviceMaskedRef, volumeSize);

        double maskedMean = sumArray(deviceMaskedRef, volumeSize);
        maskedMean = maskedMean / maskSum;

        // Normalization factor of references
        gpu_subtract << < numBlocks1D, threadsPerBlock1D >> >(deviceMaskedRef, maskedMean, volumeSize);

        double sigmaRef = sumArray(deviceMaskedRef, volumeSize, true);
        sigmaRef = sqrt(sigmaRef);

        // Fourier transform of masked ref
        cufftComplex* maskedFFT;
        gpuErrchk(cudaMalloc((void **)&maskedFFT, fftSize*sizeof(cufftComplex)));
        cufftExecR2C(planR2C, deviceMaskedRef, maskedFFT);

        cudaFree(deviceMaskedRef);

        // Convolution of masked reference and subtomo
        complexMultiplication << < numBlocks1D, threadsPerBlock1D >> >(maskedFFT, fftSubtomoOut, fftSize);

        cufftExecC2R(planC2R, maskedFFT, numerator);
        division << <numBlocks1D, threadsPerBlock1D >> > (numerator, volumeSize, volumeSize);

        // Fourier transform of mask
        cufftExecR2C(planR2C, deviceMaskRot, maskedFFT);

        // Mean of subtomo under mask
        gpu_complexMultiplication << < numBlocks1D, threadsPerBlock1D >> >(maskedFFT, fftSubtomoOut, meanMasked, fftSize);

        cufftExecC2R(planC2R, meanMasked, meanSubtomo);
        division << <numBlocks1D, threadsPerBlock1D >> > (meanSubtomo, volumeSize, volumeSize);

        // Mean intensity of subtomo under mask
        complexMultiplication << < numBlocks1D, threadsPerBlock1D >> >(maskedFFT, conjugateSquare, fftSize);

        cufftExecC2R(planC2R, maskedFFT, intensitySubtomo);
        division << <numBlocks1D, threadsPerBlock1D >> > (intensitySubtomo, volumeSize, volumeSize);

        // Calculate denominator (of eq 5 in paper)
        gpu_computeFLCC << < numBlocks, threadsPerBlock >> >(deviceCCVolume, numerator, intensitySubtomo, meanSubtomo, maskSum, sigmaRef, volumeDim);
        gpu_multiplication << <numBlocks1D, threadsPerBlock1D >> >(deviceCCVolume, deviceCCMaskRot, volumeSize);
        
        cudaFree(maskedFFT);

    }
    else
    {
        double normFactor = computeNormFactor(fftRef, fftSize);
        normFactor = volumeSize / sqrt(normFactor);

        // Calculate cross correlation and apply rotated ccmask
        gpu_conjugateMultiplication << <numBlocks1D, threadsPerBlock1D >> >(fftRef, fftSubtomoOut, normFactor, fftSize);

        // Calculate shifted fourier transform of rotated reference
        cufftExecC2R(planC2R, fftRef, deviceRefRot);
        division << <numBlocks1D, threadsPerBlock1D >> > (deviceRefRot, volumeSize, volumeSize);

        dim3 halfDim(floor(volumeDim.x / 2.0f), floor(volumeDim.y / 2.0f), floor(volumeDim.z / 2.0f));

        gpu_applyCCMask << <numBlocks, threadsPerBlock >> >(deviceRefRot, deviceCCMaskRot, deviceCCVolume, (float)volumeSize, volumeDim, halfDim);
    }

    gpuErrchk(cudaMemcpy(&ccVolume[0], deviceCCVolume, volumeSize*sizeof(float), cudaMemcpyDeviceToHost));
    cudaFree(deviceCCVolume);
    cudaFree(fftRef);

}

void GPURoutines::rotateTexture(Quaternion& rotation)
{
    float4 q;
    q.x = rotation.i;
    q.y = rotation.j;
    q.z = rotation.k;
    q.w = rotation.w;

    gpu_rotateTexture << <numBlocks, threadsPerBlock >> >(deviceRefRot, deviceMaskRot, deviceCCMaskRot, q, volumeDim.x, volumeDim.y, volumeDim.z, rotationType);
    
   // gpu_rotateVolume << <numBlocks, threadsPerBlock >> >(rrrRef, deviceRefRot, rrrMask, deviceMaskRot, rrrCCMask, deviceCCMaskRot, volumeDim.x, volumeDim.y, volumeDim.z, q, rotationType);

}

void GPURoutines::rotateSubtomogramAndWedge(vector<float>& subtomo, vector<float>& wedge, Quaternion& rotation,bool binarizeMask)
{
    float4 q;
    q.x = rotation.i;
    q.y = rotation.j;
    q.z = rotation.k;
    q.w = rotation.w;

    if (!subtomoAllocated)
    {
        cudaMalloc((void**)&fftSubtomoInNew, volumeSize*sizeof(float));
        gpuErrchk(cudaMemcpy(fftSubtomoInNew, &subtomo[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));
        subtomoAllocated = true;
    }

    copyParamsSubtomo.srcPtr = make_cudaPitchedPtr((void *)fftSubtomoInNew, volumeExtent.width*sizeof(float), volumeExtent.width, volumeExtent.height);

    // Copy 3D array from device to device
    gpuErrchk(cudaMemcpy3D(&copyParamsSubtomo));

    // Bind array to 3D texture
    gpuErrchk(cudaBindTextureToArray(texSubtomo, deviceSubtomo, channelDesc));

   // float* volumeInput;
  //  float* wedgeInput;
  //  float* dummyIn;
  //  float* dummyOut;
  //  gpuErrchk(cudaMalloc((void **)&volumeInput, volumeSize*sizeof(float)));
  //  gpuErrchk(cudaMemcpy(volumeInput, &subtomo[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));

//    gpuErrchk(cudaMalloc((void **)&wedgeInput, volumeSize*sizeof(float)));
  //  gpuErrchk(cudaMemcpy(wedgeInput, &wedge[0], volumeSize*sizeof(float), cudaMemcpyHostToDevice));

	if (binarizeMask)
		gpu_rotateTexture << <numBlocks, threadsPerBlock >> >(deviceSubtomoRot, deviceWedgeRot, deviceCCMaskRot, q, volumeDim.x, volumeDim.y, volumeDim.z, 1);
	else
		gpu_rotateTexture << <numBlocks, threadsPerBlock >> >(deviceSubtomoRot, deviceWedgeRot, deviceCCMaskRot, q, volumeDim.x, volumeDim.y, volumeDim.z, 4);
 
 //   gpu_rotateVolume << <numBlocks, threadsPerBlock >> >(fftSubtomoInNew, deviceSubtomoRot, wedgeInput, deviceWedgeRot, dummyIn, dummyOut, volumeDim.x, volumeDim.y, volumeDim.z, q, 1);

    gpuErrchk(cudaFree(fftSubtomoInNew));
    subtomoAllocated = false;
   // gpuErrchk(cudaFree(wedgeInput));

}

void GPURoutines::getRotatedSubtomoAndWedge(vector<float>& subtomo, vector<float>& wedge)
{
    cudaMemcpy(&subtomo[0], deviceSubtomoRot, volumeSize*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&wedge[0], deviceWedgeRot, volumeSize*sizeof(float), cudaMemcpyDeviceToHost);
}

void GPURoutines::getRotatedReferenceAndCCMask(vector<float>& ref, vector<float>& ccMask)
{
    cudaMemcpy(&ref[0], deviceRefRot, volumeSize*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(&ccMask[0], deviceCCMaskRot, volumeSize*sizeof(float), cudaMemcpyDeviceToHost);
}

void GPURoutines::getFFTResult(vector<float>& result)
{
    cudaMemcpy(&result[0], fftSubtomoInNew, volumeSize*sizeof(float), cudaMemcpyDeviceToHost);
}

void GPURoutines::getMaskedSubtomogram(vector<float>& outReal, vector<float>& outImag)
{
    vector<float2> output(outReal.size());

    cudaMemcpy(&output[0], fftSubtomoOut, outReal.size()*sizeof(cufftComplex), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < outReal.size(); i++)
    {
        outReal[i] = output[i].x;
        outImag[i] = output[i].y;
    }
}

void GPURoutines::getShiftFilter(vector<float>& filterReal, vector<float>& filterImag)
{
    //size_t fftz = floor((float)volumeDim.z / 2) + 1;
   
    vector<float2> filter(volumeDim.x*volumeDim.y*fftz);
    cudaMemcpy(&filter[0], shiftFilter, volumeDim.x*volumeDim.y*fftz*sizeof(cufftComplex), cudaMemcpyDeviceToHost);

    for (size_t i = 0; i < filterReal.size(); i++)
    {
        filterReal[i] = filter[i].x;
        filterImag[i] = filter[i].y;
    }
}

void GPURoutines::getShiftFilter(vector<float2>& filter)
{
    //size_t fftz = floor((float)volumeDim.z / 2) + 1;

    cudaMemcpy(&filter[0], shiftFilter, volumeDim.x*volumeDim.y*fftz*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
}

void GPURoutines::getReferenceAndCCMask(vector<float>& reference, vector<float>& ccmask)
{
    cudaMemcpy3DParms copyParamsReference = { 0 };
    copyParamsReference.dstPtr = make_cudaPitchedPtr((void *)&reference[0], volumeExtent.width*sizeof(float), volumeExtent.width, volumeExtent.height);
    copyParamsReference.extent = volumeExtent;
    copyParamsReference.kind = cudaMemcpyDeviceToHost;
    copyParamsReference.srcArray = deviceRef;

    // Copy 3D array from host to device
    gpuErrchk(cudaMemcpy3D(&copyParamsReference));

    copyParamsReference.dstPtr = make_cudaPitchedPtr((void *)&ccmask[0], volumeExtent.width*sizeof(float), volumeExtent.width, volumeExtent.height);
    copyParamsReference.srcArray = deviceCCMask;

    // Copy 3D array from host to device
    gpuErrchk(cudaMemcpy3D(&copyParamsReference));
}
