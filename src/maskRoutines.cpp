#include "maskRoutines.h"
#include "vectorRoutines.h"

using namespace std;
using namespace vr;

void MaskRoutines::createWedgeMask(vector<float>& wedgeMask, int dim, float minAngle, float maxAngle)
{
    double minAngleRad = minAngle*M_PI / 180.0;
    double maxAngleRad = maxAngle*M_PI / 180.0;

    double tan_min = tan((-M_PI / 2.0) - minAngleRad);
    double tan_max = tan(M_PI / 2.0 - maxAngleRad);

    //wedgeMask.resize((size_t)dim*dim*dim);
    fill(wedgeMask.begin(), wedgeMask.end(), 1.0f);

    float halfDim = floor(dim / 2);
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

void MaskRoutines::createWedgeMask(vector<float>& wedgeMask, vector<float>& shifted, vector<float>& shiftedReduced, int dim, float minAngle, float maxAngle, bool maskEdges, bool shift)
{

    double minAngleRad = minAngle*M_PI / 180.0;
    double maxAngleRad = maxAngle*M_PI / 180.0;

    double tan_min = tan((-M_PI / 2.0) - minAngleRad);
    double tan_max = tan(M_PI / 2.0 - maxAngleRad);

    float radius;

    if (maskEdges)
        radius = floor((dim - 1.0f) / 2.0f);
    else
        radius = 2 * dim; // no masking

    fill(wedgeMask.begin(), wedgeMask.end(), 1.0f);

    // center for the shuffle
    size_t ms = (size_t)floor((float)dim / 2.0f);
    vector<size_t> idx;
    size_t fftz = floor(dim / 2.0f) + 1;

    if (shift)
    {
        fill(shifted.begin(), shifted.end(), 1.0f);

        idx.resize(dim);
        for (size_t i = 0; i < dim; i++)
            idx[i] = i;

        std::rotate(idx.begin(), idx.begin() + ms, idx.end());

        shiftedReduced.resize(dim * dim * fftz);
        fill(shiftedReduced.begin(), shiftedReduced.end(), 1.0f);
    }

    float halfDim = floor(dim / 2);
    for (float z = -halfDim; z < -halfDim + dim ; z++)
    {
        if (z == 0)
        {
            if (shift)// edge mask
            {
                for (float x = -halfDim; x < -halfDim + dim ; x++)
                for (int y = -halfDim; y < -halfDim + dim ; y++)
                {
                    float distance = sqrt(x*x + y*y + z*z);
                    if (distance >= radius)
                    {
                        size_t voxel_index = (size_t)(idx[(x + halfDim)] + idx[(y + halfDim)] * dim + idx[(z + halfDim)] * dim  * dim );
                        shifted[voxel_index] = 0.0f;

                        if (idx[(x + halfDim)]<fftz)
                        {
                            voxel_index = (size_t)(idx[(x + halfDim)] + idx[(y + halfDim)] * fftz + idx[(z + halfDim)] * fftz*dim );
                            shiftedReduced[voxel_index] = 0.0f;
                        }
                    }
                }
            }
            continue;
        }

        for (float x = -halfDim; x < -halfDim + dim ; x++)
        {
            if ((tan_max > x / z) && (tan_min < x / z))
            {
                for (int y = 0; y < dim; y++)
                {
                    wedgeMask[(size_t)((x + halfDim) + y*dim + (z + halfDim)*dim  * dim)] = 0.0f;

                    if (shift)
                    {
                        shifted[(size_t)(idx[(x + halfDim)] + idx[y] * dim  + idx[(z + halfDim)] * dim  * dim )] = 0.0f;

                        if (idx[(x + halfDim)]<fftz)
                        {
                            size_t voxel_index = (size_t)(idx[(x + halfDim)] + idx[y] * fftz + idx[(z + halfDim)] * fftz*dim );
                            shiftedReduced[voxel_index] = 0.0f;
                        }
                    }
                }
            }
            else if (shift)// edge mask
            {
                for (int y = -halfDim; y < -halfDim + dim ; y++)
                {
                    float distance = sqrt(x*x + y*y + z*z);
                    if (distance >= radius)
                    {
                        size_t voxel_index = (size_t)(idx[(x + halfDim)] + idx[(y + halfDim)] * dim  + idx[(z + halfDim)] * dim * dim );
                        shifted[voxel_index] = 0.0f;

                        if (idx[(x + halfDim)]<fftz)
                        {
                            voxel_index = (size_t)(idx[(x + halfDim)] + idx[(y + halfDim)] * fftz + idx[(z + halfDim)] * fftz*dim );
                            shiftedReduced[voxel_index] = 0.0f;
                        }
                    }
                }
            }
        }
    }
}


vector<float> MaskRoutines::createSphericalMask(size_t x, size_t y, size_t z, float radius, float sigma)
{
    float cx = floor((float)x / 2.0f);
    float cy = floor((float)y / 2.0f);
    float cz = floor((float)z / 2.0f);

    vector<float> mask;
    mask.resize(x*y*z);
    fill(mask.begin(), mask.end(), 0.0f);

    for (float i = 0.0f; i < (float)x; i += 1.0f)
    {
        for (float j = 0.0f; j < (float)y; j += 1.0f)
        {
            for (float k = 0.0f; k < (float)z; k += 1.0f)
            {
                size_t voxel_index = (size_t)i + (size_t)j*x + (size_t)k*x*y;
                float distance = sqrt((cx - i)*(cx - i) + (cy - j)*(cy - j) + (cz - k)*(cz - k));
                if (distance < radius)
                {
                    mask[voxel_index] = 1.0f;
                }
                else if (sigma != 0)
                {
                    float gaussFallOff = exp(-((distance - radius) / sigma)*((distance - radius) / sigma));
                    if (gaussFallOff >= exp(-2))
                        mask[voxel_index] = gaussFallOff;
                }
            }
        }
    }

    return mask;
}

vector<float> MaskRoutines::createSphericalMask(size_t x, size_t y, size_t z)
{
    float radius = floor((min(min(x, y), z) - 1.0f) / 2.0f);

    return createSphericalMask(x, y, z, radius, 0);
}


vector<complex<double>> MaskRoutines::shiftVolume(vector<float>& volume, float x, float y, float z, float sx, float sy, float sz)
{
    float ox = -floor(x / 2.0f);
    float oy = -floor(y / 2.0f);
    float oz = -floor(z / 2.0f);

    sx /= x;
    sy /= y;
    sz /= z;

    vector<complex<double>> filter;
    filter.resize(x*y*z);

    complex<double> c(0, 1);
    for (float i = 0; i < x; i++)
    {
        for (float j = 0; j < y; j++)
        {
            for (float k = 0; k < z; k++)
            {
                double expValue = (i + ox)*sx + (j + oy)*sy + (k + oz)*sz;
                filter[i + j*x + k*x*y] = exp(-2.0*M_PI*c*expValue);
            }
        }
    }

    FFTRoutines::complex3DTransform(x, y, z, volume, filter);

	return	filter;
}

void MaskRoutines::symmetrizeVolume(vector<float>& volume, size_t x, size_t y, size_t z, int symmetry)
{
    if (symmetry == 1)
        return;

    vector<size_t> dim={x,y,z};
    Volume tempVolume(volume,dim);

    MaskRoutines::symmetrizeVolume(tempVolume,x,y,z,symmetry);

    copy(tempVolume.Data().begin(), tempVolume.Data().end(), volume.begin());

}

void MaskRoutines::symmetrizeVolume(Volume& volume, size_t x, size_t y, size_t z, int symmetry)
{
    if (symmetry == 1)
        return;

    //vector<float> rotVolume(volume.size());
    vector<float> tempVolume(volume.Data().size());
    copy(volume.Data().begin(), volume.Data().end(), tempVolume.begin());


    float phi_step = 360.0f / symmetry;
    for (int i = 1; i < symmetry; i++)
    {
        Quaternion q( phi_step*(i), 0, 0);
        volume.rotate(q);
        tempVolume = tempVolume + volume.Data(Volume::DataToProcess::rotatedD);
    }

    vec_div(tempVolume, (float)symmetry);

    volume.update(tempVolume);
}

void MaskRoutines::rotatePoint(float& x, float& y, float& z, float phi, float psi, float theta)
{
    double rPhi = phi / 180.0 * M_PI;
    double rPsi = psi / 180.0 * M_PI;
    double rTheta = theta / 180.0 * M_PI;

    float ox = x;
    float oy = y;
    float oz = z;

    float a = cos(rPsi);
    float b = sin(rPsi);
    float c = cos(rTheta);
    float d = sin(rTheta);
    float e = cos(rPhi);
    float f = sin(rPhi);

    x = ox*(a*e - b*c*f) - oy*(a*f + b*c*e) + oz * b*d;
    y = ox*(b*e + a*c*f) + oy*(a*c*e - b*f) - oz * a *d;
    z = ox*d*f + oy* d*e + oz*c;
}

vector<float> MaskRoutines::rotatePoint(vector<float>& point, vector<float>& angles)
{
    float px=point[0];
    float py=point[1];
    float pz=point[2];

    rotatePoint(px,py,pz,angles[0],angles[1],angles[2]);

    vector<float> p={px,py,pz};
    return p;

}
