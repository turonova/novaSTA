#include "wedgeMask.h"
#include "maskRoutines.h"
#include "wedgeMaskGenerator.h"
#include "emFile.h"
#include <vector>
#include <algorithm>

WedgeMask::WedgeMask(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC) : Volume(aDim)
{
    size = aSize;
    useRosemanCC = aUseRosemanCC;
    wedgeBaseName = aWedgeName;

    shifted.resize(size);
}

WedgeMask* WedgeMask::create(ParameterSetup& iParams, size_t aSize, vector<size_t>& aDim)
{
    switch (iParams.WedgeType())
    {
        case bowtie:
        {
            return new wm_Bowtie(aSize, aDim, iParams.WedgeMaskName(), iParams.UseRosemanCC());
        }
        case star:
        {
            return new wm_Star(aSize, aDim, iParams.WedgeMaskName(), iParams.UseRosemanCC());
        }
        case inputMasks:
        {
            return new wm_InputMasks(aSize, aDim, iParams.WedgeMaskName(), iParams.UseRosemanCC(), (int) iParams.TomoDigits());
        }
        case singleMask:
        {
            return new wm_SingleMask(aSize, aDim, iParams.WedgeMaskName(), iParams.UseRosemanCC());
        }
    }
}

void WedgeMask::shiftAndReduce()
{
    size_t fftz = floor(dim[2] / 2.0f) + 1;
    shiftedReduced.resize(dim[0] * dim[1] * fftz);

    FFTRoutines::fftshift(data, shifted, dim[0], dim[1], dim[2], -1, 0);
  
    for (size_t x = 0; x < fftz; x++)
        for (size_t y = 0; y < dim[1]; y++)
            for (size_t z = 0; z < dim[2]; z++)
            {
                shiftedReduced[x + y*fftz + z*dim[1] * fftz] = shifted[x + y*dim[0] + z*dim[0] * dim[1]];
            }
}

void WedgeMask::loadMask(string maskName)
{
    emFile::read(maskName, data);
}

void WedgeMask::binarizeRotated()
{
    for (size_t i = 0; i < size; i++)
    {
        if (rotated[i] >= 0.5f)
            rotated[i] = 1.0f;
        else
            rotated[i] = 0.0f;
    }
}

vector<float>& WedgeMask::rotateWedgeMask(Quaternion q)
{
    rotate(q);
    return rotated;
}

void WedgeMask::filterWithMask(vector<float>& edgeMask)
{
    vr::vec_mult(shifted, edgeMask);

    size_t fftz = floor(dim[2] / 2.0f) + 1;
    shiftedReduced.resize(dim[0] * dim[1] * fftz);

    for (size_t i = 0; i < fftz; i++)
    {
        for (size_t j = 0; j < dim[1]; j++)
        {
            for (size_t k = 0; k < dim[2]; k++)
            {
                shiftedReduced[i + j*fftz + k*fftz*dim[1]] = shifted[i + j*dim[0] + k*dim[0] * dim[1]];
            }
        }
    }
}



vector<float>& WedgeMask::shiftedMask(bool reducedSize)
{
    if (!reducedSize)
        return shifted;
    else
        return shiftedReduced;
}

void WedgeMask::filter(vector<float> filterMask)
{
    for (size_t i = 0; i < dim[0] * dim[1] * fftz; i++)
    {
        shiftedReduced[i] = shiftedReduced[i] * filterMask[i];
    }
}



wm_InputMasks::wm_InputMasks(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC, int aTomoDigits) : WedgeMask(aSize, aDim, aWedgeName, aUseRosemanCC)
{
    tomoDigits = aTomoDigits;
}


bool wm_InputMasks::updateMask(unsigned int& currentWedge, unsigned int newWedge, bool maskEdges, bool shiftMask, const vector<float>& wedgeFilter)
{
    if (currentWedge == newWedge)
        return false;

    currentWedge = newWedge;
    loadMask(WedgeMaskGenerator::createNameWithZeroDigits(wedgeBaseName + "_", currentWedge, ".em", tomoDigits));

    if (shiftMask)
        shiftAndReduce();

    if (useRosemanCC && wedgeFilter.size() != 0)
        filter(wedgeFilter);

    return true;
}





wm_SingleMask::wm_SingleMask(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC) : WedgeMask(aSize, aDim, aWedgeName, aUseRosemanCC)
{
    filtered = false;
    loadMask(aWedgeName);
    shiftAndReduce();
}


bool wm_SingleMask::updateMask(unsigned int& currentWedge, unsigned int newWedge, bool maskEdges, bool shiftMask, const vector<float>& wedgeFilter)
{
    currentWedge = newWedge;

    if (!filtered && useRosemanCC && wedgeFilter.size() != 0)
    {
        filter(wedgeFilter);
        filtered = true;
    }

    return false;
}

wm_Star::wm_Star(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC) : WedgeMask(aSize, aDim, aWedgeName, aUseRosemanCC)
{
}



bool wm_Star::updateMask(unsigned int& currentWedge, unsigned int newWedge, bool maskEdges, bool shiftMask, const vector<float>& wedgeFilter)
{
    if (currentWedge == newWedge)
        return false;

    currentWedge = newWedge;
    generateMask(currentWedge, maskEdges, shiftMask);

    if (useRosemanCC && wedgeFilter.size() != 0)
        filter(wedgeFilter);

    return true;

}

void wm_Star::generateMask(unsigned int id, bool maskEdges, bool shift)
{
    //TODO
}




wm_Bowtie::wm_Bowtie(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC) : WedgeMask(aSize, aDim, aWedgeName, aUseRosemanCC)
{
    initWedgeList(aWedgeName);
}


vector<float>& wm_Bowtie::rotateWedgeMask(Quaternion q)
{
    rotate(q);
    binarizeRotated();
    return rotated;
}

void wm_Bowtie::initWedgeList(string filename)
{
    emFile* wedgeListFile = new emFile(filename);

    for (size_t i = 0; i < wedgeListFile->Data().size(); i = i + wedgeListFile->Header().x)
    {
        wedgeList.insert(pair<unsigned int, pair<float, float>>(wedgeListFile->Data()[i], pair<float, float>(wedgeListFile->Data()[i + 1], wedgeListFile->Data()[i + 2])));
    }

    delete wedgeListFile;
}

bool wm_Bowtie::updateMask(unsigned int& currentWedge, unsigned int newWedge, bool maskEdges, bool shiftMask, const vector<float>& wedgeFilter)
{
    if (currentWedge == newWedge)
        return false;

    currentWedge = newWedge;
    generateMask(currentWedge, maskEdges, shiftMask);

    if (useRosemanCC && wedgeFilter.size()!=0)
        filter(wedgeFilter);

    return true;
}

void wm_Bowtie::generateMask(unsigned int id, bool maskEdges, bool shift)
{
    pair<float, float> angles = wedgeList.find(id)->second;

    double minAngle=angles.first;
    double maxAngle=angles.second;

    double minAngleRad = minAngle*M_PI / 180.0;
    double maxAngleRad = maxAngle*M_PI / 180.0;

    double tan_min = tan((-M_PI / 2.0) - minAngleRad);
    double tan_max = tan(M_PI / 2.0 - maxAngleRad);

    float radius;

    if(maskEdges)
        radius = floor((min(min(dim[0], dim[1]), dim[2]) - 1.0f) / 2.0f);
    else
        radius = 2*dim[0]; // no masking

    fill(data.begin(), data.end(), 1.0f);

    // center for the shuffle
    size_t ms = (size_t)floor((float)dim[0] / 2.0f);
    vector<size_t> idx;
    size_t fftz = floor(dim[2] / 2.0f) + 1;

    if(shift)
    {
        //fill(shifted.begin(), shifted.begin()+dim[0]*dim[1], 0.0f);
        fill(shifted.begin(), shifted.end(), 1.0f);

        idx.resize(dim[0]);
        for (size_t i = 0; i < dim[0]; i++)
            idx[i] = i;

        std::rotate(idx.begin(), idx.begin() + ms, idx.end());

        shiftedReduced.resize(dim[0] * dim[1] * fftz);
        fill(shiftedReduced.begin(), shiftedReduced.end(), 1.0f);
    }

    float halfDim = floor(dim[0] / 2);
    for (float z = -halfDim; z < -halfDim + dim[2]; z++)
    {
        if (z == 0)
        {
            if(shift)// edge mask
            {
                for (float x = -halfDim; x < -halfDim + dim[0]; x++)
                for (int y = -halfDim; y < -halfDim +dim[1]; y++)
                {
                    float distance = sqrt(x*x + y*y + z*z);
                    if (distance >= radius)
                    {
                        size_t voxel_index = (size_t)(idx[(x + halfDim)] + idx[(y + halfDim)]*dim[0] + idx[(z + halfDim)]*dim[0]*dim[1]);
                        shifted[voxel_index] = 0.0f;

                        if(idx[(x + halfDim)]<fftz)
                        {
                            voxel_index = (size_t)(idx[(x + halfDim)] + idx[(y + halfDim)]*fftz + idx[(z + halfDim)]*fftz*dim[1]);
                            shiftedReduced[voxel_index] = 0.0f;
                        }
                    }
                }
            }
            continue;
        }

        for (float x = -halfDim; x < -halfDim + dim[0]; x++)
        {
            if ((tan_max > x / z) && (tan_min < x / z))
            {
                for (int y = 0; y < dim[1]; y++)
                {
                    data[(size_t)((x + halfDim) + y*dim[0] + (z + halfDim)*dim[0]*dim[1])] = 0.0f;

                    if(shift)
                    {
                        shifted[(size_t)(idx[(x + halfDim)] + idx[y]*dim[0] + idx[(z + halfDim)]*dim[0]*dim[1])] = 0.0f;

                        if(idx[(x + halfDim)]<fftz)
                        {
                            size_t voxel_index = (size_t)(idx[(x + halfDim)] + idx[y]*fftz + idx[(z + halfDim)]*fftz*dim[1]);
                            shiftedReduced[voxel_index] = 0.0f;
                        }
                    }
                }
            }
            else if(shift)// edge mask
            {
                for (int y = -halfDim; y < -halfDim +dim[1]; y++)
                {
                    float distance = sqrt(x*x + y*y + z*z);
                    if (distance >= radius)
                    {
                        size_t voxel_index = (size_t)(idx[(x + halfDim)] + idx[(y + halfDim)]*dim[0] + idx[(z + halfDim)]*dim[0]*dim[1]);
                        shifted[voxel_index] = 0.0f;

                        if(idx[(x + halfDim)]<fftz)
                        {
                            voxel_index = (size_t)(idx[(x + halfDim)] + idx[(y + halfDim)]*fftz + idx[(z + halfDim)]*fftz*dim[1]);
                            shiftedReduced[voxel_index] = 0.0f;
                        }
                    }
                }
            }
        }
    }
}

void wm_Bowtie::createMask(unsigned int id, bool shift)
{
    pair<float, float> angles = wedgeList.find(id)->second;
    MaskRoutines::createWedgeMask(data, dim[0], angles.first, angles.second);

    if (shift)
        FFTRoutines::fftshift(data, shifted, dim[0], dim[1], dim[2], -1, 0);
}