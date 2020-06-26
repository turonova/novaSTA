#include "wedgeMaskGenerator.h"
#include "mpiRoutines.h"
#include "mrcStack.h"
#include "emFile.h"
#include <vector>
#include <algorithm>

using namespace vr;


WedgeMaskGenerator* WedgeMaskGenerator::create(ParameterSetup& iParams, int numberOfParts, int partNumber)
{
    switch (iParams.WedgeType())
    {
        case bowtie:
        {
            return new wmg_Bowtie(iParams, numberOfParts, partNumber);
        }
        case star:
        {
            return new wmg_Star(iParams, numberOfParts, partNumber);
        }
        case inputMasks:
        {
            return new wmg_InputMasks(iParams, numberOfParts, partNumber);
        }
        case singleMask:
        {
            return new wmg_SingleMask(iParams, numberOfParts, partNumber);
        }
    }
}

WedgeMaskGenerator::WedgeMaskGenerator(ParameterSetup& iParams)
{
    params = iParams;
    maskSize = iParams.VolumeSize();
    wedgeInputName = iParams.WedgeDataName();
    wedgeOutputName=iParams.WedgeMaskName();
}

void WedgeMaskGenerator::getTomoNumbers()
{
 	for (size_t i = 0; i < inputData.size(); i = i + entryOffset)
    {
        pair<map<unsigned int,size_t>::iterator,bool> ret=tomoNumbersPos.insert(pair<unsigned int, size_t>(inputData[i],i));
		
		if(ret.second==true)
		{
			tomoNumbers.push_back(inputData[i]);
		}
    }

}

void WedgeMaskGenerator::getTomoNumbersPos(map<unsigned int,size_t>& aTomoNumbersPos)
{
	aTomoNumbersPos=tomoNumbersPos;
}

void WedgeMaskGenerator::getTomoNumbers(vector<unsigned int>& aTomoNumbers)
{

	aTomoNumbers.resize(0);

	for (vector<unsigned int>::iterator itTomo = tomoListStart; itTomo != tomoListEnd; itTomo++)
    {
       	aTomoNumbers.push_back(*itTomo);
	}
}

string WedgeMaskGenerator::createNameWithZeroDigits(string baseName, unsigned int tomoNumber, string extension, int aTomoDigits)
{
    stringstream sNumber;
    sNumber << tomoNumber;
    string sTomoNumber = sNumber.str(); //Convert MotlData record number to string

    int numOfZero = aTomoDigits - sTomoNumber.length();
    if (numOfZero < 0)
    {
        cout << "Error (Negative value) while creating tomogram name. Digits in tomoDigits parameter should be greater or equal to digits in tomograms' name." << endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        sTomoNumber.insert(0, numOfZero, '0');                                  //Add zeros; depend on current TomoDigits
        string finalName = baseName + sTomoNumber + extension;       //Add exentsion on the end
        return finalName;
    }
}

wmg_Star::wmg_Star(ParameterSetup& iParams, int numberOfParts, int partNumber) : WedgeMaskGenerator(iParams)
{

}

void wmg_Star::generateMask()
{

}

wmg_SingleMask::wmg_SingleMask(ParameterSetup& iParams, int numberOfParts, int partNumber) : WedgeMaskGenerator(iParams)
{

}

void wmg_SingleMask::generateMask()
{

}

wmg_Bowtie::wmg_Bowtie(ParameterSetup& iParams, int numberOfParts, int partNumber) : WedgeMaskGenerator(iParams)
{

}

void wmg_Bowtie::generateMask()
{

}


wmg_InputMasks::wmg_InputMasks(ParameterSetup& iParams, int numberOfParts, int partNumber) : WedgeMaskGenerator(iParams)
{
    entryOffset = 4;    // tomo number + 3 coordinates

    tomoDigits = params.TomoDigits();

    emFile::read(wedgeInputName, inputData);

    mask.resize(maskSize*maskSize*maskSize);

    getTomoNumbers();

    MPIRoutines::computeEvenDistribution<unsigned int>(tomoNumbers, tomoListStart, tomoListEnd, numberOfParts, partNumber);
}


void wmg_InputMasks::generateMask()
{
    vector<float> subtomogram(maskSize*maskSize*maskSize);
    vector<float> subtomoFFT(maskSize*maskSize*maskSize);

    vector<size_t> subtomoSize = vector<size_t>{maskSize, maskSize, maskSize};

    for (vector<unsigned int>::iterator itTomo = tomoListStart; itTomo != tomoListEnd; itTomo++)
    {
        unsigned int tomoID = *itTomo;

		size_t tomoVectorIndex=tomoNumbersPos.find(tomoID)->second;
        vector<float>::iterator it = inputData.begin()+tomoVectorIndex;

        unsigned int currentTomo = tomoID;
		
		MRCStack*  tomogram = new MRCStack(createNameWithZeroDigits(params.TomoName(), currentTomo, ".rec",tomoDigits), false, false, true);
        vector<unsigned int> resolution = tomogram->getResolution();

        vector<float> maskShifted(maskSize*maskSize*maskSize,0.0f);
        float subtomoCount = 0;

        while (currentTomo == tomoID)
        {
            
            vector<float> centerCoordinates = vector<float>{*(it + 1), *(it + 2), *(it + 3) };

			vector<int> startingCoordinates;
            
            bool skipPosition= getStartingCoordinates(centerCoordinates, startingCoordinates,resolution);

            if (!skipPosition)
			{
                vector<size_t> coordinates = vector<size_t>{startingCoordinates[0], startingCoordinates[1], startingCoordinates[2]};
                tomogram->getSubvolume(&subtomogram[0], coordinates, subtomoSize);
                
                vr::vec_normalize(subtomogram);

                FFTRoutines::forwardComplex3DTransformPS(maskSize, maskSize, maskSize, subtomogram, subtomoFFT);
                maskShifted = maskShifted + subtomoFFT;
                subtomoCount++;
            }

			it = it + entryOffset;			
			
			if(it<inputData.end())		
				currentTomo = *it;
			else			
				break;		
        }

        vr::vec_div(maskShifted, subtomoCount);
        float norm_factor = (*std::max_element(maskShifted.begin(), maskShifted.end()));
        
        vr::vec_div(maskShifted, norm_factor);

        FFTRoutines::fftshift(maskShifted, mask, maskSize, maskSize, maskSize, -1);

        size_t center_dim = maskSize / 2;
        
        mask[center_dim + center_dim*maskSize + center_dim*maskSize * maskSize] = (mask[center_dim + center_dim*maskSize + center_dim*maskSize * maskSize]
            + mask[center_dim + (center_dim + 1)*maskSize + center_dim*maskSize * maskSize]
            + mask[(center_dim + 1) + (center_dim + 1)*maskSize + center_dim*maskSize * maskSize]
            + mask[(center_dim + 1) +  center_dim *maskSize + center_dim*maskSize * maskSize]
            + mask[(center_dim + 1) + (center_dim - 1)*maskSize + center_dim*maskSize * maskSize]
            + mask[center_dim + (center_dim - 1)*maskSize + center_dim*maskSize * maskSize]
            + mask[(center_dim - 1) + (center_dim - 1)*maskSize + center_dim*maskSize * maskSize]
            + mask[(center_dim + 1) + center_dim *maskSize + center_dim*maskSize * maskSize]
            + mask[(center_dim - 1) + (center_dim + 1)*maskSize + center_dim*maskSize * maskSize])/ 9;

        emFile::write(createNameWithZeroDigits(wedgeOutputName + "_", tomoID, ".em", tomoDigits), mask, maskSize, maskSize, maskSize);


        delete tomogram;
    }
}




bool wmg_InputMasks::getStartingCoordinates(const vector<float>& centerCoord, vector<int>& startCoord, vector<unsigned int>& resolution)
{   
    // setting starting coordinates x,y,z for subtomo
    for (unsigned int i = 0; i<centerCoord.size(); i++)
    {
        int startPos = (round(centerCoord[i] - (maskSize / 2))) - 1;
        if (startPos < 0 || (startPos + maskSize) >= resolution[i])
            return true;

        startCoord.push_back(startPos);
    }

    return false;

}
