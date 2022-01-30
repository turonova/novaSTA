#include "motl.h"
#include "defines.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include "vectorRoutines.h"

using namespace std;

template void Motl::computeEvenDistribution<MotlEntry>(vector<MotlEntry>& wholeData, vector<MotlEntry>::iterator& itStart, vector<MotlEntry>::iterator& itEnd, int numberOfParts, int partNumber);

Motl::Motl(string motlName)
{
    emFile* motlFile = new emFile(motlName);

    originalSize = motlFile->Data().size();

    convertData(motlFile->Data(),data);

    numberOfParticles = data.size();
    originalNumberOfParticles = data.size();

    numberOfFeatures = 1;
    numberOParticlesPerFeature.push_back(numberOfParticles);

    dataAsVector.resize(originalSize);

    copy(motlFile->Data().begin(),motlFile->Data().end(),dataAsVector.begin());

}


Motl::Motl(string motlName, int numberOfParts, int partNumber)
{
    emFile* motlFile = new emFile(motlName);

    originalSize = motlFile->Data().size();
    originalNumberOfParticles = originalSize / 20;

    vector<MotlEntry> allSubtomograms;
    convertData(motlFile->Data(),allSubtomograms);

    vector<MotlEntry>::iterator itStart, itEnd;

    computeEvenDistribution<MotlEntry>(allSubtomograms,itStart,itEnd,numberOfParts,partNumber);

    numberOfParticles = std::distance(itStart, itEnd);

    numberOfFeatures = 1;
    numberOParticlesPerFeature.push_back(numberOfParticles);


    data.resize(numberOfParticles);

    copy(itStart,itEnd,data.begin());
}

template <typename T>
void Motl::computeEvenDistribution(vector<T>& wholeData, typename vector<T>::iterator& itStart, typename vector<T>::iterator& itEnd, int numberOfParts, int partNumber)
{

    size_t particlesPerPart = (size_t)round((float)wholeData.size()/ ((float)numberOfParts));

    float particleRest = (float)wholeData.size() - (float)particlesPerPart*numberOfParts;
    float particlesUnchanged = numberOfParts-fabs(particleRest);

    int particleChange;

    if(particleRest<0)
        particleChange=-1;
    else
        particleChange=1;

    if(partNumber<particlesUnchanged)
    {
        itStart = wholeData.begin() + particlesPerPart*partNumber;
        itEnd = itStart + particlesPerPart;
    }
    else
    {
        itStart = wholeData.begin() + particlesPerPart*particlesUnchanged+(particlesPerPart+particleChange)*(partNumber-particlesUnchanged);

        if (partNumber == (numberOfParts - 1))
            itEnd = wholeData.end();
        else
            itEnd = itStart + particlesPerPart+particleChange;

    }
}

Motl::Motl(string motlName, unsigned int featureType)
{
    int numberOfParts = MPIRoutines::getCommSize();
    int partNumber = MPIRoutines::getProcessID();

    splitDataByFeature(motlName,featureType,numberOfParts,partNumber);
}

Motl::Motl(string motlName, unsigned int splitFeatureIndex, int numberOfParts, int partNumber)
{
    splitDataByFeature(motlName,splitFeatureIndex,numberOfParts,partNumber);
}

size_t Motl::getNumberOfAllParticles()
{
    size_t numberOfAllParticles = (size_t) MPIRoutines::reduceNumber(numberOfParticles);

    return numberOfAllParticles;
}

void Motl::writeFullData(string outputName,bool sortBySubtomoNumber)
{
    vector<float> motlData = getDataAsVector();
   // int motlSize = motlData.size();

    vector<float> fullMotlData;

    MPIRoutines::gatherData(motlData,fullMotlData);

    if(MPIRoutines::getProcessID()==0)
    {
        if(sortBySubtomoNumber)
        {
            vector<MotlEntry> fullMotl;
            convertData(fullMotlData,fullMotl);
            stable_sort(fullMotl.begin(), fullMotl.end(),MotlEntry::sortBySubtomoNumber);
            convertData(fullMotl,fullMotlData);
        }

        write(fullMotlData,outputName);
    }
}

float Motl::cleanByMeanGreyValue()
{

	float meanValue=0.0f;
    vector<float> partialMeans(data.size());

    for(size_t m=0; m<data.size(); m++)
    {
    	meanValue+=data[m].subtomoMean;
        partialMeans[m]=data[m].subtomoMean;
    }

	float finalMean=MPIRoutines::reduceFloatNumber(meanValue);
	size_t allParticles=MPIRoutines::reduceNumber(numberOfParticles);

    vector<float> combinedMeans;
	MPIRoutines::gatherData(partialMeans,combinedMeans);

	float sigma;
	
	if(MPIRoutines::getProcessID()==0)
	{
		finalMean=finalMean/(float)allParticles;
		sigma=sqrt(vr::vec_variance(combinedMeans,finalMean));
	}

	MPIRoutines::bcastFloatNumber(sigma);

    for(size_t m=0; m<data.size();m++)
    {
        if( (data[m].subtomoMean<(-sigma)) || (data[m].subtomoMean>sigma) )
            data[m].classNumber=-2.0f;
    }

    cleanClass(-2.0f);

	return sigma;
}

void Motl::splitDataByFeature(string motlName, unsigned int splitFeatureIndex, int numberOfParts, int partNumber)
{
    emFile* motlFile = new emFile(motlName);

    originalSize = motlFile->Data().size();
    originalNumberOfParticles = originalSize / 20;

    vector<MotlEntry> allSubtomograms;
    convertData(motlFile->Data(),allSubtomograms);

    stable_sort(allSubtomograms.begin(), allSubtomograms.end(),MotlEntry::sortByFeature);

    vector<vector<size_t>> featureBoundaries;
    size_t currentFeature=0;
    size_t firstElement;

    for (size_t m = 0; m < allSubtomograms.size(); m++)
    {
        if(currentFeature!=allSubtomograms[m].featureNumber)
        {
            //last element for given feature
            if(currentFeature!=0)
            {
                vector<size_t> currentEntry{currentFeature,firstElement,m};
                featureBoundaries.push_back(currentEntry);
            }

            // update current feature number
            currentFeature=(size_t)allSubtomograms[m].featureNumber;
            firstElement = m;
        }

        if(m==(allSubtomograms.size()-1))
        {
            vector<size_t> currentEntry{currentFeature,firstElement,m+1};
            featureBoundaries.push_back(currentEntry);
        }

    }

    numberOfParticles = 0;

    vector<vector<size_t>>::iterator itStart, itEnd;
    computeEvenDistribution<vector<size_t>>(featureBoundaries,itStart,itEnd,numberOfParts,partNumber);


    numberOfFeatures = std::distance(itStart, itEnd);
    numberOfParticles = 0;

    for(vector<vector<size_t>>::iterator it=itStart; it<itEnd; it++)
    {
        data.insert(data.begin()+numberOfParticles,allSubtomograms.begin()+(*it)[1],allSubtomograms.begin()+(*it)[2]);
        numberOfParticles+=((*it)[2]-(*it)[1]);
        numberOParticlesPerFeature.push_back((*it)[2]-(*it)[1]);
    }
}

void Motl::cleanByDistance(float distanceThreshold)
{
    size_t particleOffset = 0;

    for(size_t f=0; f<numberOfFeatures;f++)
    {
        for(size_t p=particleOffset;p<numberOParticlesPerFeature[f]+particleOffset;p++)
        {
            //if(data[p].classNumber!=-3)
            if(data[p].classNumber>0)
            {
                for(size_t c=p+1;c<numberOParticlesPerFeature[f]+particleOffset;c++)
                {
                    if(data[c].classNumber>0)
                    {
                        float distance = data[p].computeDistance(data[c]);
                        if(distance<distanceThreshold)
                        {
                            if(data[p].metricValue<=data[c].metricValue)
                            {
                                data[p].classNumber=-3;
                                break;
                            }
                            else
                            {
                                data[c].classNumber=-3;
                            }
                        }
                    }
                }
            }
        }

        particleOffset+=numberOParticlesPerFeature[f];
    }

    cleanClass(-3);
}

void Motl::unifyClassNumber(int newClassNumber)
{
    size_t particleOffset = 0;

   /* for (size_t f = 0; f<numberOfFeatures; f++)
    {
        for (size_t p = particleOffset; p<numberOParticlesPerFeature[f] + particleOffset; p++)
        {
            data[p].classNumber = newClassNumber;
        }
        particleOffset += numberOParticlesPerFeature[f];
    }*/

    for (size_t p = 0; p<data.size(); p++)
    {
        data[p].classNumber = newClassNumber;
    }
}

void Motl::computeAngularDistance(StructureGeometry& structureGeom)
{
    size_t particleOffset = 0;

    for (size_t f = 0; f<numberOfFeatures; f++)
    {
        vector<MotlEntry>::iterator it = data.begin() + particleOffset;
        
        structureGeom.initStructure(it, numberOParticlesPerFeature[f]);
        structureGeom.computeAngularDistance();
        structureGeom.clearStructure();
        
        particleOffset += numberOParticlesPerFeature[f];
    }

    cleanClass(-3);
}

void Motl::cleanClass(float classNumber)
{
    vector<MotlEntry> cleaned;

    cleaned.reserve(data.size());

    numberOfParticles=0;

    size_t particleOffset = 0;
    for(size_t f=0; f<numberOfFeatures;f++)
    {
        size_t updatedFeatureParticles=0;
        for(size_t p=particleOffset;p<numberOParticlesPerFeature[f]+particleOffset;p++)
        {
            if(data[p].classNumber!=classNumber)
            {
                cleaned.push_back(data[p]);
                numberOfParticles++;
                updatedFeatureParticles++;
            }
        }

        particleOffset+=numberOParticlesPerFeature[f];
        numberOParticlesPerFeature[f]=updatedFeatureParticles;
    }

    data.clear();
    data.resize(numberOfParticles);

    copy(cleaned.begin(),cleaned.end(),data.begin());
}

void Motl::convertData(vector<float>& inputData, vector<MotlEntry>& outputData)
{
    for (size_t m = 0; m < inputData.size(); m = m + 20)
    {
        vector<float>::iterator itStart = inputData.begin() + m;
        outputData.push_back(MotlEntry(itStart));
    }
}

void Motl::convertData(vector<MotlEntry>& inputData, vector<float>& outputData)
{
    if(outputData.size()==0)
        outputData.resize(inputData.size()*20);

    for (size_t m = 0; m < inputData.size(); m++)
    {
        vector<float>::iterator itStart = outputData.begin() + m*20;
        inputData[m].convertToVector(itStart);
    }
}

vector<float> Motl::getDataAsVector()
{
    vector<float> outputData;

    convertData(data,outputData);

    return outputData;
}

vector<float> Motl::getDataAsVector(vector<MotlEntry>& inputData)
{
    vector<float> outputData;

    convertData(inputData,outputData);

    return outputData;
}

void Motl::getDataAsVector(vector<float>& outputData)
{
    convertData(data,outputData);
}

void Motl::read(vector<MotlEntry>& outputData, string inputName)
{
    vector<float> tempData;
    emFile::read(inputName,tempData);

    convertData(tempData,outputData);
}

void Motl::read(string inputName)
{
    emFile::read(inputName,dataAsVector);
    convertData(dataAsVector,data);

    numberOfParticles=data.size();
}

void Motl::write(vector<float>& inputData, string outputName)
{
    emFile::write(outputName,inputData,20,inputData.size()/20,1);
}

void Motl::write(vector<MotlEntry>& inputData, string outputName)
{
    vector<float> outputData;
    convertData(inputData,outputData);

    emFile::write(outputName,outputData,20,inputData.size(),1);
}

void Motl::write(string outputName)
{
    vector<float> outputData;
    convertData(data,outputData);

    emFile::write(outputName,outputData,20,data.size(),1);
}

vector<MotlEntry> Motl::getData()
{
    return data;
}

vector<MotlEntry>& Motl::Data()
{
    return data;
}

size_t Motl::getOriginalSize()
{
    return originalSize;
}

size_t Motl::getNumberOfParticles()
{
    return numberOfParticles;
}

size_t Motl::getNumberOfFeatures()
{
    return numberOfFeatures;
}

size_t Motl::getNumberOfParticlesPerFeatures(size_t featureID)
{
    return numberOParticlesPerFeature[featureID];
}

size_t Motl::getOriginalNumberOfParticles()
{
    return originalNumberOfParticles;
}

void Motl::splitIntoEvenOdd(vector<float>& evenMotl, vector<float>& oddMotl)
{
    vector<MotlEntry> even;
    vector<MotlEntry> odd;

    splitIntoEvenOdd(even,odd);

    convertData(even,evenMotl);
    convertData(odd,oddMotl);

}

void Motl::setAngles(size_t midx, float phi, float psi, float theta)
{
    data[midx].phi=phi;
    data[midx].psi=psi;
    data[midx].theta=theta;
}

void Motl::setData(const vector<MotlEntry>& inputData)
{
    data=inputData;
}


void Motl::setShifts(size_t midx, float xs, float ys, float zs)
{
    data[midx].xShift=xs;
    data[midx].yShift=ys;
    data[midx].zShift=zs;
}

void Motl::setClass(size_t midx, float classNumber)
{
    data[midx].classNumber = classNumber;
}

void Motl::setMetric(size_t midx, float metricValue)
{
    data[midx].metricValue=metricValue;
}


void Motl::splitIntoEvenOdd(vector<MotlEntry>& evenMotl, vector<MotlEntry>& oddMotl)
{
    for (size_t m = 0; m < data.size(); m++)
    {
        int particleNumber = (int)data[m].subtomoNumber;
        if( (particleNumber % 2) == 0)
        {
            evenMotl.push_back(data[m]);
        }
        else
        {
            oddMotl.push_back(data[m]);
        }
    }
}

void Motl::unbin(float binningFactor)
{
    if(binningFactor==1.0f)
        return;

    for (size_t m = 0; m < data.size(); m++)
    {
        float x = (data[m].x + data[m].xShift)*binningFactor;
        float y = (data[m].y + data[m].yShift)*binningFactor;
        float z = (data[m].z + data[m].zShift)*binningFactor;

        data[m].x=round(x);
        data[m].y=round(y);
        data[m].z=round(z);

        data[m].xShift = x - round(x);
        data[m].yShift = y - round(y);
        data[m].zShift = z - round(z);
    }
}


size_t Motl::getActiveParticles(float classNumber)
{
    size_t activeParticles = 0;

    for (size_t m = 0; m < data.size(); m++)
    {
        if (data[m].classNumber == classNumber)
            activeParticles += 1;
    }

    return activeParticles;
}

void Motl::combineEvenOdd(vector<MotlEntry>& evenMotl, vector<MotlEntry>& oddMotl, vector<MotlEntry>& output)
{
   size_t combinedSize = evenMotl.size() + oddMotl.size();

   output.resize(combinedSize);

   size_t idxEven = 0;
   size_t idxOdd = 0;

   for (size_t m = 0; m < combinedSize; m++)
   {
       int even = (int)evenMotl[idxEven].subtomoNumber;
       int odd = (int)oddMotl[idxOdd].subtomoNumber;

       if(even<odd)
       {
           output[m]=evenMotl[idxEven];
           idxEven++;

          if(idxEven>=evenMotl.size())
          {
              m++;
              output[m]=oddMotl[idxOdd];
              break;
          }
       }
       else
       {
           output[m]=oddMotl[idxOdd];
           idxOdd++;

           if(idxOdd>=oddMotl.size())
           {
               m++;
               output[m]=evenMotl[idxEven];
               break;
           }
       }
   }
}

void Motl::renumberParticles(bool renumber)
{
    if(!renumber)
        return;

    for (size_t m = 0; m < data.size(); m++)
    {
        data[m].subtomoNumber=m+1;
    }
}
