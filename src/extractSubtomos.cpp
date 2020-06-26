#include "extractSubtomos.h"
#include "vectorRoutines.h"
#include "mrcStack.h"
#include "emFile.h"
#include "defines.h"
#include "float.h"
#include <string>
#include <vector>
#include <sstream>
#include <math.h>
#include <fstream>

using namespace std;
using namespace vr;

ExtractSubtomos::ExtractSubtomos(ParameterSetup& iParams, int iProcessID, int iNumberOfProcesses)
{
    params = iParams;
    processID = iProcessID;
    numberOfProcesses = iNumberOfProcesses;

    Motl::read(motlData,params.MotlName());
};

ExtractSubtomos::ExtractSubtomos(ParameterSetup& iParams, vector<MotlEntry>& iData, int iProcessID)
{
    processID = iProcessID;
    params = iParams;
    motlData = iData;
};

ExtractSubtomos::ExtractSubtomos(ParameterSetup& iParams, Log& iLog, vector<MotlEntry>& iData, int iProcessID)
{
    processID = iProcessID;
    params = iParams;
    motlData = iData;
    logger = iLog;
};


void ExtractSubtomos::getSubtomoStats(vector<float>& subtomogram, float subtomoNumber)
{
    float minValue = FLT_MAX;
    float maxValue = -FLT_MAX;

    float mean = 0.0f;

    for(size_t i=0; i<subtomogram.size(); i++)
    {
        if(maxValue < subtomogram[i])
            maxValue = subtomogram[i];

        if(minValue > subtomogram[i])
            minValue = subtomogram[i];

        mean += subtomogram[i];
    }

    mean = mean / (float) subtomogram.size();

    float variance=0.0f;

    for(size_t i=0;i<subtomogram.size();i++)
    {
        float diff = subtomogram[i] - mean;
        variance = variance + diff*diff;
    }

    variance = variance/ (float)subtomogram.size();

    float standardDeviation = sqrt(variance);

    statsFile << currentTomoNumber << "\t";
    statsFile << subtomoNumber << "\t";
    statsFile << mean << "\t";
    statsFile << minValue << "\t";
    statsFile << maxValue << "\t";
    statsFile << variance << "\t";
    statsFile << standardDeviation << "\n";

    currentSubtomoMean = mean;

}

void ExtractSubtomos::extractSubtomosPrep()
{
    subtomoSize = params.VolumeSize();
    currentTomoNumber = 0;
}


string ExtractSubtomos::createTomoName(float recordNumber,string extension)
{
    int tomoDigits = params.TomoDigits();
    stringstream motlDataStreamString;
    motlDataStreamString << recordNumber;
    string motlDataString = motlDataStreamString.str(); //Convert MotlData record number to string

    int numOfZero = tomoDigits - motlDataString.length();
    if(numOfZero < 0)
    {
        cout << "Error (Negative value) while creating tomogram name. Digits in tomoDigits parameter should be greater or equal to digits in tomograms' name." << endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        motlDataString.insert(0,numOfZero,'0');                                     //Add zeros; depend on current TomoDigits
        string tomoName = params.TomoName() + "/" + motlDataString + "." + extension;     //Add exentsion on the end
        return tomoName;
    }
}

string ExtractSubtomos::createSubtomoName(size_t subtomoNumber)
{
    string motlSubtomoName = params.SubtomoName();
    stringstream motlNumberStringStream;
    motlNumberStringStream << subtomoNumber;
    string motlNumberString = motlNumberStringStream.str();         //Converts tomoNumber to string

    string subtomoName;
    subtomoName = motlSubtomoName + "_" + motlNumberString + ".em"; //Create [subTomoName]_[motlNumber].em format in string
    return subtomoName;
}

vector<int> ExtractSubtomos::startPosition(vector<int>& coord )
{
    vector<int> startCoordinates;

    // setting starting coordinates x,y,z for subtomo
    for(unsigned int i=0; i<coord.size(); i++)
    {
	    int startPos = (round(coord[i] - (subtomoSize/2))) - 1;
	    startCoordinates.push_back(startPos);
    }
    return startCoordinates;
}

void ExtractSubtomos::checkBoundaries(vector<size_t>& realSize, vector<int>& realShift, vector<int>& startCoord, bool& cutSubtomo)
{
    cutSubtomo = false;

    for(unsigned int i = 0; i < startCoord.size(); i++)
    {
        //subtomogram has an axis outside the resolution in + value
        if(startCoord[i]+subtomoSize >= resolution[i])
        {
            realSize[i] = subtomoSize - ((startCoord[i]+subtomoSize) - resolution[i]);
            cutSubtomo = true;
        }
        else if(startCoord[i] < 0) //subtomogram has an axis outside resolution in - value
        {
            realSize[i] = subtomoSize + startCoord[i];
            realShift[i] = startCoord[i] * (-1);
            startCoord[i] = 0;
            cutSubtomo = true;
        }
    }
}

void ExtractSubtomos::normalizeSubtomo(vector<float> &subtomogram)
{
    float mean = vr::vec_mean(subtomogram);

    float variance = vr::vec_variance(subtomogram,mean);

    float stdv = sqrt(variance);

    for(unsigned int i=0;i<subtomogram.size();i++)
    {
        subtomogram[i] = (subtomogram[i] - mean) / stdv;
    }
}


/*void ExtractSubtomos::writeSubtomogramStats()
{
    ofstream writeSubtomoStats;
    for(unsigned int i=0;i<subtomoStats.size();i++)
    {
        float recordNumber = subtomoStats[i][0];
        writeSubtomoStats.open(createTomoName(recordNumber,"stats"),ofstream::out);

        // starting from 1 as we want to skip the tomo number
        for(unsigned int j=1;j<subtomoStats[i].size();j++)
        {
            writeSubtomoStats << subtomoStats[i][j];

            if((j%6)==0)
                writeSubtomoStats << endl;
            else
                writeSubtomoStats << " ";
        }

        writeSubtomoStats.close();
    }
}*/

bool ExtractSubtomos::createSubtomogram(size_t subtomoNumber,vector<int>& coord)
{
    for(unsigned int i=0; i< coord.size(); i++)
    {
        if( coord[i] >= resolution[i]  || coord[i] < 0 )
        {
            cout << "Subtomogram " << subtomoNumber << " is out of range and will be skipped."<< endl;
            return false;
        }
    }

    vector<size_t> realSize(3);
    vector<int> realShift(3);
    fill(realSize.begin(),realSize.end(),subtomoSize);
    fill(realShift.begin(),realShift.end(),0);

    bool subtomoCut;

    vector<int> startCoord = startPosition(coord);
    vector<float> subtomogram(subtomoSize*subtomoSize*subtomoSize);  //Preparing subtomogram vector

    checkBoundaries(realSize,realShift,startCoord,subtomoCut);

    vector<size_t> finalStartCoord(3);
    for(unsigned int i=0; i< startCoord.size();i++)
        finalStartCoord[i]=(size_t)startCoord[i];

    if(subtomoCut)
    {
        fill(subtomogram.begin(),subtomogram.end(),volume->getInputMeanValue());        //Filling vector with average values
        vector<float> cutSubtomogram(realSize[0]*realSize[1]*realSize[2]);
        volume->getSubvolume(&cutSubtomogram[0],finalStartCoord,realSize);

        //Copying part of tomogram to subtomogram
        int x,y,z;
        for(z=0;z<realSize[2];z++)
        {
            for(y=0;y<realSize[1];y++)
            {
                for(x=0;x<realSize[0];x++)
                {
                    size_t index1 = (size_t)((x+realShift[0])+(y+realShift[1])*subtomoSize+(z+realShift[2])*subtomoSize*subtomoSize);
                    size_t index2 = (size_t)( x + y*realSize[0] + z*realSize[0]*realSize[1] );
                    subtomogram[index1]=cutSubtomogram[index2];
                }
            }
        }
    }
    else
    {
        volume->getSubvolume(&subtomogram[0],finalStartCoord,realSize);
    }

    normalizeSubtomo(subtomogram);
    getSubtomoStats(subtomogram,subtomoNumber);

    //Writing subtomogram data to disc
    emFile::write(createSubtomoName(subtomoNumber),subtomogram,subtomoSize,subtomoSize,subtomoSize);

    return true;

}


/*void ExtractSubtomos::tomoExist(unsigned int tomoNumber)
{
    bool isAlreadyStored=false;

    for(unsigned int i=0;i<subtomoStats.size();i++)
    {
        if(subtomoStats[i][0] == tomoNumber)
        {
            isAlreadyStored = true;
            statsIndex = i;
        }
    }

    if(!isAlreadyStored)
    {
        vector<float> subtomos;
        subtomos.push_back( (float)tomoNumber );
        subtomoStats.push_back(subtomos);
        statsIndex = subtomoStats.size()-1;
    }
}*/

void ExtractSubtomos::extractSubtomograms()
{
    unsigned int motlTomoNumber;

    unsigned int motlElements = motlData.size();

    vector<int> coord(3);

    int start,end,step;

    if(params.Rerun()==-1)
    {
        start=motlElements-1;
        end=-1;
        step = -1;
    }
    else
    {
        start=0;
        end=motlElements;
        step=1;
    }

    //Creating subtomograms based on how many motlData there are
    for(int i=start;i!=end;i+=step)
    {
        if(subtomoExist(motlData[i].subtomoNumber))
        {
            motlData[i].subtomoMean = currentSubtomoMean;
            continue;
        }

        motlTomoNumber = (unsigned int) motlData[i].tomoNumber;

        if(motlTomoNumber != currentTomoNumber)
        {
            currentTomoNumber = motlTomoNumber;
            volume = new MRCStack(createTomoName(currentTomoNumber,"rec"),false, false, true);
            resolution = volume->getResolution();
            //tomoExist(currentTomoNumber);
        }

        coord[0] = motlData[i].x;
        coord[1] = motlData[i].y;
        coord[2] = motlData[i].z;

        bool subtomoCreated = createSubtomogram(motlData[i].subtomoNumber,coord);	//Creating subtomogram

        if(!subtomoCreated)
        {
            motlData[i].classNumber = -2.0;
        }
        else
        {
            motlData[i].subtomoMean = currentSubtomoMean;
        }
    }

    if(currentTomoNumber!=0)
        delete volume;
}

bool ExtractSubtomos::subtomoExist(unsigned int subtomoNumber)
{
    if(params.Rerun()!=-1)
        return false;

    try
    {
        emFile* subtomoFile = new emFile(createSubtomoName(subtomoNumber));
        currentSubtomoMean=vr::vec_mean(subtomoFile->Data());

        delete subtomoFile;

        return true;

    }
    catch(exception& e)
    {
        return false;
    }
}

string ExtractSubtomos::getSubtomoFolder(string path)
{
    string filename=path;

    // Remove filename if present
    // Do this before extension removal incase directory has a period character.
    const size_t last_slash_idx = filename.find_last_of("\\/");
    if (std::string::npos != last_slash_idx)
    {
        filename.erase(last_slash_idx + 1,filename.length());
    }

    return filename;
}

void ExtractSubtomos::prepareStatsFile()
{
    stringstream statsFilename;
    statsFilename << getSubtomoFolder(params.SubtomoName()) << "stats" << params.SplitType() << "_" << processID << ".txt";

    if(params.Rerun()==-1)
    {
        statsFile.open(statsFilename.str(), std::ofstream::out | std::ofstream::app);
        statsFile << "\n"; //necessary in the case where the last entry of the file did not contain whole line
    }
    else
        statsFile.open(statsFilename.str(), std::ofstream::out);
}

vector<MotlEntry> ExtractSubtomos::getMotlData()
{
    return motlData;
}

void ExtractSubtomos::run()
{
    prepareStatsFile();
    extractSubtomosPrep();
    extractSubtomograms();
    //writeSubtomogramStats();

    statsFile.close();
};
