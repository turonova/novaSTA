#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "parameterSetup.h"
#include <algorithm>

ParameterSetup::ParameterSetup()
{

}

void ParameterSetup::applyEvenOddSplit(string nameUpdate)
{
    motlName.value = motlBaseName + "_" + nameUpdate;
    motlBaseName = motlBaseName + "_" + nameUpdate;


    if(numberOfIterations.value!=0 || createReference)
    {
        refOriginalName = refBaseName;
        refName.value = refBaseName + "_" + nameUpdate;
        refBaseName = refBaseName + "_" + nameUpdate;
        wedgeBaseName= refBaseName + "_wedgeSum";
    }

    splitType = nameUpdate;
}

string ParameterSetup::getFilename(string path)
{
    string filename=path;

    // Remove directory if present
    // Do this before extension removal incase directory has a period character.
    const size_t last_slash_idx = filename.find_last_of("\\/");
    if (std::string::npos != last_slash_idx)
    {
        filename.erase(0, last_slash_idx + 1);
    }

    // Remove extension if present.
    const size_t period_idx = filename.rfind('.');
    if (std::string::npos != period_idx)
    {
        filename.erase(period_idx);
    }

    return filename;
}

string ParameterSetup::createBaseName(string folder, string filename)
{
    stringstream base;
    base << folder << "/" <<filename;
    return base.str();
}


void ParameterSetup::checkParamSetup()
{
    if(!numberOfIterations.set)
        parameterNotSetError(numberOfIterations.name);

    if(numberOfIterations.value==0)
        return;

    if(coneAngle.set && coneSampling.set && inplaneAngle.set && inplaneSampling.set)
    {
        if(coneAngle.value.size()==coneSampling.value.size()
        && coneAngle.value.size()==inplaneAngle.value.size()
        && coneAngle.value.size()==inplaneSampling.value.size())
        {
            if(coneAngle.value.size()!=1 && coneAngle.value.size()!=numberOfIterations.value)
            {
                wrongParameterCombination();
            }
            else if (coneAngle.value.size()==1)
            {
                changingRotations=false;

                for(unsigned int i=1; i<numberOfIterations.value;i++)
                {
                    coneAngle.value.push_back(coneAngle.value[0]);
                    coneSampling.value.push_back(coneSampling.value[0]);
                    inplaneAngle.value.push_back(inplaneAngle.value[0]);
                    inplaneSampling.value.push_back(inplaneSampling.value[0]);
                }
            }
            else
            {
                changingRotations=true;
            }
        }
        else
        {
            wrongParameterCombination();
        }
    }
    else
        wrongParameterCombination();

	
    if(lowPass.set && highPass.set && highPass.value.size()==lowPass.value.size())
    {
        if(lowPass.value.size()!=1 && lowPass.value.size()!=numberOfIterations.value)
        {
            wrongParameterCombination();
        }
        else if (lowPass.value.size()==1)
        {
            changingBandPass=false;

            for(unsigned int i=1; i<numberOfIterations.value;i++)
            {
                lowPass.value.push_back(lowPass.value[0]);
                highPass.value.push_back(highPass.value[0]);
            }
        }
        else
        {
            changingBandPass=true;
        }
		
		for(unsigned int i = 0; i< lowPass.value.size(); i++)
		{
			if (lowPass.value[i] <= 0)
			{
				changingBandPass=true;
				if(!splitIntoEvenOdd || !createReference || !computeFSC)
				{
					cout << endl << "ERROR: Addaptive filter can be only used with splitIntoEvenOdd set to 1 and createReference set to 1 !!!" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
    }
    else
        wrongParameterCombination();
}

void ParameterSetup::setInternalParameters()
{
 
    if(rerun>0)
        realStartIndex=rerun;
    else
        realStartIndex=StartingIndex();

    if(rerun<=0)
        rotationAndFilterIndex=0;
    else
        rotationAndFilterIndex=rerun-StartingIndex();


    if(numberOfIterations.value!=0 || createReference)
    {
        refBaseName = createBaseName(FolderName(),getFilename(RefName()));
        motlBaseName = createBaseName(FolderName(),getFilename(MotlName()));
        wedgeBaseName= refBaseName + "_wedgeSum";
    }

    if(copySubtomograms)
        copySubtomoName =TempSubtomoFolder() + "/" + getFilename(SubtomoName());

    motlOriginalName = motlBaseName;
    stringstream logFile;

    unsigned int endingIndex = startingIndex.value+numberOfIterations.value;
    logFile << folderName.value << "/" << "log_ref" << startingIndex.value << "_to_ref" << endingIndex << ".txt";

    logFileName = logFile.str();

}

void ParameterSetup::wrongParameterCombination()
{
    cout << endl << "ERROR: The parameter combination is wrong!!!" << endl;
    cout << "Check the cone and in-plane parameters as well as high- and lowPass ones!" << endl;
    cout << "There has to be either one value per parameter or the number of values has to be the same as number of iterations!" << endl;
    exit(EXIT_FAILURE);
}

void ParameterSetup::parameterNotSetError(string paramName)
{
    cout << endl << "ERROR: The parameter " << paramName << " was not set!!!" << endl;
    exit(EXIT_FAILURE);
}

void separateNameAndValue(std::string& aLine, std::string& aParamName, std::string& aParamValue)
{
    for (unsigned int i = 0; i < aLine.length(); i++)
    {
        if (isspace(aLine[i]))
        {
            //do not want the space neither in name or value
            aParamName = aLine.substr(0, i);
            aParamValue = aLine.substr(i + 1, aLine.length());
            return;
        }
    }
    aParamName = aLine;
    aParamValue.clear();
}

std::string removeEmptySpacesFromString(string input)
{
    string output = input;
    output.erase(std::remove(output.begin(), output.end(), ' '), output.end());
    return output;
}

void checkEmptyString(inputParameter<string>& param)
{
    if (param.value.empty())
    {
        cout << endl << "ERROR: Following parameter was set to be empty: " << param.name << endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        param.set = true;
    }
}

void parseDimensions2D(std::string& aParamValue, int& aX, int& aY, const char separator)
{
    for (unsigned int i = 0; i < aParamValue.length(); i++)
    {
        if (aParamValue[i] == separator)
        {
            //do not want the space neither in name or value
            string tmpX = aParamValue.substr(0, i);
            string tmpY = aParamValue.substr(i + 1, aParamValue.length());
            aX = atoi(tmpX.c_str());
            aY = atoi(tmpY.c_str());
            return;
        }
    }
}

void parseFloats(std::string& aParamValue, float& aX, float& aY, const char separator)
{
    for (unsigned int i = 0; i < aParamValue.length(); i++)
    {
        if (aParamValue[i] == separator)
        {
            //do not want the space neither in name or value
            string tmpX = aParamValue.substr(0, i);
            string tmpY = aParamValue.substr(i + 1, aParamValue.length());
            aX = (float)atof(tmpX.c_str());
            aY = (float)atof(tmpY.c_str());
            return;
        }
    }
}

void parseFloats(std::string& aParamValue, vector<float>& values, const char separator)
{
    int lastStop=-1;

    for (unsigned int i = 0; i < aParamValue.length(); i++)
    {
        if (aParamValue[i] == separator)
        {
            //do not want the space neither in name or value
            string tmpX = aParamValue.substr(lastStop+1, i);
            lastStop = i;
            values.push_back((float)atof(tmpX.c_str()));

            while (aParamValue[i] == separator && i<aParamValue.length())
                i++;
        }
    }

    string tmpY = aParamValue.substr(lastStop + 1, aParamValue.length());
    string lastValue=removeEmptySpacesFromString(tmpY);

    if(!lastValue.empty())
        values.push_back((float)atof(tmpY.c_str()));
}

void ParameterSetup::parseComFile(const string tiltcom)
{

    std::string line;
    ifstream ifs(tiltcom.c_str());

    if (!ifs)
    {
        cout << "Could not open following file: " << tiltcom.c_str() << endl;
        return;
    }

    while (getline(ifs, line))
    {
        if (line[0] != '#' && line[0] != '$')
        {
            //cout << "Processing following line:" << line << endl;
            string paramName;
            string paramValue;
            separateNameAndValue(line, paramName, paramValue);
            storeValues(paramName, paramValue, ' ');
        }
    }
}

bool ParameterSetup::parameterWithoutValue(string paramName)
{
    bool paramWithoutValue = false;

    if (paramName == "-PERPENDICULAR" || paramName == "-AdjustOrigin")
    {
        paramWithoutValue = true;
    }

    return paramWithoutValue;
}

void ParameterSetup::parseCommandLine(std::vector<string> argList)
{
    for (std::vector<string>::iterator it = argList.begin(); it != argList.end(); it++)
    {
        string paramName = *it;
        string paramValue;
        std::vector<string>::iterator nextArgument = it + 1;
        if (nextArgument != argList.end() && !parameterWithoutValue(paramName))
        {
            paramValue = *nextArgument;
            it++;
        }
        else if (!parameterWithoutValue(paramName))
        {
            cout << "Missing value of the following parameter: " << *it << endl;
            return;
        }

        //if(paramName!="tilt.com")
        {
            paramName.erase(0, 1); //erase "-" at the beginning of command line arguments
            storeValues(paramName, paramValue, ',');
        }
    }
}

void ParameterSetup::printParameter(string paramName, string paramValue)
{
    stringstream additionaLine;

    additionaLine << paramName << " " << paramValue << endl;

    allParameters += additionaLine.str();
}


void ParameterSetup::storeValues(string paramName, string paramValue, char separator)
{

    printParameter(paramName,paramValue);

    if (paramName == "useGPU")
    {
        useGPU = atoi(paramValue.c_str());
    }
    else if (paramName == "createRef")
    {
        createReference = atoi(paramValue.c_str());
    }
    else if (paramName == "extractSubtomos")
    {
        extractSubtomos = atoi(paramValue.c_str());
    }
    else if (paramName == "rerun")
    {
        rerun = atoi(paramValue.c_str());
    }
    else if (paramName == "symmetry")
    {
        symmetry = (unsigned int) atoi(paramValue.c_str());
    }
    else if (paramName == "class")
    {
        subtomoClass = atoi(paramValue.c_str());
    }
    else if (paramName == "threshold")
    {
        threshold = atof(paramValue.c_str());
    }
    else if (paramName == "lowPassSigma")
    {
        lowPassSigma = atof(paramValue.c_str());
    }
    else if (paramName == "highPassSigma")
    {
        highPassSigma = atof(paramValue.c_str());
    }
    else if (paramName=="interpolationType")
    {
        interpolationType = atoi(paramValue.c_str());
    }
    else if (paramName=="interpolationSigma")
    {
        interpolationSigma = atof(paramValue.c_str());
    }
    else if (paramName == "motlBinFactor")
    {
        motlBinFactor = atof(paramValue.c_str());
    }
    else if (paramName == tomoDigits.name)
    {
        tomoDigits.value = atoi(paramValue.c_str());
        tomoDigits.set = true;
    }
    else if (paramName == motlName.name)
    {
        motlName.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(motlName);
    }
    else if (paramName == fscMask.name)
    {
        fscMask.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(fscMask);
    }
    else if (paramName == folderName.name)
    {
        folderName.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(folderName);
    }
    else if (paramName == "splitIntoEvenOdd")
    {
        splitIntoEvenOdd = atoi(paramValue.c_str());
    }
    else if (paramName == refName.name)
    {
        refName.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(refName);
    }
    else if (paramName == maskName.name)
    {
        maskName.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(maskName);
    }
    else if (paramName == ccMaskName.name)
    {
        ccMaskName.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(ccMaskName);
    }
    else if (paramName == wedgeMaskName.name)
    {
        wedgeMaskName.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(wedgeMaskName);
    }
    else if (paramName == wedgeDataName.name)
    {
        wedgeDataName.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(wedgeDataName);
    }
    else if (paramName == subtomoName.name)
    {
        subtomoName.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(subtomoName);
    }
    else if (paramName == tomoName.name)
    {
        tomoName.value = removeEmptySpacesFromString(paramValue.c_str());
        checkEmptyString(tomoName);
    }
    else if (paramName == numberOfIterations.name)
    {
        numberOfIterations.value = atoi(paramValue.c_str());
        numberOfIterations.set = true;
    }
    else if (paramName == volumeSize.name)
    {
        volumeSize.value = (size_t) atoi(paramValue.c_str());
        volumeSize.set = true;
    }
    else if (paramName == distanceThreshold.name)
    {
        distanceThreshold.value = (size_t) atof(paramValue.c_str());
        distanceThreshold.set = true;
    }
    else if (paramName == "cleanByDistance")
    {
        cleanByDistance = atoi(paramValue.c_str());
    }
    else if (paramName == "renumberParticles")
    {
        renumberParticles = atoi(paramValue.c_str());
    }
    else if  (paramName == "cleanByMeanGreyValue")
    {
        cleanByMeanGreyValue = atoi(paramValue.c_str());
    }
    else if (paramName == "cleanOutOfBoundsParticles")
    {
        cleanOutOfBoundsParticles=atoi(paramValue.c_str());
    }
    else if (paramName == "wedgeType")
    {
        wedgeType = atoi(paramValue.c_str());
    }
    else if (paramName == "createWedgeMasks")
    {
        createWedgeMasks = atoi(paramValue.c_str());
    }
    else if (paramName == lowPass.name)
    {
        parseFloats(paramValue,lowPass.value,separator);
        lowPass.set = true;
    }
    else if (paramName == highPass.name)
    {
        parseFloats(paramValue,highPass.value,separator);
        highPass.set = true;
    }
    else if (paramName == "useRosemanCC")
    {
        useRosemanCC = atoi(paramValue.c_str());
    }
    else if (paramName == "runSA")
    {
        runSA = atoi(paramValue.c_str());
    }
    else if (paramName == "phiIter")
    {
        cout << endl << "ERROR: Parameter phiIter is no longer used!!! " << endl;
        cout << "Use inplaneAngle instead!!!" << endl;
        exit(EXIT_FAILURE);
    }
    else if (paramName == "phiIncr")
    {
        cout << endl << "ERROR: Parameter phiIncr is no longer used!!! " << endl;
        cout << "Use inplaneSampling instead!!!" << endl;
        exit(EXIT_FAILURE);
    }
    else if (paramName == "angIter")
    {
        cout << endl << "ERROR: Parameter angIter is no longer used!!! " << endl;
        cout << "Use coneAngle instead!!!" << endl;
        exit(EXIT_FAILURE);
    }
    else if (paramName == "angIncr")
    {
        cout << endl << "ERROR: Parameter angIncr is no longer used!!! " << endl;
        cout << "Use coneSampling instead!!!" << endl;
        exit(EXIT_FAILURE);
    }
    else if (paramName == coneAngle.name)
    {
        parseFloats(paramValue,coneAngle.value,separator);
        coneAngle.set = true;
    }
    else if (paramName == coneSampling.name)
    {
        parseFloats(paramValue,coneSampling.value,separator);
        coneSampling.set = true;
    }
    else if (paramName == inplaneAngle.name)
    {
        parseFloats(paramValue,inplaneAngle.value,separator);
        inplaneAngle.set = true;
    }
    else if (paramName == inplaneSampling.name)
    {
        parseFloats(paramValue,inplaneSampling.value,separator);
        inplaneSampling.set = true;
    }
    else if (paramName == startingIndex.name)
    {
        startingIndex.value = atoi(paramValue.c_str());
        startingIndex.set = true;
    }
    else if (paramName == pixelSize.name)
    {
        pixelSize.value = atof(paramValue.c_str());
        pixelSize.set = true;
    }
    else if ( paramName == "copySubtomos")
    {
        copySubtomograms = atoi(paramValue.c_str());
    }
	else if ( paramName == "computeFSC")
    {
        computeFSC = atoi(paramValue.c_str());
    }
    else if ( paramName == "tempSubtomoFolder")
    {
        tempSubtomoFolder.value = removeEmptySpacesFromString(paramValue.c_str());
        tempSubtomoFolder.set = true;
    }
   /* else
    {
        if (paramName != "")
            cout << "Ignoring following parameter: " << paramName << endl;
    }*/
}

void ParameterSetup::initVariables()
{
    useGPU = false;
    createReference = false;
    extractSubtomos = false;

    symmetry = 1;
    subtomoClass = 1;
    
    threshold = 0.0f;
    
    lowPassSigma = 3.0f;
    highPassSigma = 2.0f;

    rerun = 0;
    useRosemanCC = false;
    copySubtomograms = false;

    splitIntoEvenOdd = false;
    splitType = "";

    motlBinFactor = 1.0f;

    cleanByDistance = false;
    cleanOutOfBoundsParticles = true;

    renumberParticles = false;
    cleanByMeanGreyValue = false;

    wedgeType = 0;
    createWedgeMasks = false;

    runSA = true;
    logFileName = "log.txt";

    interpolationType = 0;
    interpolationSigma = 0.3f;

	computeFSC = true;
}

ParameterSetup::ParameterSetup(std::vector<string> argList)
{
    //stringstream allParameters;

    initVariables();

    //check if tilt.com is among the parameters and if so parse it first - the command line parameters have more priority and can later on change any parameters set by tilt.com
    for (std::vector<string>::iterator it = argList.begin(); it != argList.end(); it++)
    {
        if ((*it) == "-param")
        {
            std::vector<string>::iterator nextArgument = it + 1;
            cout << "Parsing the following parameter file: " << *nextArgument << endl;
            parseComFile((*nextArgument));
        }
    }

    parseCommandLine(argList);

    if (runSA)
    {
        checkParamSetup();
        setInternalParameters();
    }
}

ParameterSetup::~ParameterSetup()
{}

void ParameterSetup::printParamFile()
{
    if(!folderName.set)
        return;

    ofstream  paramFile;
    
    string paramFileName;

    if (runSA)
    {
        stringstream paramFileNameS;
        unsigned int endingIndex = startingIndex.value+numberOfIterations.value;
        paramFileNameS << folderName.value << "/" << "params_ref" << startingIndex.value << "_to_ref" << endingIndex << ".txt";
        paramFileName = paramFileNameS.str().c_str();
    }
    else
    {
        paramFileName = "params_generated.txt";
    }

    paramFile.open(paramFileName);

    paramFile << allParameters;

    paramFile.close();
}

int ParameterSetup::Rerun()
{
    return rerun;
}

bool ParameterSetup::ExtractSubtomos()
{
    return extractSubtomos;
}

bool ParameterSetup::UseGPU()
{
    return useGPU;
}

bool ParameterSetup::CreateReference()
{
    return createReference;
}

unsigned int ParameterSetup::Symmetry()
{
    return symmetry;
}

int ParameterSetup::SubtomoClass()
{
    return subtomoClass;
}

float ParameterSetup::Threshold()
{
    return threshold;
}
float ParameterSetup::LowPassSigma()
{
    return lowPassSigma;
}

float ParameterSetup::HighPassSigma()
{
    return highPassSigma;
}

// parameters without default values

string ParameterSetup::MotlName()
{
    if (!motlName.set)
        parameterNotSetError(motlName.name);
    
    return motlName.value;
}

string ParameterSetup::FolderName()
{
    if (!folderName.set)
        parameterNotSetError(folderName.name);

    return folderName.value;
}

string ParameterSetup::RefName()
{
    if (!refName.set)
        parameterNotSetError(refName.name);

    return refName.value;
}

string ParameterSetup::MaskName()
{
    if (!maskName.set)
        parameterNotSetError(maskName.name);

    return maskName.value;
}

string ParameterSetup::CcMaskName()
{
    if (!ccMaskName.set)
        parameterNotSetError(ccMaskName.name);

    return ccMaskName.value;
}

string ParameterSetup::WedgeDataName()
{
    if (!wedgeDataName.set)
        parameterNotSetError(wedgeDataName.name);

    return wedgeDataName.value;
}

string ParameterSetup::WedgeMaskName()
{
    if (!wedgeMaskName.set)
        parameterNotSetError(wedgeMaskName.name);

    return wedgeMaskName.value;
}

string ParameterSetup::SubtomoName()
{
    if (!subtomoName.set)
        parameterNotSetError(subtomoName.name);

    return subtomoName.value;
}

string ParameterSetup::TomoName()
{
    if (!tomoName.set)
        parameterNotSetError(tomoName.name);

    return tomoName.value;
}


int ParameterSetup::NumberOfIterations()
{
    if (!numberOfIterations.set)
        parameterNotSetError(numberOfIterations.name);

    return numberOfIterations.value;
}

int ParameterSetup::StartingIndex()
{
    if (!startingIndex.set)
        parameterNotSetError(startingIndex.name);

    return startingIndex.value;
}

size_t ParameterSetup::VolumeSize()
{
    if (!volumeSize.set)
        parameterNotSetError(volumeSize.name);

    return volumeSize.value;
}


vector<float> ParameterSetup::ConeAngle()
{
    if (!coneAngle.set)
            parameterNotSetError(coneAngle.name);

    return coneAngle.value;
}
vector<float> ParameterSetup::ConeSampling()
{
    if (!coneSampling.set)
            parameterNotSetError(coneSampling.name);

    return coneSampling.value;
}

vector<float> ParameterSetup::InplaneAngle()
{
    if (!inplaneAngle.set)
            parameterNotSetError(inplaneAngle.name);

    return inplaneAngle.value;
}
vector<float> ParameterSetup::InplaneSampling()
{
    if (!inplaneSampling.set)
            parameterNotSetError(inplaneSampling.name);

    return inplaneSampling.value;
}

vector<float> ParameterSetup::LowPass()
{
    if (!lowPass.set)
        parameterNotSetError(lowPass.name);

    return lowPass.value;
}

vector<float> ParameterSetup::HighPass()
{
    if (!highPass.set)
        parameterNotSetError(highPass.name);

    return highPass.value;
}

unsigned int ParameterSetup::TomoDigits()
{
    if (!tomoDigits.set)
        parameterNotSetError(tomoDigits.name);

    return tomoDigits.value;
}

bool ParameterSetup::UseRosemanCC()
{
    return useRosemanCC;
}

bool ParameterSetup::CopySubtomograms()
{
    return copySubtomograms;
}

string ParameterSetup::TempSubtomoFolder()
{
    if (!tempSubtomoFolder.set)
        parameterNotSetError(tempSubtomoFolder.name);
        
    return tempSubtomoFolder.value;
}

unsigned int ParameterSetup::RealStartIndex()
{
    return realStartIndex;
}
unsigned int ParameterSetup::RotationAndFilterIndex()
{
    return rotationAndFilterIndex;
}

string ParameterSetup::MotlBaseName()
{
    return motlBaseName;
}
string ParameterSetup::WedgeBaseName()
{
    return wedgeBaseName;
}
string ParameterSetup::RefBaseName()
{
    return refBaseName;
}

bool ParameterSetup::ChangingBandPass()
{
    return changingBandPass;
}

bool ParameterSetup::ChangingRotations()
{
    return changingRotations;
}

string ParameterSetup::CopySubtomoName()
{
    return copySubtomoName;
}

bool ParameterSetup::SplitIntoEvenOdd()
{
    return splitIntoEvenOdd;
}

string ParameterSetup::SplitType()
{
    return splitType;
}

string ParameterSetup::RefOriginalName()
{
    return refOriginalName;
}

string ParameterSetup::MotlOriginalName()
{
    return motlOriginalName;
}

string ParameterSetup::FscMask()
{
    if (!fscMask.set)
        parameterNotSetError(fscMask.name);

    return fscMask.value;
}

float ParameterSetup::PixelSize()
{
    if (!pixelSize.set)
           parameterNotSetError(pixelSize.name);

    return pixelSize.value;
}

string ParameterSetup::LogFileName()
{
    return logFileName;
}

float ParameterSetup::MotlBinFactor()
{
    return motlBinFactor;
}

float ParameterSetup::DistanceThreshold()
{
    if (!distanceThreshold.set)
        parameterNotSetError(distanceThreshold.name);

    return distanceThreshold.value;
}

bool ParameterSetup::CleanByDistance()
{
    return cleanByDistance;
}

bool ParameterSetup::CleanOutOfBoundsParticles()
{
    return cleanOutOfBoundsParticles;
}


bool ParameterSetup::RenumberParticles()
{
    return renumberParticles;
}

bool ParameterSetup::CleanByMeanGreyValue()
{
    return cleanByMeanGreyValue;
}

unsigned int ParameterSetup::WedgeType()
{
    return wedgeType;
}

bool ParameterSetup::CreateWedgeMasks()
{
    return createWedgeMasks;
}

bool ParameterSetup::RunSA()
{
    return runSA;
}

unsigned int ParameterSetup::InterpolationType()
{
    return interpolationType;
}

float ParameterSetup::InterpolationSigma()
{
    return interpolationSigma;
}

bool ParameterSetup::ComputeFSC()
{
	return computeFSC;
}
