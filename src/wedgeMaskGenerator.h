#pragma once

#include <vector>
#include <map>
#include <unordered_set>
#include "volume.h"
#include "wedgeMask.h"
#include "fftRoutines.h"
#include "vectorRoutines.h"
#include "emFile.h"
#include "parameterSetup.h"

using namespace std;


class WedgeMaskGenerator
{
public:
    
    static WedgeMaskGenerator* create(ParameterSetup& iParams, int numberOfParts, int partNumber);

    WedgeMaskGenerator(ParameterSetup& iParams);

    virtual void generateMask() = 0;

    static string createNameWithZeroDigits(string baseName, unsigned int tomoNumber, string extension, int aTomoDigits);

	void getTomoNumbersPos(map<unsigned int,size_t>& aTomoNumbersPos);

	void getTomoNumbers(vector<unsigned int>& aTomoNumbers);

protected:

    void getTomoNumbers();
    
    ParameterSetup params;

    vector<unsigned int> tomoNumbers;
	map<unsigned int,size_t> tomoNumbersPos;
    int tomoDigits;
    vector<float> inputData;

    vector<float> mask;
    size_t maskSize;

    size_t entryOffset;

    string wedgeInputName;
    string wedgeOutputName;

    vector<unsigned int>::iterator tomoListStart;
    vector<unsigned int>::iterator tomoListEnd;
};


class wmg_Bowtie : public WedgeMaskGenerator
{
public:
    
    wmg_Bowtie(ParameterSetup& iParams, int numberOfParts, int partNumber);
    ~wmg_Bowtie(){};

    void generateMask();

private:

    void parseWedgeFile();

};


class wmg_Star : public WedgeMaskGenerator
{
public:
    wmg_Star(ParameterSetup& iParams, int numberOfParts, int partNumber);
    ~wmg_Star(){};

    void generateMask();

private:

    void parseWedgeFile();

};

class wmg_SingleMask : public WedgeMaskGenerator
{
public:
    wmg_SingleMask(ParameterSetup& iParams, int numberOfParts, int partNumber);
    ~wmg_SingleMask(){};

    void generateMask();

private:

    void parseWedgeFile();

};

class wmg_InputMasks : public WedgeMaskGenerator
{
public:
    wmg_InputMasks(ParameterSetup& iParams, int numberOfParts, int partNumber);
    ~wmg_InputMasks(){};

    void generateMask();

private:

    bool getStartingCoordinates(const vector<float>& centerCoord, vector<int>& startCoord, vector<unsigned int>& resolution);

    void parseWedgeFile();

};

