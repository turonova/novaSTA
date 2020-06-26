#pragma once

#include <string>
#include <vector>
#include <map>
#include "emFile.h"
#include "parameterSetup.h"
#include "mrcStack.h"
#include "motl.h"
#include "log.h"


using namespace std;

class ExtractSubtomos
{
public:
    ExtractSubtomos(ParameterSetup& iParams, int iProcessID, int iNumberOfProcesses);
    ExtractSubtomos(ParameterSetup& iParams, vector<MotlEntry>& iData, int iProcessID);

    ExtractSubtomos(ParameterSetup& iParams, Log& iLog, vector<MotlEntry>& iData, int iProcessID);

    ~ExtractSubtomos(){};

    void run();

    vector<MotlEntry> getMotlData();

private:

    bool createSubtomogram(size_t subtomoNumber,vector<int>& coord);
    vector<int> startPosition(vector<int>& coord );
    void checkBoundaries(vector<size_t>& realSize,vector<int> &realShift,vector<int> &startCoord,bool& cutSubtomo);
    //void tomoExist(unsigned int tomoNumber);
    void prepareStatsFile();

    string getSubtomoFolder(string path);

    bool subtomoExist(unsigned int subtomoNumber);

    void normalizeSubtomo(vector<float> &subtomogram);
    void getSubtomoStats(vector<float>& subtomogram, float subtomoNumber);
   // void writeSubtomogramStats();
    void extractSubtomograms();

    string createTomoName(float recordNumber,string extension);
    string createSubtomoName(size_t subtomoNumber);

    Log logger;
    ParameterSetup params;
    int processID;
    int numberOfProcesses;
    int subtomoSize;

    vector<MotlEntry> motlData;
    vector<unsigned int> resolution;
    MRCStack* volume;
    //vector<vector <float> > subtomoStats;
    //int statsIndex;
    void extractSubtomosPrep();

    unsigned int currentTomoNumber ;
    float currentSubtomoMean;

    ofstream statsFile;

};
