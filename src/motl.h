#pragma once

#include <math.h>
#include <string>
#include <vector>
#include "emFile.h"
#include "parameterSetup.h"
#include "mpiRoutines.h"
#include "motlEntry.h"
#include "structureGeometry.h"

using namespace std;

class Motl
{
public:
    Motl(string motlName);

    Motl(string motlName, int numberOfParts, int partNumber);
    Motl(string motlName, unsigned int featureType);
    Motl(string motlName, unsigned int splitFeatureIndex, int numberOfParts, int partNumber);
    ~Motl(){};

    vector<MotlEntry> getData();
    vector<MotlEntry>& Data();
    void setData(const vector<MotlEntry>& inputData);

    vector<float> getDataAsVector();
    void getDataAsVector(vector<float>& outputData);
    static vector<float> getDataAsVector(vector<MotlEntry>& inputData);

    size_t getOriginalSize();
    size_t getNumberOfParticles();
    size_t getNumberOfFeatures();
    size_t getNumberOfParticlesPerFeatures(size_t featureID);
    size_t getNumberOfAllParticles();
    size_t getOriginalNumberOfParticles();

    void cleanByDistance(float distanceThreshold);
    void cleanClass(float classNumber);

    static void convertData(vector<float>& inputData, vector<MotlEntry>& outputData);
    static void convertData(vector<MotlEntry>& inputData, vector<float>& outputData);

    void splitIntoEvenOdd(vector<MotlEntry>& evenMotl, vector<MotlEntry>& oddMotl);
    void splitIntoEvenOdd(vector<float>& evenMotl, vector<float>& oddMotl);

    void unbin(float binningFactor);

    static void combineEvenOdd(vector<MotlEntry>& evenMotl, vector<MotlEntry>& oddMotl, vector<MotlEntry>& output);
    static void write(vector<MotlEntry>& inputData, string outputName);
    static void write(vector<float>& inputData, string outputName);
    static void read(vector<MotlEntry>& outputData, string inputName);

    void read(string inputName);
    void write(string outputName);
    void writeFullData(string outputName,bool sortBySubtomoNumber=false);

    void setAngles(size_t midx, float phi, float psi, float theta);
    void setShifts(size_t midx, float xs, float ys, float zs);
    void setClass(size_t midx, float classNumber);
    void setMetric(size_t midx, float metricValue);

    size_t getActiveParticles(float classNumber);

    void renumberParticles(bool renumber);
    float cleanByMeanGreyValue();

    void computeAngularDistance(StructureGeometry& structureGeom);
  //  vector<size_t> getAngularStats(StructureGeometry& structureGeom);

    void unifyClassNumber(int newClassNumber);

private:

    void splitDataByFeature(string motlName, unsigned int splitFeatureIndex, int numberOfParts, int partNumber);

    template <typename T>
    void computeEvenDistribution(vector<T>& wholeData, typename vector<T>::iterator& itStart, typename vector<T>::iterator& itEnd, int numberOfParts, int partNumber);

    size_t originalSize;
    size_t originalNumberOfParticles;
    size_t numberOfParticles;

    vector<MotlEntry> data;
    vector<float> dataAsVector;

    size_t numberOfFeatures;
    vector<size_t> numberOParticlesPerFeature;
};
