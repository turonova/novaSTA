#pragma once

#include "parameterSetup.h"
#include <vector>
#include "emFile.h"
#include "volume.h"
#include "motlEntry.h"
#include "log.h"

using namespace std;

enum StructureGeometryType
{
    general,        // geometrically not defined structure
    sphere,         // structure forming spherical object, i.e. immature HIV
    tube,           // structure forming tubular objects, i.e. MinCD
    npc             // nuclear pore complex
};

class StructureGeometry
{
public:

    static StructureGeometry* create(ParameterSetup& iParams);

    StructureGeometry(ParameterSetup& iParams){};

    virtual void initStructure(vector<MotlEntry>::iterator aMotlStart, size_t numberOfParticles) = 0;
    virtual void computeAngularDistance() = 0;
    virtual void clearStructure()=0;
   
    string getDescription();
    vector<size_t> getStats();

protected:

    virtual void computeStats()=0;
    virtual string createDescription() = 0;
    
    void reduceStats(const vector<size_t>& partialStats);

    template <typename T>
    string createStringFromVector(vector<T> inputValues, string valueSeparator = ", ");

    vector<MotlEntry>::iterator motlStart;
    vector<MotlEntry>::iterator motlEnd;
    vector<MotlEntry>::iterator motl;

    vector<size_t> stats;

};

class sg_General : public StructureGeometry
{
public:

    sg_General(ParameterSetup& iParams);
    ~sg_General(){};

    void initStructure(vector<MotlEntry>::iterator aMotlStart, size_t numberOfParticles){};
    void computeAngularDistance(){};
    void clearStructure(){};
    void computeStats() {};
    string createDescription() {};
};

class sg_Sphere : public StructureGeometry
{
public:

    sg_Sphere(ParameterSetup& iParams);
    ~sg_Sphere(){};

    void initStructure(vector<MotlEntry>::iterator aMotlStart, size_t numberOfParticles){};
    void computeAngularDistance(){};
    void clearStructure(){};
    void computeStats() {};
    string createDescription() {};

};

class sg_Tube : public StructureGeometry
{
public:

    sg_Tube(ParameterSetup& iParams);
    ~sg_Tube(){};

    void initStructure(vector<MotlEntry>::iterator aMotlStart, size_t numberOfParticles){};
    void computeAngularDistance(){};
    void clearStructure(){};
    void computeStats() {};
    string createDescription() {};
};

class sg_NPC : public StructureGeometry
{
public:
    
    sg_NPC(ParameterSetup& iParams);
    ~sg_NPC();

    void initStructure(vector<MotlEntry>::iterator aMotlStart, size_t numberOfParticles);
    void computeAngularDistance();
    void clearStructure();
    void computeStats();
    string createDescription();

protected:

    void markMisalignedParticles();
    void correctMisalignedParticles();

    void getMisalignedSubunits();
    void getInterpolationSeqences(vector<int>& interpStart, vector<int>& interpEnd, vector<int>& seqLength);

    float suDistanceSmallerAngle(int su1, int su2);
    float suJointAngleNormals(vector<float>& normal1, vector<float>& normal2);

    bool checkInterpolatedRotations(vector<Quaternion>& interpQuats, int startQ, int endQ, float samplingStep, bool keepLargeAngle);

    float expectedFullDistance;
    int expectedNumberOfSubunits;
    int subunitHalf;
    int actualNubmberOfSubunits;
    float angleMeasurements;
    unsigned int badSubunitsCount;

    size_t nCorrectedSUS;
    size_t nNotCorrectedSUS;
    size_t nMinRemovedSUS;

    vector<vector<float>> normals;
    vector<Quaternion> quats;
    vector<int> subunitID; 
    vector<int> badSubunits;
    unsigned int minNumberOfSUS;

    vector<float> thresholds; 
    unsigned int cleanType;
    int newClassNumber;
};
