#include "structureGeometry.h"
#include "vectorRoutines.h"
#include "quaternion.h"
#include "mpiRoutines.h"
#include <sstream>


#include <algorithm>

using namespace std;
using namespace vr;

StructureGeometry* StructureGeometry::create(ParameterSetup& iParams)
{
    switch (iParams.StructGeomType())
    {
    case 0: return new sg_General(iParams);
            break;
    case 1: return new sg_Sphere(iParams);
            break;
    case 2: return new sg_Tube(iParams);
        break;
    case 3: return new sg_NPC(iParams);
        break;
    default:
        cout << "ERROR - specified structure geometry is not known!" << endl;
        exit(EXIT_FAILURE);
    }
}

string StructureGeometry::getDescription()
{
    computeStats();
    return createDescription();
}

vector<size_t> StructureGeometry::getStats()
{
    computeStats();
    return stats;
}

void StructureGeometry::reduceStats(const vector<size_t>& partialStats)
{
    stats.clear();
    stats.resize(partialStats.size());

    for (size_t i = 0; i < stats.size(); i++)
    {
        stats[i] = MPIRoutines::reduceSizeTNumber(partialStats[i]);
    }
}

template <typename T>
string StructureGeometry::createStringFromVector(vector<T> inputValues, string valueSeparator)
{
    stringstream  vectorString;

    vectorString << inputValues[0];

    for (size_t i = 1; i < inputValues.size(); i++)
    {
        vectorString << valueSeparator << inputValues[i];
    }

    return vectorString.str();
}


sg_General::sg_General(ParameterSetup& iParams) : StructureGeometry(iParams)
{

}

sg_Sphere::sg_Sphere(ParameterSetup& iParams) : StructureGeometry(iParams)
{

}

sg_Tube::sg_Tube(ParameterSetup& iParams) : StructureGeometry(iParams)
{

}

sg_NPC::sg_NPC(ParameterSetup& iParams) : StructureGeometry(iParams)
{
    expectedFullDistance = 45.0f;
    expectedNumberOfSubunits = 8;
    subunitHalf = 4;
    minNumberOfSUS = 3;

    nCorrectedSUS = 0;
    nNotCorrectedSUS = 0;
    nMinRemovedSUS = 0;

    thresholds = iParams.AngularDistanceThresholds();

    if (thresholds.size() == 3)
    {
        newClassNumber = (int) thresholds[2];
        thresholds.pop_back();
    }
    else
        newClassNumber = iParams.SubtomoClass();

    cleanType = iParams.ComputeAngularDistance();
}

void sg_NPC::initStructure(vector<MotlEntry>::iterator aMotlStart, size_t numberOfParticles)
{
    motlStart = aMotlStart;
    motlEnd = motlStart + numberOfParticles;

    actualNubmberOfSubunits = numberOfParticles;
    angleMeasurements = actualNubmberOfSubunits - 1.0f;
    
    normals.resize(actualNubmberOfSubunits);
    quats.resize(actualNubmberOfSubunits);
    subunitID.resize(actualNubmberOfSubunits);
    
    //for (unsigned int i = 0; i < actualNubmberOfSubunits; i++)
    unsigned int i = 0;
    for (motl  = motlStart; motl != motlEnd; motl++)
    {
        normals[i] = Quaternion::getNormalVectorFromEulerAngles((*motl).psi, (*motl).theta);
        quats[i] = Quaternion((*motl).phi, (*motl).psi, (*motl).theta);
        subunitID[i] = (int)(*motl).geometricSpec2;
        i++;
    }
}

void sg_NPC::clearStructure()
{
    normals.clear();
    quats.clear();
    subunitID.clear();
    badSubunits.clear();

}

void sg_NPC::computeAngularDistance()
{
    vector<float> normalAngularDist(actualNubmberOfSubunits, 0.0f);
    vector<float> fullAngularDist(actualNubmberOfSubunits, 0.0f);

    for (unsigned int i = 0; i < actualNubmberOfSubunits; i++)
    {
        for (unsigned int j = 0; j < i; j++)
        {
            float normalDist = suJointAngleNormals(normals[i], normals[j]);

            normalAngularDist[j] += normalDist;
            normalAngularDist[i] += normalDist;

            float fullAngle;
            float fullDist;

            Quaternion::getAngularDistance(quats[i], quats[j], fullAngle, fullDist);

            float angleWeight = suDistanceSmallerAngle(subunitID[i], subunitID[j]);
          
            float weightedAngle = fabs(expectedFullDistance - fullAngle / angleWeight);
            fullAngularDist[j] += weightedAngle;
            fullAngularDist[i] += weightedAngle;
        }
    }

    unsigned int i = 0;

    for (motl = motlStart; motl != motlEnd; motl++)
    {
        (*motl).free1 = normalAngularDist[i] / angleMeasurements;
        (*motl).free2 = fullAngularDist[i] / angleMeasurements;
        i++;
    }

    if (cleanType == 2) // cleaning based on given thresholds
    {
        markMisalignedParticles();
    }
    else if (cleanType == 3) // correcting particles
    {
        correctMisalignedParticles();
    }
}

float sg_NPC::suDistanceSmallerAngle(int su1, int su2)
{
    int subunitDist = abs(su1 - su2);

    float distanceSmallerAngle;
    if (subunitDist <= subunitHalf)
    {
        distanceSmallerAngle = subunitDist;
    }
    else
    {
        distanceSmallerAngle = expectedNumberOfSubunits%subunitDist;
    }

    return distanceSmallerAngle;
}


float sg_NPC::suJointAngleNormals(vector<float>& normal1, vector<float>& normal2)
{
    float jointAngle = acos(vr::dot_product(normal1, normal2))*180.f / M_PI;

    if (isnan(jointAngle))
    {
        jointAngle = 0.0f;
    }

    return jointAngle;
}

void sg_NPC::markMisalignedParticles()
{
    unsigned int i = 0;

    float normalAngle = thresholds[0];
    float fullAngle = thresholds[1];

    size_t nBadSUS = 0;
    size_t nNotCorr = 0;

    for (motl = motlStart; motl != motlEnd; motl++)
    {
        if ((*motl).free1 > normalAngle || (*motl).free2 > fullAngle)
        {
            (*motl).classNumber = -3;
            nNotCorr++;
            nBadSUS++;
        }
    }

    if ((actualNubmberOfSubunits - nBadSUS) < minNumberOfSUS)
    {
        for (motl = motlStart; motl != motlEnd; motl++)
        {
            (*motl).classNumber = -3;
            nMinRemovedSUS++;
        }
    }
    else
    {
        nNotCorrectedSUS += nNotCorr;
    }
}

void sg_NPC::getMisalignedSubunits()
{
    unsigned int i = 0;

    float normalAngle = thresholds[0];
    float fullAngle = thresholds[1];
   
    badSubunits.resize(actualNubmberOfSubunits,0);
    badSubunitsCount = 0;

    for (motl = motlStart; motl != motlEnd; motl++)
    {
        if ((*motl).free1 > normalAngle || (*motl).free2 > fullAngle)
        {
            (*motl).free3 = 1;
            badSubunits[i] = 1;
            badSubunitsCount++;
        }
        else
        {
            (*motl).free3 = 0;
        }
        i++;
    }
}

void sg_NPC::correctMisalignedParticles()
{
    getMisalignedSubunits();
    
    if ((actualNubmberOfSubunits - badSubunitsCount) < minNumberOfSUS)
    {
        for (motl = motlStart; motl != motlEnd; motl++)
        {
            (*motl).classNumber = -3;
            nMinRemovedSUS++;
        }
        return;
    }

    vector<int> interpStart;
    vector<int> interpEnd;
    vector<int> seqLengths;

    getInterpolationSeqences(interpStart, interpEnd, seqLengths);

    unsigned int nBadSUS = 0;
    unsigned int nCorr = 0;
    unsigned int nNotCorr = 0;

    for (int i = 0; i < interpStart.size(); i++)
    {
        int startID = interpStart[i];
        int endID = interpEnd[i];
        int seqLength = seqLengths[i];

        int suDistance = subunitID[endID] - subunitID[startID];

        if (suDistance < 0)
            suDistance = expectedNumberOfSubunits - subunitID[startID] + subunitID[endID];

        float fullAngle, fullDist;
        Quaternion::getAngularDistance(quats[startID], quats[endID], fullAngle, fullDist);

        if (suDistance > 4)
            fullAngle = 360.0f - fullAngle;

        bool keepLargeAngle = false;
        if (fullAngle > 180)
            keepLargeAngle = true;

        float samplingStep = 1.0f / suDistance;

        vector<Quaternion> interpQuats;
        Quaternion::slerp(interpQuats, quats[startID], quats[endID], samplingStep, 1.0, keepLargeAngle);

        bool correctSUS = true;

        if (suDistance == 4)
            correctSUS = checkInterpolatedRotations(interpQuats, startID, endID, samplingStep, keepLargeAngle);

        vector<int> correctID;

        for (int j = 1; j < (seqLength+1); j++)
        {
            int idx = (startID + j) % actualNubmberOfSubunits;
            correctID.push_back(idx);
        }


        for (unsigned int j = 0; j < correctID.size(); j++)
        {
            int angleID;

            if (correctID.size() == (interpQuats.size()+2))
                angleID = j+1;
            else
            {
                angleID = subunitID[correctID[j]] - subunitID[startID] ;
                if (angleID < 0)
                    angleID = angleID + expectedNumberOfSubunits;
                // angleID=(subunitID[correctID[i]] - subunitID[startID] + expectedNumberOfSubunits)%+ expectedNumberOfSubunits;
            }

            vector<float> eulerAngles = interpQuats[angleID].getEulerAngles();
            motl = motlStart + correctID[j];
            (*motl).phi = eulerAngles[0];
            (*motl).psi = eulerAngles[1];
            (*motl).theta = eulerAngles[2];

            if (!correctSUS)
            {
                (*motl).classNumber = -3;
                nBadSUS++;
               // nNotCorrectedSUS++;
                nNotCorr++;
            }
            else
            {
               // nCorrectedSUS++;
                nCorr++;
                (*motl).classNumber = newClassNumber;
            }
        }
    }

    if ((actualNubmberOfSubunits - nBadSUS) < minNumberOfSUS)
    {
        for (motl = motlStart; motl != motlEnd; motl++)
        {
            (*motl).classNumber = -3;
            nMinRemovedSUS++;
        }
    }
    else
    {
        nCorrectedSUS += nCorr;
        nNotCorrectedSUS += nNotCorr;
    }
}

bool sg_NPC::checkInterpolatedRotations(vector<Quaternion>& interpQuats, int startQ, int endQ, float samplingStep, bool keepLargeAngle)
{

    vector<Quaternion> interpQuatsOppDir;
    Quaternion::slerp(interpQuatsOppDir, quats[startQ], quats[endQ], samplingStep, 1.0, !keepLargeAngle);
   
    float fullAngle, fullDist;
    float angle1, angle2;

    if (startQ > 0)
    {
        int suDist = subunitID[startQ] + 1 - subunitID[startQ - 1];
        
        Quaternion::getAngularDistance(interpQuats[1], quats[startQ - 1], fullAngle, fullDist);
        angle1 = fullAngle / suDist;

        Quaternion::getAngularDistance(interpQuatsOppDir[1], quats[startQ - 1], fullAngle, fullDist);
        angle2 = fullAngle / suDist;
    }
    else if (endQ < (actualNubmberOfSubunits-1))
    {
        int suDist = subunitID[endQ + 1] - subunitID[endQ] + 1;

        Quaternion::getAngularDistance(interpQuats[interpQuats.size() - 2], quats[endQ + 1], fullAngle, fullDist);
        angle1 = fullAngle / suDist;

        Quaternion::getAngularDistance(interpQuatsOppDir[interpQuats.size() - 2], quats[endQ + 1], fullAngle, fullDist);
        angle2 = fullAngle / suDist;
    }
    
    angle1 = fabs(angle1 - expectedFullDistance);
    angle2 = fabs(angle2 - expectedFullDistance);

    if (fabs(angle1 - angle2) < 20)
        return false;

    if (angle1>angle2)
    {
        for (unsigned int i = 0; i < interpQuats.size(); i++)
            interpQuats[i] = interpQuatsOppDir[i];
    }

    return true;
}

void sg_NPC::getInterpolationSeqences(vector<int>& interpStart, vector<int>& interpEnd, vector<int>& seqLength)
{

    if (badSubunitsCount == 0)
        return;


    for (int i = 0; i < actualNubmberOfSubunits; i++)
    {
        int seqCounter = 1;

        if (badSubunits[i] == 1)
        {
            int k = (i - 1 + actualNubmberOfSubunits) % actualNubmberOfSubunits;

            while (badSubunits[k] == 1)
            {   
                k = (k - 1 + actualNubmberOfSubunits) % actualNubmberOfSubunits;
                seqCounter++;
            }
           
            if (std::find(interpStart.begin(), interpStart.end(), k) == interpStart.end())
                interpStart.push_back(k);

            k = (i + 1) % actualNubmberOfSubunits;

            while (badSubunits[k] == 1)
            {
                k = (k + 1) % actualNubmberOfSubunits;
                seqCounter++;
            }
                

            if (std::find(interpEnd.begin(), interpEnd.end(), k) == interpEnd.end())
            {
                interpEnd.push_back(k);
                seqLength.push_back(seqCounter);
            }
               
        }

    }

    if (interpStart.size() != interpEnd.size())
    {
        cout << endl << "ERROR: Failed to compute interpolation sequence!!!" << endl;
        exit(EXIT_FAILURE);
    }
}

void sg_NPC::computeStats()
{
    vector<size_t> partialStats;

    partialStats.push_back(nMinRemovedSUS);
    partialStats.push_back(nCorrectedSUS);
    partialStats.push_back(nNotCorrectedSUS);

    reduceStats(partialStats);
}

string sg_NPC::createDescription()
{
    stringstream  description;

    if (cleanType == 0)
    {
        description << "";
    }
    else if (cleanType == 1)
    {
        description << "\n\nAngular distance for each particle was computed and stored in motl at rows 14 and 15.";
    }
    else if (cleanType == 2)
    {
        size_t removedParticles = stats[0] + stats[2];
        description << "\n\nParticles in motl were cleaned by angular distance with thresholds" << createStringFromVector(thresholds) << ".\n";
        description << "In total " << removedParticles << " particles were removed.";
    }
    else if (cleanType == 3)
    {
        size_t removedParticles = stats[0] + stats[2];
        description << "\n\nParticles in motl were corrected or cleaned by angular distance with thresholds " << createStringFromVector(thresholds) << ".\n";
        description << "In total " << removedParticles << " particles were removed (" << stats[0] << " did not reach the minimum criterion, " << stats[2]  << " could not be reliably corrected).\n";
        description << "In total " << stats[1] << " particles were corrected and their class number was set to " << newClassNumber << ".";
    }

    return description.str();
}

