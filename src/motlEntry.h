#pragma once

#include <math.h>
#include <string>
#include <vector>


using namespace std;

struct MotlEntry
{
    float metricValue;      // typically CC

    float geometricSpec1;   // e.g. for pores: 1==side,2==top,3=bottom
    float geometricSpec2;   // e.g. for pores: number of assymetric subunit

    float subtomoNumber;
    float tomoNumber;

    float featureNumber;
    float subtomoMean;

    float x;
    float y;
    float z;

    float xShift;
    float yShift;
    float zShift;

    float free1;
    float free2;
    float free3;

    float phi;
    float psi;
    float theta;

    float classNumber;

    MotlEntry(){}

    MotlEntry(vector<float>::iterator& data)
    {
        metricValue = *data;
        geometricSpec1 = *(data + 1);
        geometricSpec2 = *(data + 2);
        subtomoNumber = *(data + 3);
        tomoNumber = *(data + 4);

        featureNumber = *(data + 5);

        subtomoMean = *(data + 6);

        x = *(data + 7);
        y = *(data + 8);
        z = *(data + 9);

        xShift = *(data + 10);
        yShift = *(data + 11);
        zShift = *(data + 12);

        free1 = *(data + 13);
        free2 = *(data + 14);
        free3 = *(data + 15);

        phi = *(data + 16);
        psi = *(data + 17);
        theta = *(data + 18);

        classNumber = *(data + 19);
    }

    void convertToVector(vector<float>::iterator& data)
    {
        *data = metricValue;

        *(data + 1) = geometricSpec1;
        *(data + 2) = geometricSpec2;
        *(data + 3) = subtomoNumber;
        *(data + 4) = tomoNumber;

        *(data + 5) = featureNumber;
        *(data + 6) = subtomoMean;

        *(data + 7) = x;
        *(data + 8) = y;
        *(data + 9) = z;

        *(data + 10) = xShift;
        *(data + 11) = yShift;
        *(data + 12) = zShift;

        *(data + 13) = free1;
        *(data + 14) = free2;
        *(data + 15) = free3;

        *(data + 16) = phi;
        *(data + 17) = psi;
        *(data + 18) = theta;

        *(data + 19) = classNumber;
    }

    bool operator < (const MotlEntry& entry) const
    {
        return (featureNumber < entry.featureNumber);
    }

    static bool sortByFeature(const MotlEntry& entry1,const MotlEntry& entry2)
    {
        return (entry1.featureNumber < entry2.featureNumber);
    }

    static bool sortBySubtomoNumber(const MotlEntry& entry1,const MotlEntry& entry2)
    {
        return (entry1.subtomoNumber < entry2.subtomoNumber);
    }


    float computeDistance(MotlEntry& entry)
    {
        float distance= (x+xShift-entry.x-entry.xShift)*(x+xShift-entry.x-entry.xShift) +
                        (y+yShift-entry.y-entry.yShift)*(y+yShift-entry.y-entry.yShift) +
                        (z+zShift-entry.z-entry.zShift)*(z+zShift-entry.z-entry.zShift);

        return sqrt(distance);
    }
};
