#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <typeinfo>

using namespace std;

struct emHeader{
    
    char machineCode;
    char genPurpose;
    char notUsed;
    char dataType;

    int x;
    int y;
    int z;

    char extraData[496];
};

class emFile
{
public:
    emFile(string fileName);

    emFile(int x, int y, int z, vector<float>& inData);

    emFile(int x, int y, int z);

    ~emFile();

    void write(string outputName);
    vector<float>& Data();
    emHeader& Header();

    static void write(string outputName, vector<float>& inData, size_t dimx, size_t dimy, size_t dimz);
    static void read(string fileName, vector<float>& aData);

private:

    template <typename voxelType>
    void readSwappedData(ifstream& inputFile, size_t numberOfVoxels);
    
    template <typename voxelType>
    void readData(ifstream& inputFile, size_t numberOfVoxels);

    
    emHeader header;
    vector<float> data;


};
