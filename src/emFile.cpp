#include "emFile.h"
#include "exception.h"
#include <algorithm>

using namespace std;

template <class T>
void swapBytes(T* objp)
{
    unsigned char *memp = reinterpret_cast<unsigned char*>(objp);
    std::reverse(memp, memp + sizeof(T));
}

emFile::emFile(string fileName)
{

	ifstream inputFile;
    
	inputFile.open(fileName.c_str(), std::ios::binary);

    if (!inputFile.good())
    {
        throw ExceptionFileOpen(fileName);
    } 

    char byteHeader[512];

    inputFile.read((char*)&byteHeader, sizeof(emHeader));

    header.machineCode = byteHeader[0];

    // big endian coding - requires swap
    if(header.machineCode == 0 || header.machineCode == 3 || header.machineCode == 5)
    {
        header.genPurpose = byteHeader[1];
        header.notUsed = byteHeader[2];
        header.dataType = byteHeader[3];

        inputFile.read(reinterpret_cast<char*>(&header.x), sizeof(int));
        swapBytes<int>(&header.x);

        inputFile.read(reinterpret_cast<char*>(&header.y), sizeof(int));
        swapBytes<int>(&header.y);

        inputFile.read(reinterpret_cast<char*>(&header.z), sizeof(int));
        swapBytes<int>(&header.z);

        inputFile.read(reinterpret_cast<char*>(&header.extraData), 496);

        data.resize(header.x*header.y*header.z);
        
        switch (header.dataType)
        {
            case 1: readSwappedData<char>(inputFile, header.x*header.y*header.z);
                break;
            case 2: readSwappedData<short int>(inputFile, header.x*header.y*header.z);
                break;
            case 4: readSwappedData<long int>(inputFile, header.x*header.y*header.z);
                break;
            case 5: readSwappedData<float>(inputFile, header.x*header.y*header.z);
                break;
            case 9: readSwappedData<double>(inputFile, header.x*header.y*header.z);
                break;
            default: cout << "Data type of the following EM file is not supported: " << fileName << endl;
                break;
        }
        
    }
    else if (header.machineCode == 4)  // sun machine - not supported
    {
        cout << "The format of following EM file is 4 and thus not supported:" << fileName << endl;
        inputFile.close();
        exit(EXIT_FAILURE);
    }
    else //little endian, can be read normally
    {
        inputFile.clear();
        inputFile.seekg(0);

		inputFile.read((char*)&header, sizeof(emHeader));

        data.resize(header.x*header.y*header.z);

        switch (header.dataType)
        {
            case 1: readData<char>(inputFile, header.x*header.y*header.z);
                break;
            case 2: readData<short int>(inputFile, header.x*header.y*header.z);
                break;
            case 4: readData<long int>(inputFile, header.x*header.y*header.z);
                break;
            case 5: readData<float>(inputFile, header.x*header.y*header.z);
                break;
            case 9: readData<double>(inputFile, header.x*header.y*header.z);
                break;
            default: cout << "Data type of the following EM file is not supported: " << fileName << endl;
                break;
        }
    }

    inputFile.close();
}


emFile::emFile(int x, int y, int z, vector<float>& inData)
{
    header.dataType = 5;
    header.machineCode = 6;

    header.notUsed = 0;
    header.genPurpose = 0;

    header.x = x;
    header.y = y;
    header.z = z;

    fill_n(header.extraData, 496, 0);

    data.resize(inData.size());
    copy(std::begin(inData), std::end(inData), std::begin(data));
}

emFile::emFile(int x, int y, int z)
{
    header.dataType = 5;
    header.machineCode = 6;

    header.notUsed = 0;
    header.genPurpose = 0;

    header.x = x;
    header.y = y;
    header.z = z;

    fill_n(header.extraData, 496, 0);

    data.resize((size_t)x*y*z);
}

emFile::~emFile()
{

}

void emFile::read(string fileName, vector<float>& aData)
{
    emFile* file = new emFile(fileName);
    aData.resize(file->Data().size());

    copy(file->Data().begin(), file->Data().end(), aData.begin());

    delete file;
}

emHeader& emFile::Header()
{
    return header;
}

void emFile::write(string outputName)
{
    std::ofstream outfile;

    try
    {
        outfile.open(outputName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

        if (!outfile.good())
        {
            throw new ExceptionFileOpen(outputName);
        }

        header.dataType = 5;
        header.machineCode = 6;

        outfile.write(reinterpret_cast<char*>(&header), sizeof(emHeader));
        
        outfile.write((const char*)&data[0], data.size()*sizeof(float));

        outfile.close();

    }
    catch (ExceptionFileOpen& e)
    {
        cout << e.what() << endl;
    }
    catch (ExceptionFileFormat& e)
    {
        cout << e.what() << endl;
    }
    catch (...)
    {
        cout << "Error while writing volume to disk!" << endl;
    }
}

void emFile::write(string outputName, vector<float>& inData, size_t dimx, size_t dimy, size_t dimz)
{
    std::ofstream outfile;

    try
    {
        outfile.open(outputName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

        if (!outfile.good())
        {
            throw new ExceptionFileOpen(outputName);
        }

        emHeader tempHeader;

        tempHeader.dataType = 5;
        tempHeader.machineCode = 6;

        tempHeader.notUsed = 0;
        tempHeader.genPurpose = 0;

        tempHeader.x = dimx;
        tempHeader.y = dimy;
        tempHeader.z = dimz;

        fill_n(tempHeader.extraData, 496, 0);

        outfile.write(reinterpret_cast<char*>(&tempHeader), sizeof(emHeader));

        outfile.write((const char*)&inData[0], inData.size()*sizeof(float));

        outfile.close();

    }
    catch (ExceptionFileOpen& e)
    {
        cout << e.what() << endl;
    }
    catch (ExceptionFileFormat& e)
    {
        cout << e.what() << endl;
    }
    catch (...)
    {
        cout << "Error while writing volume to disk!" << endl;
    }
}

template <typename voxelType>
void emFile::readSwappedData(ifstream& inputFile, size_t numberOfVoxels)
{
    voxelType voxelValue;

    for (size_t i = 0; i<numberOfVoxels; i++)
    {
        inputFile.read(reinterpret_cast<char*>(&voxelValue), sizeof(voxelType));
        swapBytes<voxelType>(&voxelValue);
        data[i] = (float)voxelValue;
    }
}

template <typename voxelType>
void emFile::readData(ifstream& inputFile, size_t numberOfVoxels)
{
    voxelType *newData = new voxelType[numberOfVoxels];

    inputFile.read((char*)(newData), numberOfVoxels*sizeof(voxelType));

    for (size_t i = 0; i<numberOfVoxels; i++)
    {
        data[i] = (float)newData[i];
    }

    delete[] newData;
}

vector<float>& emFile::Data()
{
    return data;
}
