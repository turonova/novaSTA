#pragma once

#include <string>
#include "defines.h"
#include "mrcHeader.h"
#include <fstream>
#include <vector>
#include <typeinfo>
#include "vectorRoutines.h"

using namespace std;

class MRCStack
{
public:

	MRCStack(string aFileName, bool getRange, bool swapYZ, bool keepOpen);

	~MRCStack();

	void readProjections(float* aData, unsigned int numberOfProjections);
	void readProjections(float* aData, unsigned int numberOfProjections, unsigned int aIndex);
	void readAllProjections(float* aData);
	void getSubvolume(float* aData,vector<size_t> startingCoordinates, vector<size_t> subvolumeSize);

	void writeOutHeader();
	float getStackMin();
	float getStackMax();
	float getStackMean();
	float getInputMeanValue();
	unsigned int getNumberOfProjections();
	size_t getProjectionSize();
	vector<unsigned int> getResolution();

	MRCHeader getStackHeader();
	char* getExtraData();

	ifstream mStackFile;

private:

	MRCHeader mHeaderMRC;
	char* extraData;

	size_t mProjectionSize;
	string mStackName;

	template <typename voxelType>
	void readSlices(float* aData, size_t elementsToRead);

	template <typename voxelType>
	void readSubvolume(float* aData, vector<size_t> startingCoordinates, vector<size_t> subvolumeSize);

	void prepareFilePosition(size_t offset = 0);

	template <typename voxelType>
	void getDataRange(ifstream& infile);

	float mDataMin;				//min density value in the stack
	float mDataMax;				//max denisty value in the stack
	float mDataMean;
	float mInputMeanValue;
	unsigned int mNumberOfProjections;
	size_t mSizeOfVoxelType;

	bool mKeepOpen;

	vector<unsigned int> resolution;

    template<typename T>
    void swap(T& a, T& b)
    {
        T tmpValue = b;
        b = a;
        a = tmpValue;
    }

};
