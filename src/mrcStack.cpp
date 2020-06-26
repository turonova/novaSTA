#include <stdio.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include <math.h>
#include <cstring>
#include <algorithm>

#include "exception.h"
#include "mrcStack.h"

using namespace std;


	MRCStack::MRCStack(string aFileName, bool getRange, bool swapYZ, bool keepOpen)
	{
		mKeepOpen = keepOpen;

		mStackName = aFileName;
		mStackFile.open(mStackName.c_str(), std::ios::binary);

		if (!mStackFile.good())
		{
			 throw ExceptionFileOpen(mStackName);
		}

		mStackFile.read((char*)&mHeaderMRC, sizeof(MRCHeader));
		mStackFile.ignore(mHeaderMRC.extra);

		if(swapYZ)
		{
			swap<int>(mHeaderMRC.ny,mHeaderMRC.nz);
			swap<int>(mHeaderMRC.my,mHeaderMRC.mz);
			swap<int>(mHeaderMRC.nyStart,mHeaderMRC.nzStart);
			swap<float>(mHeaderMRC.cellDimY,mHeaderMRC.cellDimZ);
			swap<float>(mHeaderMRC.cellAngleY,mHeaderMRC.cellAngleZ);
			swap<int>(mHeaderMRC.mapR,mHeaderMRC.mapS);
		}

		mInputMeanValue = mHeaderMRC.dMean;
		mNumberOfProjections = mHeaderMRC.nz;
		mProjectionSize = (size_t)mHeaderMRC.nx * (size_t)mHeaderMRC.ny;
		resolution.resize(3);
		resolution[0] = mHeaderMRC.nx;
		resolution[1] = mHeaderMRC.ny;
		resolution[2] = mHeaderMRC.nz;

		//writeOutHeader();

		if(getRange)
		{
			switch(mHeaderMRC.mode)
			{
				case 0: getDataRange<unsigned char>(mStackFile);
						break;
				case 1: getDataRange<signed short>(mStackFile);
						break;
				case 2: getDataRange<float>(mStackFile);
						break;
				case 6: getDataRange<unsigned short>(mStackFile);
						break;
				default: throw ExceptionFileFormat();
						break;
			}
		}

		//header initialization - holds only for parallel geometry !!!

		switch(mHeaderMRC.mode)
		{
			case 0: mSizeOfVoxelType=sizeof(unsigned char);
					break;
			case 1: mSizeOfVoxelType=sizeof(signed short);
					break;
			case 2: mSizeOfVoxelType=sizeof(float);
					break;
			case 6: mSizeOfVoxelType=sizeof(unsigned short);
					break;
			default: throw ExceptionFileFormat();
					break;
		}

		if(mHeaderMRC.extra!=0)
		{
			extraData=new char[mHeaderMRC.extra];
			mStackFile.seekg(sizeof(MRCHeader));
			mStackFile.read(extraData, mHeaderMRC.extra);
		}
		else
		{
			extraData=new char[1];	//dummy
		}

		if(mKeepOpen)
		{
			mStackFile.seekg(sizeof(MRCHeader) + mHeaderMRC.extra);
		}
		else
			mStackFile.close();

	}

	MRCStack::~MRCStack()
	{
		if(mKeepOpen)
			mStackFile.close();

		delete[] extraData;
	}

	void MRCStack::writeOutHeader()
	{
		 cout << "Number of columns, rows, sections"  << mHeaderMRC.nx << ", " << mHeaderMRC.ny << ", " << mHeaderMRC.nz << endl;
		 cout << "Map mode" << mHeaderMRC.mode << endl;
		 cout << "Start columns, rows, sections"  << mHeaderMRC.nxStart << ", " << mHeaderMRC.nyStart << ", " << mHeaderMRC.nzStart << endl;
		 cout << "Grid size in x, y, z" << mHeaderMRC.mx << ", " << mHeaderMRC.my << ", " << mHeaderMRC.mz << endl;
		 cout << "Cell dimensions in x, y, z" << mHeaderMRC.cellDimX << ", " << mHeaderMRC.cellDimY << ", " << mHeaderMRC.cellDimZ << endl;
		 cout << "Cell angles in x, y, z (degrees)"  << mHeaderMRC.cellAngleX << ", " << mHeaderMRC.cellAngleY << ", " << mHeaderMRC.cellAngleZ << endl;
		 cout << "Axis corresponding to columns, rows and sections (1,2,3 for X,Y,Z)" << mHeaderMRC.mapC << ", " << mHeaderMRC.mapR << ", " << mHeaderMRC.mapS << endl;
		 cout << "Origin on x, y, z"<<  mHeaderMRC.originX << ", " << mHeaderMRC.originY << ", " << mHeaderMRC.originZ << endl;
		 cout << "Minimum density" <<   mHeaderMRC.dMin << endl;
		 cout << "Maximum density" <<   mHeaderMRC.dMax << endl;
		 cout << "Mean density" <<  mHeaderMRC.dMean << endl;
		 cout << "Tilt angles - original" <<   mHeaderMRC.tiltangles[0] << ", " << mHeaderMRC.tiltangles[1] << ", " << mHeaderMRC.tiltangles[2] << endl;
		 cout << "Tilt angles - current" <<   mHeaderMRC.tiltangles[3] << ", " << mHeaderMRC.tiltangles[4] << ", " << mHeaderMRC.tiltangles[5] << endl;
		 cout << "Space group" <<   mHeaderMRC.ISPG << endl;
		 cout << "Number of bytes used for symmetry data" <<  mHeaderMRC.NSYMBT << endl;
		 cout << "Number of bytes in extended header" <<  mHeaderMRC.extra << endl;
		 cout << "Creator ID" <<  mHeaderMRC.creatorID << endl;
		 cout << "ID type" <<   mHeaderMRC.idtype << endl;
		 cout << "Lens" <<   mHeaderMRC.lens << endl;
		 cout <<  mHeaderMRC.nLabel << "labels:" << endl;
		 for(int i=0; i<mHeaderMRC.nLabel;i++)
			 cout << mHeaderMRC.labels[i] << endl;
	}

	void MRCStack::readProjections(float* aData, unsigned int numberOfProjections, unsigned int sliceNumber)
	{
		prepareFilePosition(sliceNumber*mProjectionSize*mSizeOfVoxelType);

		switch(mHeaderMRC.mode)
		{
			case 0: readSlices<unsigned char>(aData,mProjectionSize*(size_t)numberOfProjections);
					break;
			case 1: readSlices<signed short>(aData,mProjectionSize*(size_t)numberOfProjections);
					break;
			case 2: readSlices<float>(aData,mProjectionSize*(size_t)numberOfProjections);
					break;
			case 6: readSlices<unsigned short>(aData,mProjectionSize*(size_t)numberOfProjections);
					break;
			default: throw ExceptionFileFormat();
					break;
		}


		if(!mKeepOpen)
			mStackFile.close();
	}

	void MRCStack::getSubvolume(float* aData,vector<size_t> startingCoordinates, vector<size_t> subvolumeSize)
	{
	    switch(mHeaderMRC.mode)
        {
            case 0: readSubvolume<unsigned char>(aData,startingCoordinates,subvolumeSize);
                    break;
            case 1: readSubvolume<signed short>(aData,startingCoordinates,subvolumeSize);
                    break;
            case 2: readSubvolume<float>(aData,startingCoordinates,subvolumeSize);
                    break;
            case 6: readSubvolume<unsigned short>(aData,startingCoordinates,subvolumeSize);
                    break;
            default: throw ExceptionFileFormat();
                    break;
        }
	}

	void MRCStack::readAllProjections(float* aData)
	{
		prepareFilePosition();

		switch(mHeaderMRC.mode)
		{
			case 0: readSlices<unsigned char>(aData,mProjectionSize*mNumberOfProjections);
					break;
			case 1: readSlices<signed short>(aData,mProjectionSize*mNumberOfProjections);
					break;
			case 2: readSlices<float>(aData,mProjectionSize*mNumberOfProjections);
					break;
			case 6: readSlices<unsigned short>(aData,mProjectionSize*mNumberOfProjections);
					break;
			default: throw ExceptionFileFormat();
					break;
		}

		if(!mKeepOpen)
			mStackFile.close();
	}

	void MRCStack::prepareFilePosition(size_t offset)
	{
		if(mKeepOpen)
		{
			if (!mStackFile.good())
			{
				throw ExceptionFileOpen(mStackName);
			}
		}
		else
		{
			mStackFile.open(mStackName.c_str(), std::ios::binary);

			if (!mStackFile.good())
			{
				throw ExceptionFileOpen(mStackName);
			}

			mStackFile.seekg(sizeof(MRCHeader) + mHeaderMRC.extra + offset);

		}
	}

    /*template <typename voxelType>
    float MRCStack::convertValue(voxelType rawValue, float minValue)
    {
        float rawFloat = (float)rawValue;
        if(mLogarithmizeData)
        {
            return logf(mDataMax / (rawFloat - minValue));
        }
        else
        {
            return rawFloat;
        }
    }*/

	template <typename voxelType>
	void MRCStack::readSlices(float* aData, size_t elementsToRead)
	{
        voxelType *data = new voxelType[elementsToRead];

        mStackFile.read((char*)(data), elementsToRead*sizeof(voxelType));

        for(size_t i=0;i<(size_t)elementsToRead;i++)
        {
            aData[i]=(float)data[i];
        }
        delete[] data;
	}

	template <typename voxelType>
	void MRCStack::readSubvolume(float* aData, vector<size_t> startingCoordinates, vector<size_t> subvolumeSize)
	{
	    //mStackFile.open(mStackName.c_str(), std::ios::binary);

        if (!mStackFile.good())
        {
            throw ExceptionFileOpen(mStackName);
        }

	    voxelType *data = new voxelType[subvolumeSize[0]];
        size_t startPosition = startingCoordinates[0]+startingCoordinates[1]*mHeaderMRC.nx+startingCoordinates[2]*mProjectionSize;

        for(size_t z=0; z<subvolumeSize[2];z++)
        {
            for(size_t y=0; y<subvolumeSize[1];y++)
            {
                size_t offset = startPosition + y*mHeaderMRC.nx+z*mProjectionSize;
                mStackFile.seekg(sizeof(MRCHeader) + mHeaderMRC.extra + offset*sizeof(voxelType));
                mStackFile.read((char*)(data), subvolumeSize[0]*sizeof(voxelType));

                for(size_t i=0;i<(size_t)subvolumeSize[0];i++)
                {
                    aData[i+ y*subvolumeSize[0]+z*subvolumeSize[0]*subvolumeSize[1]]=(float)data[i];
                }

            }
        }

        delete[] data;

        //mStackFile.close();
	}


	template <typename voxelType>
	void MRCStack::getDataRange(ifstream& infile)
	{
		mDataMax = -FLT_MAX;
		mDataMin = FLT_MAX;
		mDataMean = 0.0;

		for(size_t i=0; i<(size_t)mHeaderMRC.ny * (size_t)mHeaderMRC.nx * (size_t)mHeaderMRC.nz;  i++)
		{
			voxelType value;
			infile.read(reinterpret_cast<char*>(&value), sizeof(voxelType));
			mDataMax=std::max((float)value,mDataMax);
			mDataMin=std::min((float)value,mDataMin);
			mDataMean+=(float)value;
		}

		mDataMean=mDataMean/((float)mHeaderMRC.ny * (float)mHeaderMRC.nx * (float)mHeaderMRC.nz);
	}

	vector<unsigned int> MRCStack::getResolution()
	{
	    return resolution;
	}


	float MRCStack::getStackMin()
	{
		return mDataMin;
	}

	float MRCStack::getStackMax()
	{
		return mDataMax;
	}

	float MRCStack::getInputMeanValue()
	{
		return mInputMeanValue;
	}

	float MRCStack::getStackMean()
	{
		return mDataMean;
	}

	unsigned int MRCStack::getNumberOfProjections()
	{
		return mHeaderMRC.nz;
	}

	MRCHeader MRCStack::getStackHeader()
	{
		return mHeaderMRC;
	}

	char* MRCStack::getExtraData()
	{
		return extraData;
	}

	size_t MRCStack::getProjectionSize()
	{
		return mProjectionSize;
	}

