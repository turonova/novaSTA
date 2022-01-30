#include "mpiRoutines.h"


template void MPIRoutines::computeEvenDistribution<unsigned int>(vector<unsigned int>& wholeData, vector<unsigned int>::iterator& itStart, vector<unsigned int>::iterator& itEnd, int numberOfParts, int partNumber);


void MPIRoutines::bcastFloatNumber(float& number, MPI_Comm activeComm)
{
	MPI_Bcast(&number,1,MPI_FLOAT,0,activeComm);
}

void MPIRoutines::bcastBoolean(bool& value, MPI_Comm activeComm)
{
	MPI_Bcast(&value,1,MPI_C_BOOL,0,activeComm);
}

int MPIRoutines::getCommSize(MPI_Comm activeComm)
{
    int numberOfProcesses;
    MPI_Comm_size(activeComm, &numberOfProcesses);
    return numberOfProcesses;
}

int MPIRoutines::getProcessID(MPI_Comm activeComm)
{
    int processID;
    MPI_Comm_rank(activeComm, &processID);
    return processID;
}

void MPIRoutines::setBarrier(MPI_Comm activeComm)
{
    MPI_Barrier(activeComm);
}

void MPIRoutines::gatherData(const vector<float>& partialData, vector<float>& combinedData, MPI_Comm activeComm)
{
    int processID;
    int numberOfProcesses;

    MPI_Comm_rank(activeComm, &processID);
    MPI_Comm_size(activeComm, &numberOfProcesses);

    gatherData(partialData,combinedData,processID,numberOfProcesses,activeComm);
}

void MPIRoutines::reduceData(vector<float>& partialData, vector<float>& reducedData, size_t dataSize, MPI_Comm activeComm)
{
    int processID;
    int numberOfProcesses;

    MPI_Comm_rank(activeComm, &processID);
    MPI_Comm_size(activeComm, &numberOfProcesses);

    reduceData(partialData,reducedData,dataSize,processID,numberOfProcesses,activeComm);

}

float MPIRoutines::reduceFloatNumber(float partialNumber, MPI_Comm activeComm)
{
	int processID;
    int numberOfProcesses;

    MPI_Comm_rank(activeComm, &processID);
    MPI_Comm_size(activeComm, &numberOfProcesses);
	
	float reducedNumber;
    MPI_Reduce(&partialNumber, &reducedNumber, 1 , MPI_FLOAT, MPI_SUM, 0, activeComm);

    return reducedNumber;
}



size_t MPIRoutines::reduceSizeTNumber(size_t partialNumber, MPI_Comm activeComm)
{
    int processID;
    int numberOfProcesses;

    MPI_Comm_rank(activeComm, &processID);
    MPI_Comm_size(activeComm, &numberOfProcesses);

    return reduceSizeTNumber(partialNumber, processID, numberOfProcesses, activeComm);
}

size_t MPIRoutines::reduceSizeTNumber(size_t partialNumber, int processID, int numberOfProcesses, MPI_Comm activeComm)
{
    size_t reducedNumber;
    MPI_Reduce(&partialNumber, &reducedNumber, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, activeComm);

    return reducedNumber;
}

int MPIRoutines::reduceNumber(int partialNumber, MPI_Comm activeComm)
{
    int processID;
    int numberOfProcesses;

    MPI_Comm_rank(activeComm, &processID);
    MPI_Comm_size(activeComm, &numberOfProcesses);

    return reduceNumber(partialNumber,processID,numberOfProcesses, activeComm);
}

int MPIRoutines::reduceNumber(int partialNumber,int processID, int numberOfProcesses, MPI_Comm activeComm)
{
    int reducedNumber;
    MPI_Reduce(&partialNumber, &reducedNumber, 1 , MPI_INT, MPI_SUM, 0, activeComm);

    return reducedNumber;
}

void MPIRoutines::gatherData(const vector<float>& partialData, vector<float>& combinedData,int processID, int numberOfProcesses, MPI_Comm activeComm)
{
    vector<int> partialSizes;
    vector<int> displacements;

    if (processID== 0)
    {
        partialSizes.resize(numberOfProcesses);
        displacements.resize(numberOfProcesses);
    }

    int partialSize =(int)partialData.size();
    MPI_Gather(&partialSize, 1, MPI_INT, partialSizes.data(), 1, MPI_INT, 0, activeComm);

    if(processID==0)
    {
       displacements[0] = 0;
       for(int i=1; i< numberOfProcesses;i++)
           displacements[i] = displacements[i-1] + partialSizes[i-1];
    }

    int combinedSize = MPIRoutines::reduceNumber(partialSize,processID,numberOfProcesses,activeComm);

    if (processID== 0)
    {
       combinedData.resize(combinedSize);
    }

    MPI_Gatherv(partialData.data(), partialSize, MPI_FLOAT, combinedData.data(), partialSizes.data(),displacements.data(), MPI_FLOAT, 0, activeComm);

}

void MPIRoutines::reduceData(vector<float>& partialData, vector<float>& reducedData, size_t dataSize,int processID, int numberOfProcesses, MPI_Comm activeComm)
{
    if (processID== 0)
    {
        reducedData.resize(dataSize);
    }

    MPI_Reduce(partialData.data(), reducedData.data(), dataSize, MPI_FLOAT, MPI_SUM, 0, activeComm);
}

template <typename T>
void MPIRoutines::computeEvenDistribution(vector<T>& wholeData, typename vector<T>::iterator& itStart, typename vector<T>::iterator& itEnd, int numberOfParts, int partNumber)
{

    size_t particlesPerPart = (size_t)std::round((float)wholeData.size() / ((float)numberOfParts));

    float particleRest = (float)wholeData.size() - (float)particlesPerPart*numberOfParts;
    float particlesUnchanged = numberOfParts - fabs(particleRest);

    int particleChange;

    if (particleRest<0)
        particleChange = -1;
    else
        particleChange = 1;

    if (partNumber<particlesUnchanged)
    {
        itStart = wholeData.begin() + particlesPerPart*partNumber;
        itEnd = itStart + particlesPerPart;
    }
    else
    {
        itStart = wholeData.begin() + particlesPerPart*particlesUnchanged + (particlesPerPart + particleChange)*(partNumber - particlesUnchanged);

        if (partNumber == (numberOfParts - 1))
            itEnd = wholeData.end();
        else
            itEnd = itStart + particlesPerPart + particleChange;

    }
}
