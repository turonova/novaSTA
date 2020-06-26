#pragma once

#include <mpi.h>
#include <vector>
#include <cmath>

using namespace std;

class MPIRoutines
{
public:

    static int getCommSize(MPI_Comm activeComm = MPI_COMM_WORLD);
    static int getProcessID(MPI_Comm activeComm = MPI_COMM_WORLD);

    static void setBarrier(MPI_Comm activeComm = MPI_COMM_WORLD);

    static void reduceData(vector<float>& partialData, vector<float>& reducedData, size_t dataSize, int processID, int numberOfProcesses, MPI_Comm activeComm = MPI_COMM_WORLD);
    static void reduceData(vector<float>& partialData, vector<float>& reducedData, size_t dataSize, MPI_Comm activeComm = MPI_COMM_WORLD);

    static void gatherData(const vector<float>& partialData, vector<float>& combinedData,int processID, int numberOfProcesses, MPI_Comm activeComm = MPI_COMM_WORLD);
    static void gatherData(const vector<float>& partialData, vector<float>& combinedData, MPI_Comm activeComm = MPI_COMM_WORLD);

    static int reduceNumber(int partialNumber,int processID, int numberOfProcesses, MPI_Comm activeComm = MPI_COMM_WORLD);
    static int reduceNumber(int partialNumber, MPI_Comm activeComm = MPI_COMM_WORLD);
    static float reduceFloatNumber(float partialNumber, MPI_Comm activeComm = MPI_COMM_WORLD);

    static void bcastFloatNumber(float& number, MPI_Comm activeComm = MPI_COMM_WORLD);
	static void bcastBoolean(bool& value, MPI_Comm activeComm = MPI_COMM_WORLD);
	
    template <typename T>
    static void computeEvenDistribution(vector<T>& wholeData, typename vector<T>::iterator& itStart, typename vector<T>::iterator& itEnd, int numberOfParts, int partNumber);
};
