#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <mpi.h>
#include "sa_Alignment.h"
#include "wedgeMaskGenerator.h"
#include "extractSubtomos.h"
#include "fftRoutines.h"
#include "emFile.h"
#include <algorithm>
#include "defines.h"
#include "exception.h"

using namespace std;

string dateAndTimeStamp(time_t startTime)
{
    time_t endTime;
    //struct tm * timeinfo;
    //char buffer[80];

    time(&endTime);

    double seconds_total = difftime(endTime, startTime);
    unsigned int hours = ((unsigned int)seconds_total) / 3600;
    unsigned int minutes = (((int)seconds_total) % 3600) / 60;
    unsigned int seconds = (((int)seconds_total) % 3600) % 60;
    //timeinfo = localtime(&rawtime);
    stringstream timeDifference;


    if (hours>0)
        timeDifference << hours << "h ";

    if (minutes > 0 || hours > 0)
        timeDifference << minutes << "m ";

    timeDifference << seconds << "s";

    //strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
    //std::string stamp(buffer);

    return timeDifference.str();
}

int main(int argc, const char* argv[]) {
    
    int initialized, finalized;

    MPI_Initialized(&initialized);
    if (!initialized)
        MPI_Init(NULL, NULL);


    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    time_t startTime;
    time(&startTime);

    vector<std::string> argList;

    for (int i = 1; i<argc; i++)
    {
        argList.push_back(argv[i]);
    }


    ParameterSetup params(argList);

    Log* logger = new Log(params.LogFileName(),world_rank,world_rank);
    logger->initLogFile();

    if(world_rank==0)
    {
        //initLogFile(params.LogFileName());
        params.printParamFile();
    }

    if (params.CreateWedgeMasks())
    {
        logger->printOut("Generating wedge masks...");
        WedgeMaskGenerator* wmGen = WedgeMaskGenerator::create(params, world_size, world_rank);
        wmGen->generateMask();
        delete wmGen;
        logger->printOut("...done");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm evenOddComm;

    if (params.RunSA())
    {
        sa_Alignment* newRec = sa_Alignment::create(params, *(logger), world_rank,world_size,evenOddComm);
        newRec->run();
        delete newRec;
    }

    logger->finishLogFile();

    free(evenOddComm);

    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();

    delete logger;

    return 0;
    /*
    if (params.Algorithm() == "3dctf")
    {
        cout << "Starting WBP reconstruction with 3D-CTF correction." << endl << endl;
        CTF3d* newRec = new CTF3d(params);
        newRec->run();
        cout << "Reconstruction finished successfully! \n Required time: " << dateAndTimeStamp(startTime) << endl;
        delete newRec;
    }
    else if (params.Algorithm() == "defocus")
    {
        cout << "Starting computation of defocus." << endl << endl;
        Defocus* newRec = new Defocus(params);
        newRec->run();
        cout << "Finished successfully! \n Required time: " << dateAndTimeStamp(startTime) << endl;
        delete newRec;
    }
    else if (params.Algorithm() == "filterProjections")
    {
        cout << "Starting filtering projections." << endl << endl;
        FilterProjections* newRec = new FilterProjections(params);
        newRec->run();
        cout << "Filtering of projections finished successfully! \n Required time: " << dateAndTimeStamp(startTime) << endl;
        delete newRec;
    }
    else if (params.Algorithm() == "ctfCorrection")
    {
        cout << "Starting CTF correction." << endl << endl;
        CTFCorrection* newRec = new CTFCorrection(params);
        newRec->run();
        cout << "CTF correction finished successfully! \n Required time: " << dateAndTimeStamp(startTime) << endl;
        delete newRec;
    }
    else
    {
        cout << "An algorithm was not specified! Use -Algorithm defocus/ctfCorrection/filterProjections/3dctf option to specify the algorithm you want to use." << endl << endl;
    }

    return 0;*/
}
