#include "log.h"
#include "exception.h"

using namespace std;


Log::Log(string fileName)
{
    logFileName = fileName;
    processID = 0;
    globalProcessID = 0;
}

Log::Log(string fileName, int iProcessID, int iGlobalProcessID)
{
    logFileName = fileName;
    processID = iProcessID;
    globalProcessID = iGlobalProcessID;
}


Log Log::operator=(const Log& iLog)
{
    processID=iLog.getProcessID();
    globalProcessID=iLog.getGlobalProcessID();
    logFileName = iLog.getLogFileName();

    return *this;
}


Log::~Log()
{
    if(logFile.is_open())
        logFile.close();
}

void Log::initLogFile()
{
    if(globalProcessID!=0)
        return;

    try
    {

        logFile.open(logFileName.c_str());


        if (!logFile.good())
        {
            throw ExceptionFileOpen(logFileName);
        }

        printDivider();
        printOut("#  TomSA v0.9");
        printDivider();
        logFile << endl;

        time_t rawTime;
        time(&rawTime);


        printOutTime("Processing started on ");

        logFile.close();
    }
    catch (ExceptionFileOpen& e)
    {
        cout << e.what() << endl;
    }
}

void Log::finishLogFile()
{
    if(globalProcessID!=0)
        return;

    try
    {
        logFile.open(logFileName, std::ofstream::out | std::ofstream::app);


        if (!logFile.good())
        {
            throw ExceptionFileOpen(logFileName);
        }

       // printDivider();
        printOutTime("\n\n  Processing ended on ");
      //  printDivider();

        logFile.close();
    }
    catch (ExceptionFileOpen& e)
    {
        cout << e.what() << endl;
    }
}

void Log::openLogFile()
{
    if(processID!=0)
        return;

    try
    {
        logFile.open(logFileName, std::ofstream::out | std::ofstream::app);

        if (!logFile.good())
        {
            throw ExceptionFileOpen(logFileName);
        }
    }
    catch (ExceptionFileOpen& e)
    {
        cout << e.what() << endl;
    }

}

void Log::closeLogFile()
{
    if(logFile.is_open())
        logFile.close();
}

void Log::printOut(string message,  bool useProcessID)
{
    if(!useProcessID && globalProcessID!=0)
        return;
    else if (useProcessID && processID!=0)
        return;

    logFile << message << endl;
}


void Log::printOutTime(string message, bool useProcessID)
{
    if(!useProcessID && globalProcessID!=0)
        return;
    else if (useProcessID && processID!=0)
        return;

    time_t rawTime;
    time(&rawTime);

    logFile << message << ctime(&rawTime) << endl;
}


int Log::getProcessID() const
{
    return processID;
}

int Log::getGlobalProcessID() const
{
    return globalProcessID;
}

string Log::getLogFileName() const
{
    return logFileName;
}

void Log::setProcessID(int iProcessID)
{
    processID = iProcessID;
}

void Log::setGlobalProcessID(int iGlobalProcessID)
{
    globalProcessID = iGlobalProcessID;
}

void Log::printDivider(bool useProcessID)
{
    if(!useProcessID && globalProcessID!=0)
        return;
    else if (useProcessID && processID!=0)
        return;

    logFile << "######################################################" << endl;
}

