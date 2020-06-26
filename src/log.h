#pragma once

#include <iostream>
#include <fstream>

using namespace std;


class Log
{
public:

    Log(){};
    Log(string fileName);
    Log(string fileName, int iProcessID, int iGlobalProcessID);
    Log(const Log&){};

    Log operator=(const Log& iLog);

    ~Log();

    void initLogFile();
    void finishLogFile();

    int getProcessID() const;
    int getGlobalProcessID() const;
    string getLogFileName() const;

    void setProcessID(int iProcessID);
    void setGlobalProcessID(int iGlobalProcessID);

    void printOut(string message, bool useProcessID = false);
    void printOutTime(string message, bool useProcessID = false);

    void openLogFile();
    void closeLogFile();

private:

    //Log(const Log&);

    void printDivider(bool useProcessID = false);

    int processID;
    int globalProcessID;

    ofstream logFile;
    string logFileName;
};
