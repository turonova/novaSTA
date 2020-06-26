#pragma once

#include <vector>
#include <string>
#include <sstream>

using namespace std;

template <typename paramType>
struct inputParameter
{
    paramType value;
    bool set;
    string name;

    inputParameter(string iName)
    {
        name = iName;
        set = false;
    }
};

class ParameterSetup
{

public:
    ParameterSetup(std::vector<string> argList);
    ParameterSetup();

    ~ParameterSetup();

    void printParamFile();
    void applyEvenOddSplit(string nameUpdate);

    bool UseGPU();
    bool CreateReference();
    bool ExtractSubtomos();
    unsigned int Symmetry();
    int SubtomoClass();
    float Threshold();
    float LowPassSigma();
    float HighPassSigma();

    // parameters without default values
    string MotlName();
    string FolderName();
    string RefName();
    string MaskName();
    string CcMaskName();
    string WedgeMaskName();
    string WedgeDataName();
    string SubtomoName();
    string TomoName();
    string FscMask();

    int NumberOfIterations();
    int StartingIndex();

    vector<float> ConeAngle();
    vector<float> ConeSampling();

    vector<float> InplaneAngle();
    vector<float> InplaneSampling();

    vector<float> LowPass();
    vector<float> HighPass();

    size_t VolumeSize();

    unsigned int TomoDigits();

    int Rerun();
    
    bool UseRosemanCC();

    bool CopySubtomograms();
    string TempSubtomoFolder();
	
    unsigned int RealStartIndex();
    unsigned int RotationAndFilterIndex();

    string MotlBaseName();
    string WedgeBaseName();
    string RefBaseName();

    string RefOriginalName();
    string MotlOriginalName();

    bool ChangingBandPass();
    bool ChangingRotations();

    string CopySubtomoName();

    bool SplitIntoEvenOdd();
    string SplitType();

    float PixelSize();

    string LogFileName();

    float MotlBinFactor();

    float DistanceThreshold();
    bool CleanByDistance();
    bool CleanOutOfBoundsParticles();

    bool RenumberParticles();
    bool CleanByMeanGreyValue();

    unsigned int WedgeType();
    bool CreateWedgeMasks();

    bool RunSA();

    unsigned int InterpolationType();

    float InterpolationSigma();

	bool ComputeFSC();
	
private:

    void initVariables();

    void parseComFile(const string tiltcom);
    void parseCommandLine(vector<string> argList);

    void storeValues(string paramName, string paramValue, char separator);
    bool parameterWithoutValue(string paramName);

    void printParameter(string paramName, string paramValue);

    void checkParamSetup();
    void wrongParameterCombination();
    void parameterNotSetError(string paramName);

    string getFilename(string path);
    string createBaseName(string folder, string filename);

    void setInternalParameters();

    // parameters with default values
    bool useGPU;
    bool createReference;
    bool extractSubtomos;

    unsigned int symmetry;
    int subtomoClass;
    float threshold;
    float lowPassSigma; 
    float highPassSigma;
    int rerun;
    bool useRosemanCC;
    bool copySubtomograms;
    string allParameters;
    string copySubtomoName;

    bool splitIntoEvenOdd;
    string splitType;

    // Internal parameters - set based on user's input

    unsigned int realStartIndex;
    unsigned int rotationAndFilterIndex;
    bool changingBandPass;
    bool changingRotations;

    float motlBinFactor;

    string motlBaseName;
    string wedgeBaseName;
    string refBaseName;

    string refOriginalName;
    string motlOriginalName;

    string logFileName;

    bool cleanByDistance;
    bool cleanOutOfBoundsParticles;
    bool renumberParticles;
    bool cleanByMeanGreyValue;

    unsigned int wedgeType;
    bool createWedgeMasks;

    bool runSA;

    unsigned int interpolationType;
    float interpolationSigma;

	bool computeFSC;
	
    // parameters without default values
    inputParameter<string> motlName = inputParameter<string>("motl");
    inputParameter<string> folderName = inputParameter<string>("folder");
    inputParameter<string> refName = inputParameter<string>("ref");
    inputParameter<string> maskName = inputParameter<string>("mask");
    inputParameter<string> ccMaskName = inputParameter<string>("ccMask");
    inputParameter<string> wedgeMaskName = inputParameter<string>("wedgeList");
    inputParameter<string> wedgeDataName = inputParameter<string>("wedgeInputSpec");
    inputParameter<string> subtomoName = inputParameter<string>("subtomograms");
    inputParameter<string> tomoName = inputParameter<string>("tomograms");
    inputParameter<string> fscMask = inputParameter<string>("fscMask");

    inputParameter<int> numberOfIterations = inputParameter<int>("iter");
    inputParameter<int> startingIndex = inputParameter<int>("startIndex");

    inputParameter<vector<float>> coneAngle = inputParameter<vector<float>>("coneAngle");
    inputParameter<vector<float>> coneSampling = inputParameter<vector<float>>("coneSampling");

    inputParameter<vector<float>> inplaneAngle = inputParameter<vector<float>>("inplaneAngle");
    inputParameter<vector<float>> inplaneSampling = inputParameter<vector<float>>("inplaneSampling");

    inputParameter<vector<float>> lowPass = inputParameter<vector<float>>("lowPass");
    inputParameter<vector<float>> highPass = inputParameter<vector<float>>("highPass");

    inputParameter<size_t> volumeSize = inputParameter<size_t>("subtomoSize");
    inputParameter<unsigned int> tomoDigits = inputParameter<unsigned int>("tomoDigits");
    
    inputParameter<string> tempSubtomoFolder = inputParameter<string>("tempSubtomoFolder");
    inputParameter<float> pixelSize = inputParameter<float>("pixelSize");

    inputParameter<float> distanceThreshold = inputParameter<float>("distanceThreshold");
};
