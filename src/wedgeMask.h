#pragma once

#include <vector>
#include <map>
#include "volume.h"
#include "parameterSetup.h"

using namespace std;

enum WedgeType
{
    bowtie,         // the original one specified as min and max angle in wedge list 3xNumberOfTomos (file name path/wedge_list_name.em)
    star,          // precise "star" mask created by a DDA - specified in wedge list NumberOfAnglesxNumberOfTomos - less tilts should be filled in with zeros (file name path/wedge_list_name.em)
    inputMasks,     // input wedge mask - different masks for different tomograms (file name path/wedge_mask_name_#tomo.em specified as path/wedge_mask_name)
    singleMask     // input wedge mask - one mask for all tomogramgs (file name path/wedge_mask_name.em)
};

class WedgeMask: public Volume
{
public:
    
    WedgeMask(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC);

    static WedgeMask* create(ParameterSetup& iParams, size_t aSize, vector<size_t>& aDim);
    virtual bool updateMask(unsigned int& currentWedge, unsigned int newWedge, bool maskEdges, bool shiftMask, const vector<float>& wedgeFilter = vector<float>()) = 0;

    vector<float>& rotateWedgeMask(Quaternion q);
    //void createMask(unsigned int id, bool shift = false);
    
    vector<float>& shiftedMask(bool reducedSize = false);
    void filterWithMask(vector<float>& edgeMask);

    
    void filter(vector<float> filterMask);

 
protected:

    void binarizeRotated();
    void shiftAndReduce();
    void loadMask(string maskName);

    size_t size;

    vector<float> shiftedReduced;
    bool useRosemanCC;
    string wedgeBaseName;

    ParameterSetup params;
};



class wm_Bowtie : public WedgeMask
{
public:    
    wm_Bowtie(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC);
    ~wm_Bowtie(){};

    bool updateMask(unsigned int& currentWedge, unsigned int newWedge, bool maskEdges, bool shiftMask, const vector<float>& wedgeFilter = vector<float>());

    void createMask(unsigned int id, bool shift);

    vector<float>& rotateWedgeMask(Quaternion q);
private:

    void initWedgeList(string filename);
    void generateMask(unsigned int id, bool maskEdges, bool shift = false);

    map<unsigned int, pair<float, float>> wedgeList;
};


class wm_Star : public WedgeMask
{
public:
    wm_Star(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC);
    ~wm_Star(){};

    bool updateMask(unsigned int& currentWedge, unsigned int newWedge, bool maskEdges, bool shiftMask, const vector<float>& wedgeFilter = vector<float>());

private:

    void generateMask(unsigned int id, bool maskEdges, bool shift = false);

};

class wm_SingleMask : public WedgeMask
{
public:
    wm_SingleMask(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC);
    ~wm_SingleMask(){};

    bool updateMask(unsigned int& currentWedge, unsigned int newWedge, bool maskEdges, bool shiftMask, const vector<float>& wedgeFilter = vector<float>());

private:
    bool filtered;

};

class wm_InputMasks : public WedgeMask
{
public:
    wm_InputMasks(size_t aSize, vector<size_t>& aDim, string aWedgeName, bool aUseRosemanCC, int tomoDigits);
    ~wm_InputMasks(){};

    bool updateMask(unsigned int& currentWedge, unsigned int newWedge, bool maskEdges, bool shiftMask, const vector<float>& wedgeFilter = vector<float>());
    
private:

    int tomoDigits;

};