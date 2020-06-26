#include "sa_gpu_Alignment.h"
#include <sstream> 
#include <string>
#include <float.h>
#include <string.h>
#include "maskRoutines.h"
#include "fftRoutines.h"
#include "motl.h"
#include "extractSubtomos.h"
#include "cuComplex.h"

using namespace vr;


sa_gpu_Alignment::sa_gpu_Alignment(ParameterSetup& iParams, Log& iLog, int iProcessID, int iNumberOfProcesses, MPI_Comm iActiveComm, int iGlobalProcessID)
{
    params = iParams;
    processID = iProcessID;
    numberOfProcesses = iNumberOfProcesses;
    activeComm = iActiveComm;
    globalProcessID=iGlobalProcessID;
    
    logger = iLog;
    logger.setGlobalProcessID(globalProcessID);
    logger.setProcessID(processID);
};

sa_gpu_Alignment::~sa_gpu_Alignment()
{
    delete wedgeMask;

    if (params.NumberOfIterations() > 0)
    {
        delete ref;
        delete mask;
        delete ccMask;
    }
    
    if(params.SplitIntoEvenOdd())
        delete fsc;
        
    logger.closeLogFile();
}

void sa_gpu_Alignment::computePartialAverage(unsigned int iteration)
{
    unsigned int currentTomoWedge = 0;

    vector<float> subtomo(volSize);

    gpu = new GPURoutines(dimx, dimy, dimz, false, processID);

    gpu->allocateDeviceMemory();
    gpu->prepareSubtomoTexture("linear");

    for (size_t m = 0; m < motl->Data().size(); m++)
    {
        // Check if motl should be procesed based on class
        if (motl->Data()[m].classNumber != params.SubtomoClass())
            continue;

        // Check wedge and generated mask if needed
        // do not mask edges regardless of roseman !!!
        bool mask_updated=wedgeMask->updateMask(currentTomoWedge,motl->Data()[m].tomoNumber,false,false);
        
        if (mask_updated)
            gpu->copyAndBindWedgeTexture(wedgeMask->Data());
    

        // Create a quaternion based on motl
        Quaternion rotation(-motl->Data()[m].psi, -motl->Data()[m].phi, -motl->Data()[m].theta);

        // Read in subtomogram
        emFile::read(createName(subtomoName, motl->Data()[m].subtomoNumber, ".em"), subtomo);
        
        gpu->shiftSubtomogram(subtomo, -motl->Data()[m].xShift, -motl->Data()[m].yShift, -motl->Data()[m].zShift);
        gpu->rotateSubtomogramAndWedge(subtomo, wedgeMask->Data(), rotation, binarizeWedgeMask);

        numberOfSubtomos += 1;
    }

    gpu->getRotatedSubtomoAndWedge(partialAverage, partialWedgeMask);

    delete gpu;
}

void sa_gpu_Alignment::prepareSubtomogram(Particle* originalSubtomo, vector<float>& rotatedMask)
{

}

//void sa_gpu_Alignment::updateResults()
//{
//    partialAverage = partialAverage + subtomogram->shiftAndRotate(shifts,quat);
//    partialWedgeMask = partialWedgeMask + wedgeMask->rotateWedgeMask(quat);
//    numberOfSubtomos += 1;
//}



void sa_gpu_Alignment::align(unsigned int iteration)
{

    if(!runAlignment(iteration))
        return;

    unsigned int currentTomoWedge = 0;
    
    Particle* subtomo = new Particle(dim);

    vector<float> ccf(volSize);
    vector<float> shift(3);

    gpu = new GPURoutines(dimx, dimy, dimz, params.UseRosemanCC(),processID);

    gpu->allocateDeviceMemory();
    gpu->prepareSubtomoTexture("linear");
    gpu->prepareReferenceTexture(ref->Data(), mask->Data(), ccMask->Data(), "linear");

    for (size_t m = 0; m < motl->Data().size(); m++)
    {
        // Check if motl should be procesed based on class
        if (motl->Data()[m].classNumber != params.SubtomoClass())
            continue;

        // Check wedge and generated mask if needed
        bool mask_updated=wedgeMask->updateMask(currentTomoWedge, motl->Data()[m].tomoNumber, !params.UseRosemanCC(), true, bandPassMaskReduced);        
        
        if (mask_updated)
            gpu->copyAndBindWedgeTexture(wedgeMask->Data());

        //Parse Euler angles from motl
        Quaternion originalRotation(motl->Data()[m].phi,motl->Data()[m].psi,motl->Data()[m].theta);

        // Parse shifts from motl file
        vector<float> oShift = { motl->Data()[m].xShift, motl->Data()[m].yShift, motl->Data()[m].zShift };

        // Initialize CCC
        float ccc = -1.0f;
        
        // Initialize shift vector and new angles
        Quaternion finalRotation;
        vector<float> tshift = { 0.0f, 0.0f, 0.0f };

        // Read in subtomogram
        subtomo->getData(createName(subtomoName, motl->Data()[m].subtomoNumber, ".em"));
        
        //Mask subtomogram
        if (params.UseRosemanCC())
        {
            gpu->maskSubtomogram(subtomo->Data(), bandPassMaskReduced);
        }
        else
        {
            subtomo->maskAndNormalize(mask->rotateAndShift(originalRotation, oShift), volSignal);
            gpu->maskSubtomogram(subtomo->DataMasked(), bandPassMaskReduced);
        }

        vector<float> newAngles;
        // Loop over all rotations
        for(size_t phiId = 0; phiId<phiAngles.size(); phiId++)
        {
            for(size_t rot=0; rot<rotations.size(); rot++)
            {
                Quaternion coneRotation=Quaternion::mult(rotations[rot],originalRotation);
                newAngles=coneRotation.getEulerAngles();
                Quaternion newRotation(motl->Data()[m].phi+phiAngles[phiId],newAngles[1],newAngles[2]);
                
                //rotateVolume(ref,rotated_ref,newRotation);
                gpu->rotateTexture(newRotation);

                gpu->computeCC(ccf, wedgeMask->shiftedMask(true));

                float cccTemp = findSubpixelPeak(ccf,shift);


                // If CCC is greater than current ccc, update
                // ccc, Euler angles, and shifts
                if (cccTemp > ccc)
                {
                    ccc = cccTemp;
                    finalRotation=newRotation;

                    if (shift[0] < 0.0f && shift[1] < 0.0f && shift[2] < 0.0f)  // If the position of the peak is less than zero
                    {
                        tshift = oShift; // Set shift to old shifts
                    }
                    else
                    {
                        tshift = shift - volCenter; // Else set shift to peak position minus center position
                    }
                }
            }
        }// rotations

        newAngles= finalRotation.getEulerAngles();
        
        updateConvergenceStats(oShift,tshift,originalRotation, finalRotation, motl->Data()[m].phi, newAngles[0]);
        
        
        motl->setAngles(m,newAngles[0],newAngles[1],newAngles[2]);
        motl->setShifts(m,tshift[0],tshift[1],tshift[2]);
        motl->setMetric(m,ccc);

        //kick off bad particles
        if( ccc <= params.Threshold() )
            motl->setClass(m,-1.0f);   //bad CCF : kick into class -1
        else
        {
            motl->setClass(m,1.0f);    //good particles -> class one
            
            gpu->shiftSubtomogram(subtomo->Data(), -tshift[0], -tshift[1],-tshift[2]);
            Quaternion iq = finalRotation.inverseQuaternion();
            gpu->rotateSubtomogramAndWedge(subtomo->Data(), wedgeMask->Data(), iq,binarizeWedgeMask);
            numberOfSubtomos += 1;
        }
    
    }

    gpu->getRotatedSubtomoAndWedge(partialAverage, partialWedgeMask);

    delete gpu;
    delete subtomo;
}

