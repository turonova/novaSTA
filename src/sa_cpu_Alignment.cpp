#include <sstream> 
#include <string>
#include <float.h>
#include <string.h>
#include "sa_cpu_Alignment.h"
#include "maskRoutines.h"
#include "fftRoutines.h"
#include "motl.h"
#include "extractSubtomos.h"

using namespace vr;

sa_cpu_Alignment::sa_cpu_Alignment(ParameterSetup& iParams, Log& iLog, int iProcessID, int iNumberOfProcesses, MPI_Comm iActiveComm,int iGlobalProcessID)
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

sa_cpu_Alignment::~sa_cpu_Alignment()
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

void sa_cpu_Alignment::computePartialAverage(unsigned int iteration)
{
    unsigned int currentTomoWedge = 0;
    
    Particle* subtomo = new Particle(dim, params.InterpolationType(), params.InterpolationSigma());
    subtomo->initWithValue(0.0f);


    for (size_t m = 0; m < motl->Data().size(); m++)
    {
        // Check if motl should be procesed based on class
        if (motl->Data()[m].classNumber != params.SubtomoClass())
            continue;

        // Check wedge and generated mask if needed
        wedgeMask->updateMask(currentTomoWedge, motl->Data()[m].tomoNumber, false, false);

        // Read in subtomogram
        subtomo->getData(createName(subtomoName, motl->Data()[m].subtomoNumber, ".em"));

        vector<float> shift = { motl->Data()[m].xShift, motl->Data()[m].yShift, motl->Data()[m].zShift };

        Quaternion rotation(-motl->Data()[m].psi, -motl->Data()[m].phi, -motl->Data()[m].theta);
        updateResults(subtomo, -shift, rotation );

    }

    partialAverage = subtomo->Data(Volume::DataToProcess::rotatedD);
    partialWedgeMask = subtomo->Data(Volume::DataToProcess::rotatedMaskD);

    delete subtomo;
}

void sa_cpu_Alignment::prepareSubtomogram(Particle* originalSubtomo, Quaternion& rotation, vector<float>& shift)
{
    if (params.UseRosemanCC())
    {
        originalSubtomo->realForwardFFT(Volume::DataToProcess::originalD);
        originalSubtomo->maskFFT2(bandPassMaskReduced);
    }
    else 
    {
        originalSubtomo->maskAndNormalize(mask->rotateAndShift(rotation, shift), volSignal);
        originalSubtomo->realForwardFFT(Volume::DataToProcess::maskedD);
        originalSubtomo->maskFFT(bandPassMaskReduced);
    }

}

void sa_cpu_Alignment::align(unsigned int iteration)
{

    if(!runAlignment(iteration))
        return;

    unsigned int currentTomoWedge = 0;
    
    Particle* subtomo = new Particle(dim,params.InterpolationType(),params.InterpolationSigma());
    subtomo->initWithValue(0.0f);
    ref->initWithValue(0.0f);

	if(params.CcMaskIsSphere())
		ref->initWithCCMask(ccMask->Data());
	
	
    vector<float> ccf(volSize);
    vector<float> shift(3);

    unsigned int rotationType = getRotationType();

    for (size_t m = 0; m < motl->Data().size(); m++ )
    {
        // Check if motl should be procesed based on class
        if (motl->Data()[m].classNumber != params.SubtomoClass())
            continue;

        // Check wedge and generated mask if needed
        wedgeMask->updateMask(currentTomoWedge, motl->Data()[m].tomoNumber, !params.UseRosemanCC(), true, bandPassMaskReduced);
       
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

        prepareSubtomogram(subtomo, originalRotation, oShift);

        vector<float> newAngles;
        // Loop over all rotations
        for(size_t phiId = 0; phiId<phiAngles.size(); phiId++)
        {
            for(size_t rot=0; rot<rotations.size(); rot++)
            {

                Quaternion coneRotation=Quaternion::mult(rotations[rot],originalRotation);
                newAngles=coneRotation.getEulerAngles();
                Quaternion newRotation(motl->Data()[m].phi+phiAngles[phiId],newAngles[1],newAngles[2]);

                //stringstream printAngles;
                //printAngles << motl[m+16]+phiAngles[phiId] << "\t" << newAngles[1] << "\t" << newAngles[2];
               // printOut(printAngles,processID);

                ref->rotate(newRotation, mask->Data(), ccMask->Data(), rotationType);

                computeSubtomoRefCC(subtomo, ccf);

                vec_mult(ccf, ref->Data(Volume::DataToProcess::rotatedCCMaskD) );

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
            updateResults(subtomo, -tshift, finalRotation.inverseQuaternion());
        }
    
    }

    partialAverage = subtomo->Data(Volume::DataToProcess::rotatedD);
    partialWedgeMask = subtomo->Data(Volume::DataToProcess::rotatedMaskD);

    writeIntermediateOutputs(iteration);

    delete subtomo;
}

void sa_cpu_Alignment::computeRosemanCC(Particle* subtomo, vector<float>& ccf)
{
    // Calculate shifted fourier transform of rotated reference
    ref->realForwardFFT(Volume::DataToProcess::rotatedD);

    // Apply bandpass filtered wedge and edge mask
    ref->maskFFT2(wedgeMask->shiftedMask(true), false);

    FFTRoutines::inverseReal3DTransform(dimx, dimy, dimz, ref->DataFFT(), ccf);

    // sum rotated Mask
    double maskSum = vec_sum(ref->Data(Volume::DataToProcess::rotatedMaskD));

    // mask and normalize reference
    vector<float> maskedRef = ref->Data(Volume::DataToProcess::rotatedMaskD)*ccf;

    float maskedMean = vec_sum(maskedRef) / maskSum; // maskedRef.size();
    double sigmaRef = 0.0;

    for (size_t i = 0; i < maskedRef.size(); i++)
    {
        if (maskedRef[i] != 0.0f)
        {
            maskedRef[i] = maskedRef[i] - maskedMean;
            sigmaRef = sigmaRef + (double)(maskedRef[i] * maskedRef[i]);
        }
    }

    // Normalization factor of references
    sigmaRef = sqrt(sigmaRef);

    // Fourier transform of masked ref
    vector<complex<float>> maskedFFT(ref->DataFFT().size());
    FFTRoutines::forwardReal3DTransform(dimx, dimy, dimz, maskedRef, maskedFFT);

    // Convolution of masked reference and subtomo
    vec_mult(maskedFFT, subtomo->DataFFT());
    vector<float> numerator(maskedRef.size());

    FFTRoutines::inverseReal3DTransform(dimx, dimy, dimz, maskedFFT, numerator);

    // Calculate demoninator

    // Fourier transform of mask
    FFTRoutines::forwardReal3DTransform(dimx, dimy, dimz, ref->Data(Volume::DataToProcess::rotatedMaskD), maskedFFT);

    // Mean of subtomo under mask
    vector<float> meanSubtomo(maskedRef.size());
    vector<complex<float>> meanMasked = maskedFFT*subtomo->DataFFT();

    FFTRoutines::inverseReal3DTransform(dimx, dimy, dimz, meanMasked, meanSubtomo);

    // Mean intensity of subtomo under mask
    vec_mult(maskedFFT, subtomo->DataFFTSquaredConj());
    vector<float> intensitySubtomo(maskedRef.size());
    FFTRoutines::inverseReal3DTransform(dimx, dimy, dimz, maskedFFT, intensitySubtomo);

    // Calculate denominator (of eq 5 in paper)

    size_t shift = floor(dimx / 2.0f) + 1;

    for (size_t x = 0; x < dimx; x++)
    for (size_t y = 0; y < dimy; y++)
    for (size_t z = 0; z < dimz; z++)
    {
        size_t inputIndex = volSize - 1 - (x + y*dimx + z*dimx*dimy);
        size_t outputIndex = (x + shift) % dimx + ((y + shift) % dimy)*dimx + ((z + shift) % dimz)*dimx*dimy;
        double denominator = sqrt(intensitySubtomo[inputIndex] - meanSubtomo[inputIndex] * meanSubtomo[inputIndex] / maskSum)*sigmaRef;
        if (denominator != 0.0)
            ccf[outputIndex] = numerator[inputIndex] / denominator;
        else
            ccf[outputIndex] = -1.0f;
    }
}

void sa_cpu_Alignment::computeSubtomoRefCC(Particle* subtomo,vector<float>& ccf)
{
    if (params.UseRosemanCC())
    {
        computeRosemanCC(subtomo, ccf);
        return;
    }
    
    // Calculate shifted fourier transform of rotated reference
    ref->realForwardFFT(Volume::DataToProcess::rotatedD);

    // Apply bandpass filtered wedge and edge mask
    ref->maskFFT(wedgeMask->shiftedMask(true));

    // Calculate cross correlation and apply rotated ccmask
    vec_conj_mult(ref->DataFFT(), subtomo->DataFFT());

    FFTRoutines::inverseReal3DTransform(dimx, dimy, dimz, ref->DataFFT(), ccf);
    FFTRoutines::fftshift(ccf, dimx, dimy, dimz, 1, 0);

    vec_div(ccf, volSize);

}

void sa_cpu_Alignment::updateResults(Particle* subtomogram, vector<float> shifts, Quaternion quat)
{
    subtomogram->shiftAndRotate(shifts,quat,wedgeMask->Data(),binarizeWedgeMask);
    numberOfSubtomos += 1;
}

void sa_cpu_Alignment::divideByMeanAmplitude(vector<complex<double>>& data)
{
    double normFactor = 0.0;
    for (size_t i = 0; i< data.size(); i++)
    {
        normFactor = normFactor + (data[i].real()*data[i].real() + data[i].imag()*data[i].imag());
    }
    normFactor = volSize / sqrt(normFactor);
    vec_mult(data, normFactor); 
}
