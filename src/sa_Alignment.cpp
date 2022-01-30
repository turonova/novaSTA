#include "sa_Alignment.h"
#include "sa_cpu_Alignment.h"
#include "sa_gpu_Alignment.h"

#include <sstream>
#include <string>
#include <float.h>
#include <string.h>
#include "sa_cpu_Alignment.h"
#include "maskRoutines.h"
#include "fftRoutines.h"
#include "motl.h"
#include "extractSubtomos.h"
#include "exception.h"
#include "mpiRoutines.h"

using namespace vr;

sa_Alignment* sa_Alignment::create(ParameterSetup& iParams, Log& iLog, int iProcessID, int iNumberOfProcesses, MPI_Comm iActiveComm)
{
    int commProcessID;
    int commNumberOfProcesses;

    
    if(iParams.SplitIntoEvenOdd())
    {
        if(iProcessID==0)
            splitIntoEvenOdd(iParams);

        MPI_Barrier(MPI_COMM_WORLD);

        int color;

        if(iProcessID<iNumberOfProcesses/2)
        {
            color=0;
            iParams.applyEvenOddSplit("even");
        }
        else
        {
            color=1;
            iParams.applyEvenOddSplit("odd");
        }

        MPI_Comm_split(MPI_COMM_WORLD, color, iProcessID, &iActiveComm);

        MPI_Barrier(iActiveComm);

        MPI_Comm_rank(iActiveComm, &commProcessID);
        MPI_Comm_size(iActiveComm, &commNumberOfProcesses);
    }
    else
    {
        commProcessID=iProcessID;
        commNumberOfProcesses = iNumberOfProcesses;
        iActiveComm=MPI_COMM_WORLD;
    }

    if (iParams.UseGPU())
        return new sa_gpu_Alignment(iParams, iLog, commProcessID, commNumberOfProcesses, iActiveComm,iProcessID);
    else
        return new sa_cpu_Alignment(iParams, iLog, commProcessID, commNumberOfProcesses, iActiveComm,iProcessID);
}

void sa_Alignment::splitIntoEvenOdd(ParameterSetup& iParams)
{
    // Reruning -> separate motls should exist already
    if(iParams.RealStartIndex()!=iParams.StartingIndex())
        return;

    string motlName = createName(iParams.MotlName(), iParams.StartingIndex(), ".em");

    Motl tempMotl(motlName);

    tempMotl.renumberParticles(iParams.RenumberParticles());

    vector<float> evenMotl;
    vector<float> oddMotl;

    tempMotl.splitIntoEvenOdd(evenMotl,oddMotl);

    emFile::write(createName(iParams.MotlBaseName() +"_even", iParams.StartingIndex(), ".em"),evenMotl,20,evenMotl.size()/20,1);
    emFile::write(createName(iParams.MotlBaseName()+"_odd", iParams.StartingIndex(), ".em"),oddMotl,20,oddMotl.size()/20,1);

}

void sa_Alignment::run()
{
    logger.openLogFile();

    initVaribales();
    prepareMotl();
    prepareMasks();
    extractSubtomos();

    MPI_Barrier(MPI_COMM_WORLD);
    preprocessMotl();

    MPI_Barrier(MPI_COMM_WORLD);
    
    copySubtomos();
    
    MPI_Barrier(MPI_COMM_WORLD);

    createStartingReference();

    MPI_Barrier(activeComm);

    unsigned int i;

    for(i=realStartinIndex; i<(params.StartingIndex()+params.NumberOfIterations()); i++)
    {
        MPIRoutines::bcastBoolean(terminateEarly);   
		
		if(terminateEarly)
		{
			logger.printOut("\n\n  FSC curve did not cross 0.143 threshold.");
			logger.printOut("  The processing will now terminate.");
			break;
		}
		
		logger.printOut("\n\n Starting iteration #" + to_string(i));

        numberOfSubtomos = 0;
        setBandPassFilter();
        setRotations();
        prepareReference(i);

        align(i+1);
        computeFinalAverage(i+1);

        MPI_Barrier(MPI_COMM_WORLD);

        combineEvenOdd(i+1);

        MPI_Barrier(MPI_COMM_WORLD);
        computeFinalStats(i+1);

        clean(i+1);

        rotationAndFilterIndex++;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    postprocessMotl(i);

};

float sa_Alignment::computeAddaptiveFilter()
{
    MPIRoutines::bcastFloatNumber(currentFscValue);    
		
	float lowPass;
	if(currentFscValue>EPS)
	{
		lowPass = ceil(params.VolumeSize()*params.PixelSize()/currentFscValue) -  fabs(params.LowPass()[rotationAndFilterIndex]);
	}
	else
	{
		logger.printOut("FSC curve did not cross 0.143 threshold and thus addaptive filter cannot be computed. The processing will terminate now.");
		exit(EXIT_SUCCESS);
	}
	
    return lowPass;
}

void sa_Alignment::computeFSC(vector<float>& evenVolume, vector<float>& oddVolume,unsigned int iteration)
{
    if(!params.ComputeFSC())
        return;
	
	fsc->computeFSC(evenVolume,oddVolume);
    fsc->printFSCCurve(createName(params.RefOriginalName(), iteration,"_fsc.txt"));

    currentFscValue = fsc->getFSCValue(0.143f);

    logger.printOut("\nFSC at 0.5 is: " + to_string(fsc->getFSCValue(0.5f)));
    logger.printOut("FSC at 0.143 is: " + to_string(fsc->getFSCValue(0.143f)));
	
	if(currentFscValue<EPS)
	{
		terminateEarly=true;
	}
}



void sa_Alignment::postprocessMotl(unsigned int iteration)
{
    if (!params.PostprocessMotl())
        return;

    string motlName=createName(params.MotlOriginalName(), iteration, ".em");
    string origMotlName = createName(params.MotlOriginalName(), iteration, "_before_post_processing.em");

    Motl* ppMotl = new Motl(motlName, 6);

    ppMotl->writeFullData(origMotlName);
    
    size_t nPaticlesBeforePostprocess = ppMotl->getNumberOfAllParticles();
    size_t nParticlesAfterPostprocess;

    if (params.ComputeAngularDistance() > 0)
    {
        StructureGeometry* geom = StructureGeometry::create(params);
        ppMotl->computeAngularDistance(*geom);
        logger.printOut(geom->getDescription());
        logger.printOut("Number of particles is " + to_string(ppMotl->getNumberOfAllParticles()));
        delete geom;
    }

    if(params.CleanByDistance())
    {
        ppMotl->cleanByDistance(params.DistanceThreshold());

        logger.printOut("\n\nParticles in motl were cleaned by distance threshold of " + to_string(params.DistanceThreshold()) + " voxels.");
        logger.printOut("Number of particles after the distance cleaning is " + to_string(ppMotl->getNumberOfAllParticles()));
    }

    if (params.UnifyClassNumber() > 0)
    {
        ppMotl->unifyClassNumber(params.UnifyClassNumber());
        logger.printOut("The class of all particles in the motl was set to " + to_string(params.UnifyClassNumber()));
    }

    ppMotl->writeFullData(motlName);

    logger.printOut("The motl before post-processing was saved as " + origMotlName);
    logger.printOut("The final motl was saved as " + motlName);

    delete ppMotl;
}

void sa_Alignment::copySubtomos()
{
    if(!params.CopySubtomograms())
    {
        subtomoName = params.SubtomoName();
        return;
    }
    
   // string subtomo = getFilename(params.SubtomoName());
    subtomoName = params.CopySubtomoName();
    
    
    logger.printOut("Copying subtomograms from " +  params.SubtomoName() + " to " + subtomoName);

    
    for (size_t m = 0; m < motl->Data().size(); m++)
    {
        string sourceName = createName(params.SubtomoName(), motl->Data()[m].subtomoNumber, ".em");
        string destName = createName(subtomoName, motl->Data()[m].subtomoNumber, ".em");
        ifstream  src(sourceName, std::ios::binary);
        ofstream  dst(destName,   std::ios::binary);

        dst << src.rdbuf();
        src.close();
        dst.close();   
    }

}

void sa_Alignment::updateConvergenceStats(vector<float>& oldShift, vector<float>& newShift, Quaternion& oldRotation, Quaternion& newRotation, float oldPhi, float newPhi)
{
    for(unsigned int i=0; i<3; i++)
    {
        convergenceStats[i]+=(newShift[i]-oldShift[i])*(newShift[i]-oldShift[i]);
    }

    float angle, distance;

    Quaternion::getAngularDistance(newRotation, oldRotation, angle, distance);

    convergenceStats[3] += angle*angle;
    convergenceStats[4] += distance*distance;
    convergenceStats[5] += (newPhi-oldPhi)*(newPhi-oldPhi);
}

unsigned int sa_Alignment::getRotationType()
{
    unsigned int rotationType;
    
    if (params.UseRosemanCC() && params.CcMaskIsSphere())
        rotationType = 6;
    else  if (params.UseRosemanCC())
        rotationType = 3;
    else if (params.CcMaskIsSphere())
        rotationType = 5;
    else
        rotationType=2;
   
    return rotationType;
}

void sa_Alignment::initVaribales()
{
    terminateEarly = false;
	currentFscValue = 0.0f;
	convergenceStats.resize(6,0.0f);

    dimx = params.VolumeSize();
    dimy = params.VolumeSize();
    dimz = params.VolumeSize();

    volSize = dimx*dimy*dimz;

    dim = { dimx, dimy, dimz };

    volCenter.push_back(floor((float)dimx / 2.0f) + 1.0f);
    volCenter.push_back(floor((float)dimy / 2.0f) + 1.0f);
    volCenter.push_back(floor((float)dimz / 2.0f) + 1.0f);

    rotationAndFilterIndex = params.RotationAndFilterIndex();

    partialAverage.resize(volSize);
    partialWedgeMask.resize(volSize);

    numberOfSubtomos = 0;

    realStartinIndex=params.RealStartIndex();

    if(params.SplitIntoEvenOdd())
        fsc=new FSC(params);
}

void sa_Alignment::setRotations()
{
    if(params.ChangingRotations() || rotations.size()==0)
    {
        rotations.clear();
        phiAngles.clear();
        Quaternion::generateConeRotations(rotations,phiAngles,params.ConeAngle()[rotationAndFilterIndex],params.ConeSampling()[rotationAndFilterIndex], params.InplaneAngle()[rotationAndFilterIndex], params.InplaneSampling()[rotationAndFilterIndex]);
    }

    logger.printOut("Total number of rotations for this iteration: " + to_string(rotations.size()*phiAngles.size()));
    logger.printOut("Out of this " + to_string(rotations.size()) + " cone rotations and " + to_string(phiAngles.size()) + " in-plane rotations" );
}


void sa_Alignment::preprocessMotl()
{
    if(params.CleanOutOfBoundsParticles())
    {
        motl->cleanClass(-2);

        logger.printOut("\n\nParticles that were out of bounds were removed from the motl.");
        logger.printOut("Current number of particles is " + to_string(motl->getNumberOfAllParticles()));
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(params.CleanByMeanGreyValue())
    {
        float sigma = motl->cleanByMeanGreyValue();

		logger.printOut("\n\nStandard deviation of subtomograms means is: " + to_string(sigma));
        logger.printOut("Particles with mean higher or lower than std of all subtomograms were removed from motl.");
        logger.printOut("Current number of particles is " + to_string(motl->getNumberOfAllParticles()));
    }

    MPI_Barrier(MPI_COMM_WORLD);

	string updatedMotlName=createName(params.MotlOriginalName(), realStartinIndex, ".em");
    motl->writeFullData(updatedMotlName,true);

	logger.printOut("Motive list with all changes will was saved as: " +  updatedMotlName);
}

void sa_Alignment::extractSubtomos()
{
    if(!params.ExtractSubtomos())
        return;

    ExtractSubtomos* extract = new ExtractSubtomos(params,motl->Data(),processID);
    extract->run();

    motl->setData(extract->getMotlData());

    delete extract;
}

void sa_Alignment::createStartingReference()
{
    if(!params.CreateReference())
        return;

    computePartialAverage(params.StartingIndex());
    computeFinalAverage(params.StartingIndex());

    MPI_Barrier(MPI_COMM_WORLD);

    combineEvenOdd(params.StartingIndex());
}

void sa_Alignment::setBandPassFilter()
{
    //the filter is the same and was already computed in the first iteration -> return;
    if (!params.ChangingBandPass() && bandPassMaskReduced.size() != 0)
        return;

    size_t fftz = floor(dimz / 2.0f) + 1;
    size_t fftzSize = fftz*dimx*dimy;

    float lowPassValue;

    if(params.LowPass()[rotationAndFilterIndex]<=0)
    {
        lowPassValue = computeAddaptiveFilter();
    }
    else
    {
        lowPassValue = params.LowPass()[rotationAndFilterIndex];
    }

    logger.printOut("Low-pass filter for this iteration is: " + to_string(lowPassValue));
    logger.printOut("High-pass filter for this iteration is: " + to_string(params.HighPass()[rotationAndFilterIndex]));

    vector<float> lowpass = MaskRoutines::createSphericalMask(dimx, dimy, dimz, lowPassValue, params.LowPassSigma());
    vector<float> highpass = MaskRoutines::createSphericalMask(dimx, dimy, dimz, params.HighPass()[rotationAndFilterIndex], params.HighPassSigma());
    lowpass = lowpass-highpass;

    // Inverse FFT shift
    FFTRoutines::fftshift(lowpass, dimx, dimy, dimz, -1,0);


    if (bandPassMaskReduced.size() != fftzSize)
        bandPassMaskReduced.resize(fftzSize);

    for (size_t x = 0; x < fftz; x++)
        for (size_t y = 0; y < dimy; y++)
            for (size_t z = 0; z < dimz; z++)
                bandPassMaskReduced[x + y*fftz + z*fftz*dimy] = lowpass[x + y*dimx + z*dimx*dimy];

}

void sa_Alignment::prepareReference(unsigned int index)
{
    // Read in reference

    string refName;

    if(index==params.StartingIndex() && !params.CreateReference())
        refName=createName(params.RefName(), index, ".em");
    else
        refName=createName(params.RefBaseName(), index, ".em");

    Volume unmaskedRef(refName,dim);

    // Copy the initial reference to the working folder for future check ups
    if(processID==0 && index==params.StartingIndex()  && !params.CreateReference() )
        emFile::write(createName(params.RefBaseName(), index, "_initial.em"),unmaskedRef.Data(),dimx,dimy,dimz);

    MaskRoutines::symmetrizeVolume(unmaskedRef, dimx, dimy, dimz, params.Symmetry());

    if(params.Symmetry()!=1)
        emFile::write(createName(params.RefBaseName(), index, "_symmetrized.em"),unmaskedRef.Data(),dimx,dimy,dimz);

    // Mask the reference with the alignment mask

    ref = new Particle(unmaskedRef.Data(), dim);

    if (!params.UseRosemanCC())
        ref->maskAndNormalize(mask->Data(), volSignal,false);

}


void sa_Alignment::prepareMasks()
{

    // Init wedge mask
    vector<size_t> dim = { dimx, dimy, dimz };
    wedgeMask = WedgeMask::create(params,volSize, dim);
    if (params.WedgeType() == 0)
        binarizeWedgeMask = true;
    else
        binarizeWedgeMask = false;


    // Only reference creation, no masks necessary
    if(params.NumberOfIterations()==0)
        return;

    // Read in masks
    ccMask = new Mask(params.CcMaskName(),dim);
    mask = new Mask(params.MaskName(),dim);
    
    if (params.UseRosemanCC())
    {
        mask->binarize();
    }

    volSignal = vec_sum(mask->Data());
}

void sa_Alignment::prepareMotl()
{
    string motlName;

    if(realStartinIndex==params.StartingIndex())
        motlName=createName(params.MotlName(), realStartinIndex, ".em");
    else
        motlName=createName(params.MotlBaseName(), realStartinIndex, ".em");

    motl = new Motl(motlName,numberOfProcesses,processID);

    subtomosPerMotl = motl->getNumberOfParticles();

    logger.printOut("Number of particles in motl " + motlName + " is " + to_string(motl->getOriginalNumberOfParticles()), true);


    if(params.MotlBinFactor()!=1.0f)
    {
        motl->unbin(params.MotlBinFactor());
        logger.printOut("Coordinates in motl were scaled by factor " + to_string(params.MotlBinFactor()), true);
    }
}

void sa_Alignment::clean(unsigned int iteration)
{
    remove(createName(params.MotlBaseName(), iteration, processID, ".em").c_str());
    remove(createName(params.RefBaseName(), iteration, processID, ".em").c_str());
    remove(createName(params.WedgeBaseName(), iteration, processID, ".em").c_str());
}

bool sa_Alignment::runAlignment(unsigned int iteration)
{
    if(params.Rerun()!=0)
        return true;

    try
    {
        motl->read(createName(params.MotlBaseName(), iteration, processID, ".em"));
        emFile::read(createName(params.RefBaseName(), iteration, processID, ".em"),partialAverage);
        emFile::read(createName(params.WedgeBaseName(), iteration, processID, ".em"),partialWedgeMask);

        numberOfSubtomos = motl->getActiveParticles(params.SubtomoClass());
    }
    catch(exception& e)
    {
        return true;
    }

    return false;

}

void sa_Alignment::computeFinalStats(unsigned int iteration)
{
    if(params.NumberOfIterations()==0)
        return;

    vector<float> finalConvergenceStats;

    if(globalProcessID==0)
    {
        finalConvergenceStats.resize(convergenceStats.size());
    }

    int totalNumberOfSubtomos;

    MPI_Reduce(convergenceStats.data(), finalConvergenceStats.data(), convergenceStats.size(), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&numberOfSubtomos, &totalNumberOfSubtomos, 1 , MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    for(unsigned int i=0; i<convergenceStats.size();i++)
    {
        convergenceStats[i] = 0.0f;
    }

    if (globalProcessID != 0)
        return;

    for(unsigned int i=0; i<finalConvergenceStats.size();i++)
    {
        finalConvergenceStats[i] = sqrt(finalConvergenceStats[i]/(float)totalNumberOfSubtomos);
    }

    logger.printOut("\nRMSE x shift: " + to_string(finalConvergenceStats[0]));
    logger.printOut("RMSE y shift: " + to_string(finalConvergenceStats[1]));
    logger.printOut("RMSE z shift: " + to_string(finalConvergenceStats[2]));
    logger.printOut("RMSE rotation: " + to_string(finalConvergenceStats[3]));
    logger.printOut("RMSE angular distance: " + to_string(finalConvergenceStats[4]));
    logger.printOut("RMSE in-plane rotation: " + to_string(finalConvergenceStats[5]));
}


void sa_Alignment::computeFinalAverage(unsigned int iteration)
{
    vector<float> finalAverage;
    vector<float> finalWedge;
    vector<float> finalMotl;


    MPIRoutines::reduceData(partialAverage, finalAverage, volSize,processID,numberOfProcesses,activeComm);
    MPIRoutines::reduceData(partialWedgeMask, finalWedge, volSize,processID,numberOfProcesses,activeComm);

    int totalNumberOfSubtomos=MPIRoutines::reduceNumber(numberOfSubtomos, processID, numberOfProcesses, activeComm);
    //MPI_Reduce(&numberOfSubtomos, &totalNumberOfSubtomos, 1 , MPI_INT, MPI_SUM, 0, activeComm);

    MPIRoutines::gatherData(motl->getDataAsVector(),finalMotl,processID,numberOfProcesses,activeComm);

    if (processID != 0)
        return;

    // Divide summed wedge-weight by number of input subtomograms
    vec_div(finalAverage,totalNumberOfSubtomos);
    vec_div(finalWedge,totalNumberOfSubtomos);

    vector<complex<double>> fftAverage(volSize);
    vector<float> shiftWedge(volSize);
    vector<float> shiftMask(volSize);

    FFTRoutines::forwardComplex3DTransform(dimx,dimy,dimz,finalAverage, fftAverage);
    FFTRoutines::fftshift(finalWedge,shiftWedge,dimx,dimy,dimz,-1);
    vr::vec_div(fftAverage,shiftWedge,EPS);
    vector<float> finalMask = MaskRoutines::createSphericalMask(dimx, dimy, dimz, floor(dimx/2)-3, 0.0f);
    FFTRoutines::fftshift(finalMask, shiftMask, dimx, dimy, dimz, -1);
    vr::vec_mult(fftAverage,shiftMask);
    FFTRoutines::inverseComplex3DTransform(dimx,dimy,dimz,fftAverage, finalAverage);

    //emFile::write(createName(params.MotlBaseName(), iteration, ".em"),finalMotl,20,finalMotlSize/20,1);
    Motl::write(finalMotl,createName(params.MotlBaseName(), iteration, ".em"));
    emFile::write(createName(params.RefBaseName(), iteration, ".em"),finalAverage,dimx,dimy,dimz);
    emFile::write(createName(params.WedgeBaseName(), iteration, ".em"),finalWedge,dimx,dimy,dimz);

}

void sa_Alignment::combineEvenOdd(unsigned int iteration)
{
    if(globalProcessID!=0 || !params.SplitIntoEvenOdd())
        return;

    vector<float> volEven,volOdd;

    emFile::read(createName(params.RefOriginalName()+"_odd", iteration,".em"),volOdd);
    emFile::read(createName(params.RefOriginalName()+"_even", iteration,".em"),volEven);

    vector<float> combinedAverage;
    combinedAverage=volEven+volOdd;
    vec_div(combinedAverage,2.0f);

    emFile::write(createName(params.RefOriginalName(), iteration, ".em"),combinedAverage,dimx,dimy,dimz);

    vector<MotlEntry> motlEven,motlOdd;

    Motl::read(motlOdd,createName(params.MotlOriginalName()+"_odd", iteration,".em"));
    Motl::read(motlEven,createName(params.MotlOriginalName()+"_even", iteration,".em"));

    vector<MotlEntry> combinedMotl;

    Motl::combineEvenOdd(motlEven,motlOdd,combinedMotl);
    Motl::write(combinedMotl,createName(params.MotlOriginalName(), iteration, ".em"));
	
	computeFSC(volEven,volOdd,iteration);	
}

void sa_Alignment::writeIntermediateOutputs(unsigned int iteration)
{
    motl->write(createName(params.MotlBaseName(), iteration, processID, ".em"));
    emFile::write(createName(params.RefBaseName(), iteration, processID, ".em"), partialAverage, dimx, dimy, dimz);
    emFile::write(createName(params.WedgeBaseName(), iteration, processID, ".em"), partialWedgeMask, dimx, dimy, dimz);

}

string sa_Alignment::createName(string basename, unsigned int number, string typeExtension)
{
    stringstream output;
    output << basename << "_" << number << typeExtension;

    return output.str();
}

string sa_Alignment::createName(string basename, unsigned int iteration, unsigned int number, string typeExtension)
{
    stringstream output;
    output << basename << "_" << iteration << "_" << number << typeExtension;

    return output.str();
}

float sa_Alignment::findSubpixelPeak(vector<float>& volume, vector<float>& shift)
{
    float peak = -FLT_MAX;
    size_t x = 0;
    size_t y = 0;
    size_t z = 0;

    for (size_t k = 0; k < dimz; k++)
    {
        for (size_t j = 0; j < dimy; j++)
        {
            for (size_t i = 0; i < dimx; i++)
            {
                if (volume[i + j*dimx + k*dimx*dimy] > peak)
                {
                    peak = volume[i + j*dimx + k*dimx*dimy];
                    x = i;
                    y = j;
                    z = k;
                }
            }
        }
    }

    shift[0] = x + 1.0f;
    shift[1] = y + 1.0f;
    shift[2] = z + 1.0f;

    if (x == 0 || x == (dimx - 1) || y == 0 || y == (dimy - 1) || z == 0 || z == (dimz - 1))
        return peak;

    float dx = (volume[(x + 1) + y*dimx + z*dimx*dimy] - volume[(x - 1) + y*dimx + z*dimx*dimy]) / 2.0f;
    float dy = (volume[x + (y + 1)*dimx + z*dimx*dimy] - volume[x + (y - 1)*dimx + z*dimx*dimy]) / 2.0f;
    float dz = (volume[x + y*dimx + (z + 1)*dimx*dimy] - volume[x + y*dimx + (z - 1)*dimx*dimy]) / 2.0f;

    float dxx = volume[(x + 1) + y*dimx + z*dimx*dimy] + volume[(x - 1) + y*dimx + z*dimx*dimy] - 2.0f * peak;
    float dyy = volume[x + (y + 1)*dimx + z*dimx*dimy] + volume[x + (y - 1)*dimx + z*dimx*dimy] - 2.0f * peak;
    float dzz = volume[x + y*dimx + (z + 1)*dimx*dimy] + volume[x + y*dimx + (z - 1)*dimx*dimy] - 2.0f * peak;

    float dxy = (volume[(x + 1) + (y + 1)*dimx + z*dimx*dimy] + volume[(x - 1) + (y - 1)*dimx + z*dimx*dimy] - volume[(x + 1) + (y - 1)*dimx + z*dimx*dimy] - volume[(x - 1) + (y + 1)*dimx + z*dimx*dimy]) / 4.0f;
    float dxz = (volume[(x + 1) + y*dimx + (z + 1)*dimx*dimy] + volume[(x - 1) + y*dimx + (z - 1)*dimx*dimy] - volume[(x + 1) + y*dimx + (z - 1)*dimx*dimy] - volume[(x - 1) + y*dimx + (z + 1)*dimx*dimy]) / 4.0f;
    float dyz = (volume[x + (y + 1)*dimx + (z + 1)*dimx*dimy] + volume[x + (y - 1)*dimx + (z - 1)*dimx*dimy] - volume[x + (y - 1)*dimx + (z + 1)*dimx*dimy] - volume[x + (y + 1)*dimx + (z - 1)*dimx*dimy]) / 4.0f;

    float det = dxx * (dzz * dyy - dyz * dyz) - dxy * (dzz * dxy - dyz * dxz) + dxz * (dyz * dxy - dyy * dxz);

    if (fabs(det) == 0.0f )
    {
        return peak;
    }


    float m11 =  (dzz * dyy - dyz * dyz) / det;
    float m12 = -(dzz * dxy - dyz * dxz) / det;
    float m13 =  (dyz * dxy - dyy * dxz) / det;

    float m21 = -(dzz * dxy - dxz * dyz) / det;
    float m22 =  (dzz * dxx - dxz * dxz) / det;
    float m23 = -(dyz * dxx - dxy * dxz) / det;

    float m31 =  (dyz * dxy - dxz * dyy) / det;
    float m32 = -(dyz * dxx - dxz * dxy) / det;
    float m33 =  (dyy * dxx - dxy * dxy) / det;


    float x1 = -(m11 * dx + m12 * dy + m13 * dz);
    float x2 = -(m21 * dx + m22 * dy + m23 * dz);
    float x3 = -(m31 * dx + m32 * dy + m33 * dz);

    if ((fabs(x1) > 1.0f) || (fabs(x2) > 1.0f) || (fabs(x3) > 1.f))
    {
        return peak;
    }
    else
    {
        peak = peak + dx*x1 + dy*x2 + dz*x3 + dxx*(x1 *x1)/2.0f
            + dyy*(x2 *x2)/2.0f + dzz*(x3 *x3)/2.0f + x1*x2*dxy
            + x1*x3*dxz + x2*x3*dyz;

        shift[0] = x + x1 + 1;
        shift[1] = y + x2 + 1;
        shift[2] = z + x3 + 1;

    }

    return peak;
}

