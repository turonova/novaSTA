#include "fsc.h"
#include "fftRoutines.h"
#include "maskRoutines.h"
#include "vectorRoutines.h"
#include "emFile.h"
#include "exception.h"

using namespace std;

 FSC::FSC(ParameterSetup& iParams)
 {
     emFile::read(iParams.FscMask(),mask);
     symmetry=iParams.Symmetry();
     boxSize=iParams.VolumeSize();

     pixelSize=iParams.PixelSize();

     fftz = floor(boxSize / 2.0f) + 1;
     fftSize = boxSize*boxSize*fftz;

     numberOfShells=boxSize/2;
     centerOffset = floor(boxSize / 2.0f);

     fscCurve.resize(numberOfShells);
 }


float FSC::getFSCValue(float threshold)
{
     float x1,x2,y1,y2;
     float fscValue = 0.0f;

     for(size_t i=0; i< fscCurve.size(); i++)
     {
         if(fscCurve[i]<threshold)
         {
             x1=i;
             x2=i+1.0f;
             y1=fscCurve[i-1];
             y2=fscCurve[i];

             float slope = (y2-y1)/(x2-x1);

             fscValue = (float)boxSize*pixelSize/(((threshold-y1)/slope)+x1);

             break;
         }
     }

     return fscValue;
}

void FSC::printFSCCurve(string outputName)
{
    try
    {
        std::ofstream outfile;
        outfile.open(outputName.c_str());

        if (!outfile.good())
        {
            throw new ExceptionFileOpen(outputName);
        }

        for(unsigned int i=0; i<fscCurve.size(); i++)
            outfile << fscCurve[i] << endl;

        outfile.close();

    }
    catch (ExceptionFileOpen& e)
    {
        cout << e.what() << endl;
    }
}

void FSC::computeFSC(vector<float>& vol1, vector<float>& vol2)
{
      vector<complex<float>> fftVol1;
      vector<complex<float>> fftVol2;

      fftVol1.resize(fftSize);
      fftVol2.resize(fftSize);

      MaskRoutines::symmetrizeVolume(vol1, boxSize, boxSize, boxSize, symmetry);
      MaskRoutines::symmetrizeVolume(vol2, boxSize, boxSize, boxSize, symmetry);

      vr::vec_mult(vol1,mask);
      vr::vec_mult(vol2,mask);

      FFTRoutines::forwardReal3DTransform(boxSize, boxSize, boxSize, vol1, fftVol1);
      FFTRoutines::forwardReal3DTransform(boxSize, boxSize, boxSize, vol2, fftVol2);

      vector<complex<float>> complexConjugate;

      vr::vec_conj_mult(fftVol2,fftVol1,complexConjugate);
      vr::vec_conj_mult(fftVol1,fftVol1);
      vr::vec_conj_mult(fftVol2,fftVol2);

      vector<complex<float>> sumComplexConjugate(numberOfShells,complex<float>(0.0f,0.0f));
      vector<complex<float>> sumVol1(numberOfShells,complex<float>(0.0f,0.0f));
      vector<complex<float>> sumVol2(numberOfShells,complex<float>(0.0f,0.0f));


      for(size_t i=0; i<fftz; i++)
      {
          for(size_t j=0; j<boxSize; j++)
          {
              for(size_t k=0; k<boxSize; k++)
              {
                  float dist = sqrt(i*i+(j-centerOffset)*(j-centerOffset)+(k-centerOffset)*(k-centerOffset));

                  size_t idx=floor(dist);

                  if(idx>=numberOfShells)
                      continue;

                  size_t is = (i + (size_t)centerOffset+1) % fftz;
                  size_t js = (j + (size_t)centerOffset) % boxSize;
                  size_t ks = (k + (size_t)centerOffset) % boxSize;
                  size_t input_idx = is+js*fftz+ks*boxSize*fftz;

                  if(i==0)
                  {
                      sumComplexConjugate[idx] +=complexConjugate[input_idx];
                      sumVol1[idx] += fftVol1[input_idx];
                      sumVol2[idx] += fftVol2[input_idx];
                  }
                  else
                  {
                     sumComplexConjugate[idx] +=2.0f*complexConjugate[input_idx];
                     sumVol1[idx] += 2.0f*fftVol1[input_idx];
                     sumVol2[idx] += 2.0f*fftVol2[input_idx];
                  }

              }
          }
      }



      for(size_t i=0; i<fscCurve.size(); i++)
      {
          fscCurve[i]=real(sumComplexConjugate[i])/sqrt(real(sumVol1[i]*sumVol2[i]));
      }

      if(fscCurve[0]<0.0f)
          fscCurve[0]*=(-1.0f);
}


float FSC::getFSCValue(vector<float>& fsc, float threshold, unsigned int iBoxSize, float iPixelSize)
{
    float x1,x2,y1,y2;
    float fscValue = 0.0f;

    for(size_t i=0; i< fsc.size(); i++)
    {
        if(fsc[i]<threshold)
        {
            x1=i;
            x2=i+1.0f;
            y1=fsc[i-1];
            y2=fsc[i];

            float slope = (y2-y1)/(x2-x1);

            fscValue = (float)iBoxSize*iPixelSize/(((threshold-y1)/slope)+x1);

            break;
        }
    }

    return fscValue;
}

 vector<float> FSC::computeFSC(vector<float>& vol1, vector<float>& vol2, vector<float>& fscMask, unsigned int iBoxSize, unsigned int iSymmetry, bool iRandomized)
 {
     vector<float> fsc;

     size_t fftz = floor(iBoxSize / 2.0f) + 1;
     size_t fftSize = iBoxSize*iBoxSize*fftz;
     vector<complex<float>> fftVol1;
     vector<complex<float>> fftVol2;

     fftVol1.resize(fftSize);
     fftVol2.resize(fftSize);

     MaskRoutines::symmetrizeVolume(vol1, iBoxSize, iBoxSize, iBoxSize, iSymmetry);
     MaskRoutines::symmetrizeVolume(vol2, iBoxSize, iBoxSize, iBoxSize, iSymmetry);

     vr::vec_mult(vol1,fscMask);
     vr::vec_mult(vol2,fscMask);

     FFTRoutines::forwardReal3DTransform(iBoxSize, iBoxSize, iBoxSize, vol1, fftVol1);
     FFTRoutines::forwardReal3DTransform(iBoxSize, iBoxSize, iBoxSize, vol2, fftVol2);

     vector<complex<float>> complexConjugate;

     vr::vec_conj_mult(fftVol2,fftVol1,complexConjugate);
     vr::vec_conj_mult(fftVol1,fftVol1);
     vr::vec_conj_mult(fftVol2,fftVol2);

     unsigned int numberOfShells=iBoxSize/2;
     float centerOffset = floor(iBoxSize / 2.0f);

     vector<complex<float>> sumComplexConjugate(numberOfShells,complex<float>(0.0f,0.0f));
     vector<complex<float>> sumVol1(numberOfShells,complex<float>(0.0f,0.0f));
     vector<complex<float>> sumVol2(numberOfShells,complex<float>(0.0f,0.0f));


     for(size_t i=0; i<fftz; i++)
     {
         for(size_t j=0; j<iBoxSize; j++)
         {
             for(size_t k=0; k<iBoxSize; k++)
             {
                 float dist = sqrt(i*i+(j-centerOffset)*(j-centerOffset)+(k-centerOffset)*(k-centerOffset));

                 size_t idx=floor(dist);

                 if(idx>=numberOfShells)
                     continue;

                 size_t is = (i + (size_t)centerOffset+1) % fftz;
                 size_t js = (j + (size_t)centerOffset) % iBoxSize;
                 size_t ks = (k + (size_t)centerOffset) % iBoxSize;
                 size_t input_idx = is+js*fftz+ks*iBoxSize*fftz;

                 if(i==0)
                 {
                     sumComplexConjugate[idx] +=complexConjugate[input_idx];
                     sumVol1[idx] += fftVol1[input_idx];
                     sumVol2[idx] += fftVol2[input_idx];
                 }
                 else
                 {
                    sumComplexConjugate[idx] +=2.0f*complexConjugate[input_idx];
                    sumVol1[idx] += 2.0f*fftVol1[input_idx];
                    sumVol2[idx] += 2.0f*fftVol2[input_idx];
                 }

             }
         }
     }

     fsc.resize(numberOfShells);


     for(size_t i=0; i<fsc.size(); i++)
     {
         fsc[i]=real(sumComplexConjugate[i])/sqrt(real(sumVol1[i]*sumVol2[i]));
     }

     if(fsc[0]<0.0f)
         fsc[0]*=(-1.0f);

     return fsc;
 }

/*vector<float> FSC::computeFSC(vector<float>& vol1, vector<float>& vol2, vector<float>& fscMask, unsigned int iBoxSize, unsigned int iSymmetry, bool iRandomized)
{
    vector<float> fsc;

    vector<complex<float>> fftVol1;
    vector<complex<float>> fftVol2;

    fftVol1.resize(vol1.size());
    fftVol2.resize(vol2.size());

    MaskRoutines::symmetrizeVolume(vol1, iBoxSize, iBoxSize, iBoxSize, iSymmetry);
    MaskRoutines::symmetrizeVolume(vol2, iBoxSize, iBoxSize, iBoxSize, iSymmetry);

    vr::vec_mult(vol1,fscMask);
    vr::vec_mult(vol2,fscMask);

    FFTRoutines::forwardComplex3DTransform(iBoxSize, iBoxSize, iBoxSize, vol1, fftVol1);
    FFTRoutines::forwardComplex3DTransform(iBoxSize, iBoxSize, iBoxSize, vol2, fftVol2);

    vector<complex<float>> complexConjugate;

    vr::vec_conj_mult(fftVol2,fftVol1,complexConjugate);
    vr::vec_conj_mult(fftVol1,fftVol1);
    vr::vec_conj_mult(fftVol2,fftVol2);

    unsigned int numberOfShells=iBoxSize/2;

    float centerOffset = floor(iBoxSize / 2.0f);

    vector<complex<float>> sumComplexConjugate(numberOfShells,complex<float>(0.0f,0.0f));
    vector<complex<float>> sumVol1(numberOfShells,complex<float>(0.0f,0.0f));
    vector<complex<float>> sumVol2(numberOfShells,complex<float>(0.0f,0.0f));


    for(size_t i=0; i<iBoxSize; i++)
    {
        for(size_t j=0; j<iBoxSize; j++)
        {
            for(size_t k=0; k<iBoxSize; k++)
            {
                float dist = sqrt((i-centerOffset)*(i-centerOffset)+(j-centerOffset)*(j-centerOffset)+(k-centerOffset)*(k-centerOffset));

                size_t idx=floor(dist);

                if(idx>=numberOfShells)
                    continue;

                size_t is = (i + (size_t)centerOffset) % iBoxSize;
                size_t js = (j + (size_t)centerOffset) % iBoxSize;
                size_t ks = (k + (size_t)centerOffset) % iBoxSize;

                sumComplexConjugate[idx] += complexConjugate[is+js*iBoxSize+ks*iBoxSize*iBoxSize];
                sumVol1[idx] += fftVol1[is+js*iBoxSize+ks*iBoxSize*iBoxSize];
                sumVol2[idx] += fftVol2[is+js*iBoxSize+ks*iBoxSize*iBoxSize];

            }
        }
    }

    fsc.resize(numberOfShells);


    for(size_t i=0; i<fsc.size(); i++)
    {
        fsc[i]=real(sumComplexConjugate[i])/sqrt(real(sumVol1[i]*sumVol2[i]));
    }

    return fsc;
}*/
