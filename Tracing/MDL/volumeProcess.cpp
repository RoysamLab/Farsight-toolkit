/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

/*  Volume dataset processing
 *  accept a sequence of volumes
 *  Windows version, taken in from Linux version
 *   Author: Xiaosong Yuan, RPI
 *  Modified on Sep. 29, 2005  

 *  Input parameters
 *          1. sizeExpand   
 *          2. preproess          */

#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif

#include <iostream>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

#define DATATYPEIN unsigned char
#define DATATYPEOUT unsigned char

const unsigned char  m_NumberOfHistogramBins = 128;

DATATYPEIN *volin;
int sizeX, sizeY, sizeZ;

double  OtsuThreshold (int sizeX,int sizeY,int sizeZ);

int main(int argc, char *argv[])
{

  //check that the user specified the right number of command line arguments
  if(argc < 3)
    {
    cerr << argv[0] << " <input file> <output file>" << endl;
    cerr << argv[0] << " <raw input file> <sizeX> <sizeY> <sizeZ> <output file>"
         << endl;
    return EXIT_FAILURE;
    }

  //check if the input image is .raw
  bool rawInput = false;
  string inputFileName = argv[1];
  const char *outputFileName;
  if(inputFileName.rfind(".raw") != string::npos)
    {
    //if so, the user is also required to pass in image dimensions
    if(argc < 6)
      {
      cerr << "Usage: <raw input file> <sizeX> <sizeY> <sizeZ> <output file>" << endl;
      return EXIT_FAILURE;
      }
    rawInput = true;
    sizeX = atoi(argv[2]);
    sizeY = atoi(argv[3]);
    sizeZ = atoi(argv[4]);
    outputFileName = argv[5];
    }
  else
    {
    outputFileName = argv[2];
    }

  FILE *infile;
  FILE *outfile;
  int i,j,k;
  int ii, jj, kk;
  DATATYPEOUT *volout;
  long idx;
  double threshold;
  int sizeExpand = 0;
  DATATYPEOUT blockMax;
  int timesDilate;
  int border;

  //make sure we can write to the output file
  if((outfile=fopen(outputFileName, "wb")) == NULL)
    {
    cerr << "Output file open error!" << endl;
    return EXIT_FAILURE;
    }

  //initialize the input image
  if(rawInput)
    {
    if((infile=fopen(inputFileName.c_str(), "rb"))==NULL)
      {
      cout << "Input file open error!" << endl;
      return EXIT_FAILURE;
      }

    volin = (DATATYPEIN*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEIN));

    if (fread(volin, sizeof(DATATYPEIN), sizeX*sizeY*sizeZ, infile) < (unsigned int)(sizeX*sizeY*sizeZ))
      {
      cerr << "File size is not the same as volume size" << endl;
      return EXIT_FAILURE;
      }
    }
  else
    {
    //use ITK to read all non-raw images
    typedef unsigned char     PixelType;
    const   unsigned int      Dimension = 3;
    typedef itk::Image< PixelType, Dimension >    ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    ImageType::Pointer inputImage = reader->GetOutput();
    reader->SetFileName( inputFileName.c_str()  );
    try 
      { 
      reader->Update(); 
      } 
    catch( itk::ExceptionObject & err ) 
      { 
      cerr << "ExceptionObject caught!" << endl; 
      cerr << err << endl; 
      return EXIT_FAILURE;
      } 
    ImageType::RegionType region = inputImage->GetBufferedRegion();
    ImageType::SizeType  size = region.GetSize();
    cout << "input image size: " << size << endl;
    sizeX = size[0];
    sizeY = size[1];
    sizeZ = size[2];
    
    //manually copy the values from inputImage to volin, which we will use for
    //the rest of the MDL process
    itk::ImageRegionIterator< ImageType >
      itr( inputImage, inputImage->GetBufferedRegion() );
    itr.GoToBegin();
    long int idx = 0;
    volin = (DATATYPEIN*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEIN));
    while( ! itr.IsAtEnd() )
      {
      volin[idx] = itr.Get();
      ++itr;
      ++idx;
      }
    }

  //allocate memory for the output image
  volout = (DATATYPEOUT*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEOUT));

  cout << "Volume Processing..." << endl;

  // Pre-Processing 
  // by xiao liang, using 3 sigma theory to estimate threshold;
  double meanValue =0.0, VarianceValue =  0.0;

  // threshold first
  for (k=0; k<sizeZ; k++)
    {
    for (j=0; j<sizeY; j++) 
      {
      for (i=0; i<sizeX; i++)
        {
        idx = k *sizeX*sizeY + j *sizeX + i;
        meanValue += volin[idx];
        }
      }
    }

  meanValue = meanValue/(double)(sizeX*sizeY*sizeZ);

  // threshold first
  for (k=0; k<sizeZ; k++)
    {
    for (j=0; j<sizeY; j++) 
      {
      for (i=0; i<sizeX; i++)
        {
        idx = k *sizeX*sizeY + j *sizeX + i;
        VarianceValue += (volin[idx]-meanValue)*(volin[idx]-meanValue);
        }
      }
    }

  VarianceValue =  VarianceValue/(double)(sizeX*sizeY*sizeZ);
  VarianceValue = sqrt(VarianceValue);

  double m_threshold=OtsuThreshold (sizeX,sizeY,sizeZ);
  if (m_threshold > 7 || m_threshold < 0)
    {
    threshold =(meanValue-VarianceValue/30); 
    }
  else
    {
    threshold = m_threshold;
    }

  cout << "OTSU optimal threshold " << threshold << endl;

     for (k=0; k<(sizeZ+sizeExpand*2); k++)
      for (j=0; j<sizeY; j++)
        for (i=0; i<sizeX; i++) {

          volout[k *sizeX*sizeY + j *sizeX + i] = 0; 
        }  //initial to zeros

  for (k=0; k<sizeZ; k++)
    {
    // threshold first
    for (j=0; j<sizeY; j++)
      {
      for (i=0; i<sizeX; i++)
        {
        idx = k *sizeX*sizeY + j *sizeX + i;
        if (volin[idx] < threshold) 
          {
          volin[idx] = 0;
          }
        }
      }
    }   

  // Method 2: Dilation of the object
  timesDilate = 1;
  border = 3;
  while (timesDilate >0 )
    {
    for (k=border; k<sizeZ-border; k++)
      {
      for (j=border; j<sizeY-border; j++)
        {
        for (i=border; i<sizeX-border; i++)
          {
          blockMax = volin[k *sizeX*sizeY + j *sizeX + i];
          for (kk=-1; kk<=1; kk++)
            {
            for (jj=-1; jj<=1; jj++)
              {
              for (ii=-1; ii<=1; ii++)
                {
                if(volin[(k+kk)*sizeX*sizeY + (j+jj)*sizeX + (i+ii)] > blockMax) 
                  {
                  blockMax = volin[(k+kk)*sizeX*sizeY + (j+jj)*sizeX + (i+ii)];
                  }
                }
              }
            }
          // Keep the peak of the original intensity
          if (blockMax == volin[k *sizeX*sizeY + j *sizeX + i] && blockMax != 0)
            {
            blockMax = blockMax + 1;
            //if (blockMax > 255)   blockMax = 255;
            }
          volout[k *sizeX*sizeY + j *sizeX + i] = blockMax;
          }
        }
      }

    // copy volout back to volin for the next dilation
    for (k=0; k<sizeZ; k++) 
      {
      for (j=0; j<sizeY; j++) 
        {
        for (i=0; i<sizeX; i++)
          {
          volin[k *sizeX*sizeY + j *sizeX + i] = volout[k *sizeX*sizeY + j *sizeX + i];
          }
        }
      }
    timesDilate--;
    }

  //write the output image and free memory
  fwrite(volout, sizeX*sizeY*sizeZ, sizeof(DATATYPEOUT), outfile);

  fclose(infile);
  fclose(outfile);

  free(volin);  // by xiao
  free(volout); // by xiao
  volin=NULL;
  volout=NULL;

  cout << "Done" << endl;
  return 0;
}


// This function is designed to compute the opyimal threshold using OTSU method;
// this algoritm is implemented by xiao liang based on ITK's OTSU algorithm 
double  OtsuThreshold (int sizeX,int sizeY,int sizeZ)

{

  int i,j,k;
  double m_Threshold =0;
  double totalPixels = (double)(sizeX*sizeY*sizeZ);
 

  if ( totalPixels == 0 ) { return m_Threshold; }

  unsigned char MinValue = volin[0], MaxValue =  volin[0];
  double meanValue=0.0, varianceValue=0.0;
  int idx; // just for the location0

  for (k=0; k < sizeZ; k++)
  {  // 
   for (j=0; j < sizeY; j++) 
      {
        for (i=0; i<sizeX; i++)
        {
          idx = k *sizeX*sizeY + j *sizeX + i;
          meanValue += volin[idx];
          if (volin[idx]> MaxValue) MaxValue = volin[idx];
          if (volin[idx]< MinValue) MinValue = volin[idx];

        }
      }
    }

      cout << "Max = " << (int)MaxValue << ", Min = " << (int)MinValue << endl;
       meanValue = meanValue/totalPixels;
 
     for (k=0; k<sizeZ; k++)
    {  //
      for (j=0; j<sizeY; j++) 
      {
        for (i=0; i<sizeX; i++)
        {
          idx = k *sizeX*sizeY + j *sizeX + i;
          varianceValue += (volin[idx]-meanValue)*(volin[idx]-meanValue);
        }
      }
    } 


    if ( MinValue >= MaxValue)
    {
     m_Threshold=MinValue;
       return m_Threshold;
    
    }
   
  m_Threshold = (meanValue-varianceValue/30); 
    // this step is only initialized a good experimental value for m_Threshold, because the 3D image
    // is sparse, there are lots of zero values; 

  // create a histogram
  double relativeFrequency[m_NumberOfHistogramBins];

  for ( j = 0; j < m_NumberOfHistogramBins; j++ )
    {
       relativeFrequency[j] = 0.0;
    }

  double binMultiplier = (double) m_NumberOfHistogramBins /(double) ( MaxValue - MinValue);
  cout << "binMultiplier = " << binMultiplier << endl;
 
  
  double Voxelvalue; // temp variable
  unsigned int binNumber;

  for (k=0; k<sizeZ; k++) 
  {  // 
  for (j=0; j<sizeY; j++) 
   {
  for (i=0; i<sizeX; i++)
   {
    idx = k *sizeX*sizeY + j *sizeX + i;
        Voxelvalue = volin[idx];

        if ( Voxelvalue == MinValue ) 
         {
         binNumber = 0;
         } // end if 

        else
          {
             binNumber = (unsigned int)(((Voxelvalue - MinValue) * binMultiplier ) - 1);

             if ( binNumber == m_NumberOfHistogramBins ) // in case of rounding errors
              {
                binNumber -= 1;
              }// end if 
           }// end else

       //  printf("binNumber???? = %f  MaxValue=%f \n", binNumber+1,MaxValue);
         relativeFrequency[binNumber] += 1.0;
   
  }//
  }
   }
// test 

//  for (i=0;i<m_NumberOfHistogramBins;i++)
//    printf ( "%d bin = %f,  ",i, relativeFrequency[i]);
 
  // normalize the frequencies
  double totalMean = 0.0;
  for ( j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    relativeFrequency[j] /= totalPixels;
    totalMean += (j+1) * relativeFrequency[j];
    }


  // compute Otsu's threshold by maximizing the between-class
  // variance
  double freqLeft = relativeFrequency[0];
  double meanLeft = 1.0;
  double meanRight = ( totalMean - freqLeft ) / ( 1.0 - freqLeft );

  double maxVarBetween = freqLeft * ( 1.0 - freqLeft ) * sqrt( meanLeft - meanRight );
  int maxBinNumber = 0;

  double freqLeftOld = freqLeft;
  double meanLeftOld = meanLeft;

  for ( j = 1; j < m_NumberOfHistogramBins; j++ )
    {
    freqLeft += relativeFrequency[j];
    meanLeft = ( meanLeftOld * freqLeftOld + 
                 (j+1) * relativeFrequency[j] ) / freqLeft;
    if (freqLeft == 1.0)
      {
      meanRight = 0.0;
      }
    else
      {
      meanRight = ( totalMean - meanLeft * freqLeft ) / ( 1.0 - freqLeft );
      }
    double varBetween = freqLeft * ( 1.0 - freqLeft ) * sqrt( meanLeft - meanRight );
   
    if ( varBetween > maxVarBetween )
      {
      maxVarBetween = varBetween;
      maxBinNumber = j;
      }

    // cache old values
    freqLeftOld = freqLeft;
    meanLeftOld = meanLeft; 

    } 

    m_Threshold = double( MinValue + ( maxBinNumber + 1 ) / binMultiplier );
  cout << "m_threshold = " << m_Threshold << endl;
  return m_Threshold; 
}
