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

#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif

# define Curvature_Anistropic_Diffusion  1
#include <iostream>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkOtsuThresholdImageFilter.h" 
#include "itkRescaleIntensityImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "GCBinarization/cell_binarization.h"
#include "distTransform.h"
using std::cerr;
using std::cout;
using std::endl;
using std::string;

#define DATATYPEIN unsigned char
#define DATATYPEOUT unsigned char

const unsigned char  m_NumberOfHistogramBins = 128;

DATATYPEIN *volin;
int sizeX, sizeY, sizeZ;
float *tempimage;
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
  //double space[3]={1,1,1};
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

  FILE *infile = 0;
  FILE *outfile;
  int i,j,k;
  int ii, jj, kk;
  DATATYPEOUT *volout;
  long idx;

 
  int sizeExpand = 0;
  DATATYPEOUT blockMax;
  int timesDilate;
  int border;
  unsigned short *GraphcutResults;

  cout << "Volume Processing..." << endl;
  //make sure we can write to the output file
  if((outfile=fopen(outputFileName, "wb")) == NULL)
    {
    cerr << "Output file open error!" << endl;
    return EXIT_FAILURE;
    }

 
    typedef float     PixelType;
    const   unsigned int      Dimension = 3;
    typedef itk::Image< PixelType, Dimension >    ImageType;
    ImageType::Pointer inputImage;

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
    
	inputImage = ImageType::New();

	ImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
    origin[2] = 0.;    
	inputImage->SetOrigin( origin );

	ImageType::IndexType start;
    start[0] =   0;  
    start[1] =   0;  
	start[2] =   0;  

	ImageType::SizeType  size;
    size[0]  = sizeX;  
    size[1]  = sizeY;  
	size[2]  = sizeZ;  

    ImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    inputImage->SetRegions( region );
    inputImage->Allocate();
    inputImage->FillBuffer(0);
	inputImage->Update();
	
	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	IteratorType iterator1(inputImage,inputImage->GetRequestedRegion());
	int lng = sizeX*sizeY*sizeZ;
	
	for(int i=0; i<lng; i++)
	{
	 iterator1.Set((float)volin[i]);
	 ++iterator1;	
	}
  } // end if(rawInput)

   else
    {
  
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    inputImage = reader->GetOutput();
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
	} // end of itk image read 

 
	 // ---------------Linear Mapping --------------------//
    typedef itk::RescaleIntensityImageFilter<
               ImageType, ImageType>  RescaleFilterType;
    RescaleFilterType::Pointer    rescaleFilter    = RescaleFilterType::New();
    rescaleFilter->SetInput(inputImage);
	rescaleFilter->SetOutputMinimum(  0 );
    rescaleFilter->SetOutputMaximum( 255 );
	
	try 
      { 
       rescaleFilter->Update();
      } 
    catch( itk::ExceptionObject & err ) 
      { 
      cerr << "ExceptionObject caught!" << endl; 
      cerr << err << endl; 
      return EXIT_FAILURE;
      }

	inputImage = rescaleFilter->GetOutput();
	cout << "The Linear Mapping is done\n";

    # if Curvature_Anistropic_Diffusion 
	{
     cout << "The Curvature Diffusion is doing\n"; 
	 typedef itk::CurvatureAnisotropicDiffusionImageFilter<
               ImageType, ImageType >  MCD_FilterType;

     MCD_FilterType::Pointer MCDFilter = MCD_FilterType::New();
    
     //Initialnization,  using the paper's optimal parameters
     const unsigned int numberOfIterations = 5;
     const double       timeStep = 0.0425;
     const double       conductance = 3;
     MCDFilter->SetNumberOfIterations(numberOfIterations);
     MCDFilter->SetTimeStep( timeStep );
     MCDFilter->SetConductanceParameter( conductance );
     MCDFilter->SetInput(inputImage);
	 try 
      { 
       MCDFilter->Update();
      } 
     catch( itk::ExceptionObject & err ) 
      { 
      cerr << "ExceptionObject caught!" << endl; 
      cerr << err << endl; 
      return EXIT_FAILURE;
      } 
     inputImage=MCDFilter->GetOutput(); 
     cout << "The Curvature Diffusion is done\n";
	}
	#endif

    ImageType::RegionType region = inputImage->GetBufferedRegion();
    ImageType::SizeType  size = region.GetSize();
    cout << "input image size: " << size << endl;
    sizeX = size[0];
    sizeY = size[1];
    sizeZ = size[2];
    itk::ImageRegionIterator< ImageType >
      itr( inputImage, inputImage->GetBufferedRegion() );
    itr.GoToBegin();
    idx = 0;
    volin = (DATATYPEIN*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEIN));
    while( ! itr.IsAtEnd() )
      {
      volin[idx] = itr.Get();
      ++itr;
      ++idx;
      }
  //allocate memory for the output image
  volout = (DATATYPEOUT*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEOUT));
  
  // one pre-processing  scheme 
  
  GraphcutResults = (unsigned short*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(unsigned short));
  Neuron_Binarization_3D(volin,GraphcutResults,sizeX,sizeY,sizeZ,0,1);
  
  for (k=0; k<(sizeZ+sizeExpand*2); k++)
      for (j=0; j<sizeY; j++)
        for (i=0; i<sizeX; i++) {
          volout[k *sizeX*sizeY + j *sizeX + i] = 0; 
        }  //initial to zeros
   
  
  std::cout << "Do you think we need the distance transform to make the centerline of image become bright with higher intensity?";
  
  char tmpAnswer = 'n';
  //tmpAnswer = getchar();
  if (tmpAnswer =='Y' || tmpAnswer =='y')
  {
  unsigned char *GraphcutResultsTmp;
  GraphcutResultsTmp = (unsigned char*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(unsigned char));
  for (k=0; k<(sizeZ+sizeExpand*2); k++)
      for (j=0; j<sizeY; j++)
        for (i=0; i<sizeX; i++) 
		{
          GraphcutResultsTmp[k *sizeX*sizeY + j *sizeX + i] = (unsigned char) GraphcutResults[k *sizeX*sizeY + j *sizeX + i]; 
        }  //initial to zeros

  distTransform(GraphcutResultsTmp,sizeX,sizeY,sizeZ);
  
  for (k=0; k<sizeZ; k++)
    {
    // threshold first
     for (j=0; j<sizeY; j++)
       {
       for (i=0; i<sizeX; i++)
        {
        idx = k *sizeX*sizeY + j *sizeX + i;
        volin[idx] = GraphcutResultsTmp [idx];
        } // end for 
      } // end for 
    } // end for 
  
  free(GraphcutResultsTmp);
  GraphcutResultsTmp=NULL;
  } // end if 

  else {
 
   for (k=0; k<sizeZ; k++)
    {
    // threshold first
     for (j=0; j<sizeY; j++)
       {
       for (i=0; i<sizeX; i++)
        {
        idx = k *sizeX*sizeY + j *sizeX + i;
		
        if (GraphcutResults [idx] == 0) 
          {
          volin[idx] = 0;
          }
        }
      }
    }   

 } // end else

 free(GraphcutResults);
 GraphcutResults=NULL;

  // the secpnd Pre-Processing method, corresponding to the old version MDL 
  /*
  double threshold; 
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

  //threshold =7;
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
  */
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

  //-----  Image write 
   /*
    typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	IteratorType iterator1(inputImage,inputImage->GetRequestedRegion());
	int lng = sizeX*sizeY*sizeZ;
	
	for(int i=0; i<lng; i++)
	{
	 iterator1.Set((float)volin[i]);
	 ++iterator1;	
	}
   */
  //write the output image and free memory
  fwrite(volout, sizeX*sizeY*sizeZ, sizeof(DATATYPEOUT), outfile);
  FILE *mhdfile;
  
  if((mhdfile=fopen("volume_Processed.mhd","w"))==NULL)
    {
    cerr << "output file open error!" << endl;
    return -1;
    }
  fprintf (mhdfile,"ObjectType = Image\n");
  fprintf (mhdfile,"NDims = 3\n");
  fprintf (mhdfile,"BinaryData = True\n");
  fprintf (mhdfile,"BinaryDataByteOrderMSB = False\n");
  fprintf (mhdfile,"CompressedData = False\n");
  fprintf (mhdfile,"TransformMatrix = 1 0 0 0 1 0 0 0 1\n");
  fprintf (mhdfile,"Offset = 0 0 0\n");
  fprintf (mhdfile,"CenterOfRotation = 0 0 0\n");
  fprintf (mhdfile,"AnatomicalOrientation = RAI\n");
  fprintf (mhdfile,"ElementSpacing = 1 1 1\n");
  fprintf (mhdfile,"DimSize = %d %d %d\n",sizeX,sizeY,sizeZ);
  fprintf (mhdfile,"ElementType = MET_UCHAR\n");
  fprintf (mhdfile,"ElementDataFile = volume_Processed.raw\n");
  fclose(mhdfile);

  if (rawInput)
	  fclose(infile);
  fclose(outfile);

  free(volin);  // by xiao
  free(volout); // by xiao
  volin=NULL;
  volout=NULL;

  //cout << "Done" << endl;
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
