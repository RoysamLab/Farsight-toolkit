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

#include <iostream>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

#define DATATYPEIN unsigned char
#define DATATYPEOUT unsigned char

int main(int argc, char *argv[])
{
  bool Inside=255;
  //check that the user specified the right number of command line arguments
  if(argc < 4)
    {
    cerr << argv[0] << " <input file1> <input file2> <output file> <bool>"
         << endl;
    return EXIT_FAILURE;
    }

  DATATYPEIN *volin1,*volin2;
  unsigned int sizeX, sizeY, sizeZ;

  string inputFileName1 = argv[1];
  string inputFileName2 = argv[2];

  const char *outputFileName;
  outputFileName = argv[3];
  Inside = atoi(argv[4]); 
  printf("Inside is %d ",Inside);

  FILE *outfile;
  unsigned int i,j,k;
  DATATYPEOUT *volout;
  long idx;
  int sizeExpand = 0;


  //make sure we can write to the output file
   if((outfile=fopen(outputFileName, "wb")) == NULL)
    {
    cerr << "Output file open error!" << endl;
    return EXIT_FAILURE;
    }

    //use ITK to read all non-raw images
    typedef unsigned char     PixelType;
    const   unsigned int      Dimension = 3;
    typedef itk::Image< PixelType, Dimension >    ImageType;
    typedef itk::ImageFileReader< ImageType >  ReaderType;
    ReaderType::Pointer reader1 = ReaderType::New();
    ReaderType::Pointer reader2 = ReaderType::New();
    ImageType::Pointer inputImage1 = reader1->GetOutput();
	ImageType::Pointer inputImage2 = reader2->GetOutput();
    reader1->SetFileName( inputFileName1.c_str()  );
    reader2->SetFileName( inputFileName2.c_str()  );
    try 
      { 
      reader1->Update(); 
      reader2->Update(); 
      } 
    catch( itk::ExceptionObject & err ) 
      { 
      cerr << "ExceptionObject caught!" << endl; 
      cerr << err << endl; 
      return EXIT_FAILURE;
      } 
    ImageType::RegionType region1 = inputImage1->GetBufferedRegion();
    ImageType::SizeType  size1 = region1.GetSize();
	ImageType::RegionType region2 = inputImage2->GetBufferedRegion();
    ImageType::SizeType  size2 = region2.GetSize();
    cout << "input image size: " << size1 << endl;
    sizeX = size1[0];
    sizeY = size1[1];
    sizeZ = size1[2];

	if (sizeX != size2[0] || sizeY != size2[1] || sizeZ != size2[2])
	  {
    cerr << "Output file open error!" << endl;
    return EXIT_FAILURE;
    }
	else
	{
    //manually copy the values from inputImage to volin, which we will use for
    //the rest of the MDL process
    itk::ImageRegionIterator< ImageType >
      itr1( inputImage1, inputImage1->GetBufferedRegion() );
    itr1.GoToBegin();
    long int idx = 0;
    volin1 = (DATATYPEIN*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEIN));
    while( ! itr1.IsAtEnd() )
      {
      volin1[idx] = itr1.Get();
      ++itr1;
      ++idx;
      }

	itk::ImageRegionIterator< ImageType >
    itr2( inputImage2, inputImage2->GetBufferedRegion() );
    itr2.GoToBegin();
    idx = 0;
    volin2 = (DATATYPEIN*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEIN));
    while( ! itr2.IsAtEnd() )
      {
      volin2[idx] = itr2.Get();
      ++itr2;
      ++idx;
      }
	 }

  //allocate memory for the output image
  volout = (DATATYPEOUT*)malloc(sizeX*sizeY*(sizeZ+sizeExpand*2)*sizeof(DATATYPEOUT));
  cout << "Volume Processing..." << endl;

  for (k=0; k<sizeZ; k++)
    {
    // threshold first
    for (j=0; j<sizeY; j++)
      {
      for (i=0; i<sizeX; i++)
        {
        idx = k *sizeX*sizeY + j *sizeX + i;
		if(Inside ==1 )
        volout[idx]= ( volin2[idx]*volin1[idx]/255);
		else 
        volout[idx]= ( (255-volin2[idx])*volin1[idx]/255);
        }
      }
    }   

  //write the output image and free memory
  fwrite(volout, sizeX*sizeY*sizeZ, sizeof(DATATYPEOUT), outfile);

  FILE *mhdfile;
  if((mhdfile=fopen("ObjectPreProcessing.mhd","w"))==NULL)
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
  fprintf (mhdfile,"ElementDataFile = ObjectPreProcessing.raw\n");
  fclose(mhdfile);


  


  free(volout); 
  volout=NULL;

  cout << "Done" << endl;
  return 0;
}


