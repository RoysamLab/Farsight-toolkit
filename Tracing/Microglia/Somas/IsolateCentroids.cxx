//image read/type classes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"

//label object classes
#include "itkShapeLabelObject.h"
#include "itkShapeLabelMapFilter.h" 
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"
#include "itkLabelMapToRGBImageFilter.h"

//for min/max limits
#include <float.h>
#include <limits.h>


#include <iostream>
#include <string>
using std::ofstream;
using std::endl;
using std::string;

int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " binarySomaImageInput minimumCentroidVolume centroidsTextOutput somaOutputImage" << std::endl;
    return EXIT_FAILURE;
    }

  typedef unsigned short      InputPixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< InputPixelType, Dimension >    ImageType;
  typedef itk::Image< itk::RGBPixel< unsigned char >, Dimension > ColorImageType;

  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ColorImageType >  ColorWriterType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  
  typedef unsigned long LabelType;
  typedef itk::ShapeLabelObject< LabelType, Dimension > LabelObjectType;
  typedef itk::LabelMap< LabelObjectType > LabelMapType;

  ReaderType::Pointer reader = ReaderType::New();

  const char * inputFilename  = argv[1];
  int minimumVolume = atoi(argv[2]);
  const char * textOutputFilename = argv[3];
  const char * imageOutputFilename = argv[4];

  std::cout << "minimum volume required to be considered a soma: " <<
    minimumVolume << std::endl;

  reader->SetFileName( inputFilename  );

  //we use a BinaryImageToShapeLabelMapFilter to convert the binary
  //input image into a collection of objects
  typedef itk::BinaryImageToShapeLabelMapFilter< ImageType, LabelMapType >
    ConverterType;
  ConverterType::Pointer converter = ConverterType::New();
  converter->SetInputForegroundValue(255);
  converter->SetInput(reader->GetOutput());
  LabelMapType::Pointer Somas = converter->GetOutput();
  converter->Update();

  //initialize output image
  ImageType::Pointer outputImage = ImageType::New();
  ImageType::IndexType start;
  start[0] =   0;
  start[1] =   0;
  start[2] =   0;
  ImageType::RegionType region;
  region.SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  region.SetIndex( start );
  outputImage->SetRegions( region );
  outputImage->Allocate();
  outputImage->FillBuffer(0);

  double minElongation = DBL_MAX;
  double maxElongation = DBL_MIN;
  double avgElongation = 0.0;
  double minFlatness = DBL_MAX;
  double maxFlatness = DBL_MIN;
  double avgFlatness = 0.0;
  double minRadius = DBL_MAX;
  double maxRadius = DBL_MIN;
  double avgRadius = 0.0;
  unsigned long minSize = ULONG_MAX;
  unsigned long maxSize = 0;
  unsigned long avgSize = 0;

  ofstream outfile(textOutputFilename);
  outfile.precision(1);
  unsigned int numSomas = Somas->GetNumberOfLabelObjects();
  std::cout << "Before pruning, I count " << numSomas << " distinct somas." << std::endl; 
  for(unsigned int label=1; label<= numSomas; ++label)
    {
    const LabelObjectType * labelObject = Somas->GetLabelObject(label);
    std::cout << "Blob #" << label << " has a volume of " <<
      labelObject->GetPhysicalSize() << std::endl;
    if(labelObject->GetPhysicalSize() < minimumVolume)
      {
      //skip small blobs: they aren't real somas
      Somas->RemoveLabel(label);
      continue;
      }

    const LabelObjectType::CentroidType centroid = labelObject->GetCentroid();
    ImageType::IndexType pixelIndex;
    reader->GetOutput()->TransformPhysicalPointToIndex( centroid, pixelIndex );
    outfile << std::fixed << pixelIndex[0] << " " << pixelIndex[1] << " " <<
      pixelIndex[2] << endl;

    double elongation = labelObject->GetElongation();
    if(elongation < minElongation)
      {
      minElongation = elongation;
      }
    if(elongation > maxElongation)
      {
      maxElongation = elongation;
      }

    double flatness = labelObject->GetFlatness();
    if(flatness < minFlatness)
      {
      minFlatness = flatness;
      }
    if(flatness > maxFlatness)
      {
      maxFlatness = flatness;
      }

    double radius = labelObject->GetEquivalentSphericalRadius();
    if(radius < minRadius)
      {
      minRadius = radius;
      }
    if(radius > maxRadius)
      {
      maxRadius = radius;
      }

    unsigned long size = labelObject->GetPhysicalSize();
    if(size < minSize)
      {
      minSize = size;
      }
    if(size > maxSize)
      {
      maxSize = size;
      }

    avgElongation += elongation;
    avgFlatness += flatness;
    avgRadius += radius;
    avgSize += size;
    }
  
  numSomas = Somas->GetNumberOfLabelObjects();
  std::cout << "I count " << numSomas << " distinct somas." << std::endl; 
  outfile.close();
  avgElongation /= numSomas;
  avgFlatness /= numSomas;
  avgRadius /= numSomas;
  avgSize /= numSomas;
 
  std::cout << "\t\t(min\tavg\tmax)" << std::endl;
  std::cout << "Elongation\t" << minElongation << "\t" << avgElongation << "\t"
            << maxElongation << std::endl;
  std::cout << "Flatness\t" << minFlatness << "\t" << avgFlatness << "\t"
            << maxFlatness << std::endl;
  std::cout << "Radius\t\t" << minRadius << "\t" << avgRadius << "\t"
            << maxRadius << std::endl;
  std::cout << "Size\t\t" << minSize << "\t" << avgSize << "\t" << maxSize
            << std::endl;

  typedef itk::LabelMapToRGBImageFilter< LabelMapType, ColorImageType >
    ColorerType;
  ColorerType::Pointer colorer = ColorerType::New();
  colorer->SetInput(Somas);
  
  typedef itk::LabelMapToBinaryImageFilter< LabelMapType, ImageType >
    BinarizerType;
  BinarizerType::Pointer binarizer = BinarizerType::New();
  binarizer->SetInput(Somas);
  
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( imageOutputFilename );
  writer->SetInput( binarizer->GetOutput() );
  writer->Update();
  
  string colorOutputFilename(argv[4]);
  size_t found = colorOutputFilename.rfind(".mhd");
  if (found != string::npos)
    {
    colorOutputFilename.replace(found, 4, "-color.tif");
    }
  else
    {
    found = colorOutputFilename.rfind(".tif");
    if (found != string::npos)
      {
      colorOutputFilename.replace(found, 4, "-color.tif");
      }
    }
  ColorWriterType::Pointer colorWriter = ColorWriterType::New();
  colorWriter->SetFileName( colorOutputFilename );
  colorWriter->SetInput(colorer->GetOutput());
  colorWriter->Update();

  return EXIT_SUCCESS; 
}

