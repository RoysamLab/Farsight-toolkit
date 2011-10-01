#include "itkImage.h"
#include <itkBinaryBallStructuringElement.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile outputImageFile minRadius" << std::endl;
    return EXIT_FAILURE;
    }

  typedef   unsigned char  InputPixelType;
  typedef   unsigned char  OutputPixelType;

  typedef itk::Image< InputPixelType,  3 >   InputImageType;
  typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

  typedef itk::ImageFileReader< InputImageType  >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  typedef itk::BinaryBallStructuringElement<float, 3> StructuringElementType;
  typedef itk::BinaryMorphologicalOpeningImageFilter
    < InputImageType, OutputImageType, StructuringElementType > OpenFilterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  InputImageType::Pointer inputImage = reader->GetOutput();
  reader->Update();

  //radius (in microns) of the smallest structure that we wish to keep.
  //anything smaller will be eroded away by this filter.
  float minRadius = atof( argv[3] );

  InputImageType::SpacingType spacing = inputImage->GetSpacing();
  
  StructuringElementType structuringElement;
  
  //radius of the structuring element should be minRadius microns in each
  //dimension.  This code assumes that the image spacing is set in microns.
  StructuringElementType::RadiusType radius =
    { ( minRadius / spacing[0] ), ( minRadius / spacing[1] ),
      ( minRadius / spacing[2] ) };
  structuringElement.SetRadius( radius );
  structuringElement.CreateStructuringElement();
 
  OpenFilterType::Pointer openFilter = OpenFilterType::New();
  openFilter->SetKernel( structuringElement );
  openFilter->SetInput( inputImage );

  openFilter->SetInput( reader->GetOutput() );
  
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( openFilter->GetOutput() );
  writer->Update();
 
  return EXIT_SUCCESS;
}
