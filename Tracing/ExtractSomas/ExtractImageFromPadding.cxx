#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"


int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile  outputImageFile " << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image< unsigned char,  3 > ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  const char * inputFilename  = argv[1];
  const char * outputFilename = argv[2];
  
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  
  reader->SetFileName( inputFilename  );
  writer->SetFileName( outputFilename );
  
  typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractFilterType;
  ExtractFilterType::Pointer filter = ExtractFilterType::New();
  
  reader->Update();
  ImageType::RegionType inputRegion =
           reader->GetOutput()->GetLargestPossibleRegion();

  ImageType::SizeType size = inputRegion.GetSize();
  size[0] -= 2;
  size[1] -= 2;
  size[2] -= 2;
  
  ImageType::IndexType start = inputRegion.GetIndex();
  start[0] = 1;
  start[1] = 1;
  start[2] = 1;
  
  ImageType::RegionType desiredRegion;
  desiredRegion.SetSize(  size  );
  desiredRegion.SetIndex( start );
  filter->SetExtractionRegion( desiredRegion );
  
  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );

  try 
    { 
    writer->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 

  return EXIT_SUCCESS;
}
