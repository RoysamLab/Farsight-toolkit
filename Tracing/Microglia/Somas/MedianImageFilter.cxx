#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMedianImageFilter.h"

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImageFile   outputImageFile" << std::endl;
    return EXIT_FAILURE;
    }

  typedef   unsigned char  PixelType;
  typedef itk::Image< PixelType,  3 > ImageType;

  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

  typedef itk::MedianImageFilter< ImageType, ImageType >  FilterType;

  FilterType::Pointer filter = FilterType::New();

  ImageType::SizeType indexRadius;
  indexRadius[0] = 1; // radius along x
  indexRadius[1] = 1; // radius along y
  indexRadius[2] = 1; // radius along z

  filter->SetRadius( indexRadius );

  filter->SetInput( reader->GetOutput() );
  ImageType::Pointer output = filter->GetOutput();
  filter->Update();

  output->SetOrigin(reader->GetOutput()->GetOrigin());
  output->SetSpacing(reader->GetOutput()->GetSpacing());

  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

