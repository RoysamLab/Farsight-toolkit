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

  //read input image & ensure that its spacing is set 
  ImageType::Pointer input = reader->GetOutput();
  try
    {
    input->Update();
    }
  catch (itk::ExceptionObject &err)
    {
    std::cerr << "reader Exception: " << err << std::endl;
    }

  ImageType::SpacingType spacing;
  spacing[0] = 1;
  spacing[1] = 1;
  spacing[2] = 1;

  //check if image spacing was provided on command-line
  if( argc > 3)
    {
    std::cout << "Reading image spacing from command-line" << std::endl;
    spacing[0] = atof(argv[3]);
    spacing[1] = spacing[0];
    spacing[2] = atof(argv[4]);
    input->SetSpacing(spacing);
    std::cout << "Input image's spacing now set to " << input->GetSpacing() << std::endl;
    }
  //else, if input image exhibits default spacing, ask user for correct values
  else if(input->GetSpacing() == spacing)
    {
    std::cout << "WARNING: spacing not set on input image.  Spacing should be set in microns per pixel." << std::endl;
    std::cout << "Please input spacing for X dimension now:" << std::endl;
    std::cin >> spacing[0];
    std::cout << "Y-spacing is assumed to be the same as X-spacing" << std::endl;
    spacing[1] = spacing[0];
    std::cout << "Please input spacing for Z dimension now:" << std::endl;
    std::cin >> spacing[2];
    input->SetSpacing(spacing);
    std::cout << "Input image's spacing now set to " << input->GetSpacing() << std::endl;
    }

  typedef itk::MedianImageFilter< ImageType, ImageType >  FilterType;

  FilterType::Pointer filter = FilterType::New();

  ImageType::SizeType indexRadius;
  indexRadius[0] = 1; // radius along x
  indexRadius[1] = 1; // radius along y
  indexRadius[2] = 1; // radius along z

  filter->SetRadius( indexRadius );

  filter->SetInput( input );
  ImageType::Pointer output = filter->GetOutput();
  filter->Update();

  output->SetOrigin( input->GetOrigin() );
  output->SetSpacing( input->GetSpacing() );

  writer->SetInput( output );
  writer->Update();

  return EXIT_SUCCESS;
}

