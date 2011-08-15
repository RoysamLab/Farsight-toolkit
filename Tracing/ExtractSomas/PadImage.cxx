#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPasteImageFilter.h"

int main(int argc, char **argv)
{

  if( argc < 3 )
    {
    std::cerr << "usage: " << argv[0] << " inputImage outputImage" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image< unsigned char, 3 > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::PasteImageFilter< ImageType, ImageType > PasteFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  //read in an image
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::Pointer originalImage = reader->GetOutput();
  ImageType::SizeType originalSize =
    originalImage->GetLargestPossibleRegion().GetSize();
 
  //allocate a new padded image that's slightly bigger
  ImageType::Pointer paddedImage = ImageType::New();
  
  ImageType::IndexType start;
  start[0] =   0;
  start[1] =   0;
  start[2] =   0;
  
  ImageType::SizeType  paddedSize;
  paddedSize[0]  = originalSize[0] + 2;
  paddedSize[1]  = originalSize[1] + 2;
  paddedSize[2]  = originalSize[2] + 2;
  
  ImageType::RegionType region;
  region.SetSize( paddedSize );
  region.SetIndex( start );
  
  paddedImage->SetRegions( region );
  paddedImage->Allocate();
  paddedImage->FillBuffer(255);

  //paste the original image into it
  start[0] =   1;
  start[1] =   1;
  start[2] =   1;

  PasteFilterType::Pointer pasteFilter = PasteFilterType::New();
  pasteFilter->SetSourceImage( originalImage );
  pasteFilter->SetDestinationImage( paddedImage );
  pasteFilter->SetDestinationIndex( start );
  pasteFilter->SetSourceRegion( originalImage->GetBufferedRegion() );

  //write it out
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( pasteFilter->GetOutput() );
  writer->Update();

  return 0;
}
