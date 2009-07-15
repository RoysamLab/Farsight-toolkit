#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                   Demo: Calculating  3D Morphological, Edge and Texture features                              /////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "itkImage.h"


#include"itkBinarizeFilter.h"        //Binarization Filter
#include "itkFindObjFilter.h"        //Segmentation Filter
#include "itkObj2FeatureFilter.h"    //Protein Image Feature Calculation Filter
#include "itkObj2FeatureDNAFilter.h" //DNA Image Feature Calculation Filter
#include "itkEdgeFeaturesFilter.h"   //Edege Feature Calculation Filter

#include "itkCoOccurrenceFilter.h"     // Co-Occurrence Filter
#include "itkCO2Feature3DFilter.h"     // Co-Occurrence to Features Filter
#include "itkTextureFeaturesFilter.h"  //Texture Feature Calculation Filter

#include "itkSLF9FeatureCalFilter.h"  //SLF9 Feature Calculation Filter  (for Morphological Features which are also included in SLF11 )
#include "itkSLF11FeatureCalFilter.h "//SLF11 Feature Calculation Filter (for Edge and Texture Features)

#include "itkImageSeriesReader.h"
#include "itkNumericSeriesFileNames.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVector.h"


int main( int argc, char* argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  " << std::endl;
    std::cerr << " filenam filenameDNA first last output" << std::endl;     
    return EXIT_FAILURE;
    }
  
  typedef itk::Image<short, 3>                         ImageType;
  typedef itk::Image<double, 3>                        ImageType1;
  typedef itk::Image<unsigned char, 3>                IOImageType;
  typedef itk::Image<short, 2>                   outputImageType;
  typedef itk::ImageSeriesReader<IOImageType>           ReaderType;
  typedef ImageType::IndexType                         IndexType;
  typedef itk::NumericSeriesFileNames          NameGeneratorType;
  typedef itk::RescaleIntensityImageFilter< IOImageType, ImageType > RescaleIOImageType;
  typedef itk::SLF9FeatureCalFilter<ImageType,ImageType> SLF9FeatureCalFilterType;
  typedef itk::SLF11FeatureCalFilter<ImageType,ImageType1> SLF11FeatureCalFilterType;
  typedef itk::Vector<float, 28> VectorType;
 // typedef itk::Vector<float, 2> VectorType1;

  SLF9FeatureCalFilterType::Pointer filter   = SLF9FeatureCalFilterType::New();
  SLF11FeatureCalFilterType::Pointer filter1 = SLF11FeatureCalFilterType::New();
  RescaleIOImageType::Pointer rescaler    = RescaleIOImageType::New();
  RescaleIOImageType::Pointer rescalerDNA = RescaleIOImageType::New();

 // int first = 1;
 // int last  = 14;
  int first = atoi(argv[3]);
  int last  = atoi(argv[4]);

  std::string str1("%d.tif");
  std::string fileName = argv[1] + str1;
  std::string fileNameDNA = argv[2] + str1;

  // Read Image Series
  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
  nameGenerator->SetSeriesFormat( fileName );
  nameGenerator->SetStartIndex( first );
  nameGenerator->SetEndIndex( last );
  nameGenerator->SetIncrementIndex( 1 );

  NameGeneratorType::Pointer nameGeneratorDNA = NameGeneratorType::New();
  nameGeneratorDNA->SetSeriesFormat( fileNameDNA );
  nameGeneratorDNA->SetStartIndex( first );
  nameGeneratorDNA->SetEndIndex( last );
  nameGeneratorDNA->SetIncrementIndex( 1 );
  
  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer readerDNA = ReaderType::New();
  reader->SetFileNames( nameGenerator->GetFileNames() );
  readerDNA->SetFileNames( nameGeneratorDNA->GetFileNames() );

  // Rescale Intensity
  rescaler->SetInput(reader->GetOutput());
  rescalerDNA->SetInput(readerDNA->GetOutput());
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescalerDNA->SetOutputMinimum(0);
  rescalerDNA->SetOutputMaximum(255);
  rescaler->Update();
  rescalerDNA->Update();

  // SL9 Feature Calculation Filter, returns a Feature Vector
  std::cout<<"||||Calculating Morphological Features ...................."<<std::endl;
  filter->SetInput(rescaler->GetOutput());              // Input Protein Image
  filter->SetDNAInput(rescalerDNA->GetOutput());        // Input DNA Image
 // filter->SetThreshold(atoi(argv[1]));
  filter->SetThreshold(2);
  filter->Update();
  VectorType FeatureVector = filter->GetFeatureVector();
 


  
  std::cout<<"----------------------------3D Morphological Features--------------------------"<<std::endl;
  for(int i=0; i<28; i++) 
  {
	 std::cout<<"3D-SLF9."<<i+1<<": "<<FeatureVector[i]<<std::endl;
  }
  // std::cout<<"----------------------------3D Morphological Features----------------------------"<<std::endl;

  std::cout<<"||||Calculating Edge and Texture Features......................"<<std::endl;
  filter1->SetInput(rescaler->GetOutput());
  filter1->Update();
  VectorType FeatureVector1 = filter1->GetFeatureVector();


  std::cout<<"----------------------------3D Edge Features------------------------------"<<std::endl;
  for(int i=0; i<2; i++) 
  {
	 std::cout<<"3D-SLF11."<<i + 15<<": "<<FeatureVector1[i]<<std::endl;
  }
  //  std::cout<<"----------------------------SLF9 Feature Set----------------------------"<<std::endl;


  std::cout<<"----------------------------3D Texture Features----------------------------"<<std::endl;
  for(int i=2; i< 28; i++) 
  {
	 std::cout<<"3D-SLF11."<<i + 15 <<": "<<FeatureVector1[i]<<std::endl;
  }

  std::ofstream myfile(argv[5]);
  if (myfile.is_open())
  {
  //   myfile<<"3D Morphological Features\n";
   for(int i = 0; i< 28 ; i++)
   {
   myfile << FeatureVector[i] <<"\n";
   }
  //   myfile<<"3D Edge Features\n";
   for(int i = 0; i< 2 ; i++)
   {
   myfile << FeatureVector1[i] <<"\n";
   }
  //   myfile<<"3D Texture Features\n";
   for(int i = 2; i< 28 ; i++)
   {
   myfile << FeatureVector1[i] <<"\n";
   }
  myfile.close();
  }
  else std::cout << "Unable to open file";

  return 0;
}