
#if defined(_MSC_VER)
#pragma warning (disable : 4786)
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif
#
#include <iostream>
#include <math.h>
#include "Daniela.h"
using namespace std;
Daniela::Daniela()
{
}
Daniela::Daniela(std::string imageFileName, int numberOfInitialClasses)
{
	M =0;
	N =0;
	
	this->imageFileName = imageFileName;
	
	this->numberOfInitialClasses = numberOfInitialClasses;


}
Daniela::~Daniela()
{

}


float Daniela::Set_S_Image(ImageType::Pointer inputimage)
{

	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	IteratorType It( inputimage, inputimage->GetRequestedRegion());
	ImageType::SizeType  size;
	size = inputimage->GetRequestedRegion().GetSize();
	M=size[1];	
	N=size[0];

	ImageType::IndexType requestedIndex =inputimage->GetRequestedRegion().GetIndex();


	/*unsigned char **inImg = (unsigned char **) malloc (M * sizeof(unsigned char *));
	for (int k = 0; k < M; k++)
		inImg[k] = (unsigned char *) malloc (N * sizeof(unsigned char));*/

	float mean=0.0;
	for (It.GoToBegin(); !It.IsAtEnd(); ++It)
	{

		ImageType::IndexType idx = It.GetIndex();
		//inImg[idx[1]][idx[0]]= It.Get(); 
		mean  += It.Get();

	}

	//Free array memory
	/*for (int k = 0; k < M; k++)
		free(inImg[k]);
	free(inImg);
	inImg = NULL;*/

	mean =mean/(float)(M*N);
	cout<<	mean<<", ";
	return mean;

}




ImageType::Pointer Daniela::GzetImage(ImageType::Pointer inputimage )
{

	ImageType::Pointer image = ImageType::New();	
	ImageType::RegionType outputRegion;
	ImageType::RegionType::SizeType  size;

	size[0]  = N;
	size[1]  = M;
	size[2] = 1;


	ImageType::RegionType::IndexType outputStart;

	outputStart[0] = 0;
	outputStart[1] = 0;
	outputStart[2] =0;

	outputRegion.SetSize( size );
	outputRegion.SetIndex( outputStart );
	image->SetRegions(outputRegion);
	image->Allocate();
	typedef itk::ImageRegionIterator< ImageType> IteratorType;

	IteratorType It( inputimage, inputimage->GetRequestedRegion());
	IteratorType out(image ,outputRegion);


	int index = 0;			 
	for (It.GoToBegin(),out.GoToBegin(); !It.IsAtEnd(); ++It,++out)
	{

		if (It.Get()==0)
			out.Set(255);
		else
			out.Set(0);

	}
	return image;
}
void Daniela::Read_RGBInputImage()
{
	ImageFileReaderType::Pointer imageReader = ImageFileReaderType::New();
	std::cout << "Reading file: " << imageFileName << std::endl;
	imageReader->SetFileName( imageFileName);
	try
	{
		imageReader ->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught !" << std::endl; 
		std::cout << err << std::endl; 
		
	}

	// Set input images from reader.	
	rgbImage = imageReader->GetOutput();

}


int Daniela::Get_Daniela()
{

	Read_RGBInputImage();
	int useNonContiguousLabels = 256;

	// Instantiate reader for input image.
	
	//if(outputColorSpace == std::string("HSV"))

	std::cerr << "Writing H, S, and V channels of HSV color space." << std::endl; 
	// Convert RGB image to HSV color space
	typedef itk::Accessor::RGBToHSVColorSpacePixelAccessor<unsigned char, float> RGBToHSVColorSpaceAccessorType;
	typedef itk::ImageAdaptor<RGBImageType, RGBToHSVColorSpaceAccessorType> RGBToHSVAdaptorType;
	RGBToHSVAdaptorType::Pointer rgbToHSVAdaptor = RGBToHSVAdaptorType::New();
	rgbToHSVAdaptor->SetImage(rgbImage);

	typedef itk::VectorIndexSelectionCastImageFilter<RGBToHSVAdaptorType, FloatImageType> VectorCastFilterType;
	VectorCastFilterType::Pointer vectorCastFilter = VectorCastFilterType::New();
	vectorCastFilter->SetInput(rgbToHSVAdaptor);
	vectorCastFilter->SetIndex(1); /// get The S image

	vectorCastFilter->Update();

	// Rescale float image to char image.
	typedef itk::RescaleIntensityImageFilter<FloatImageType, CharImageType> FloatRescalerType; 
	FloatRescalerType::Pointer floatRescaler = FloatRescalerType::New(); 
	floatRescaler->SetOutputMinimum(0);
	floatRescaler->SetOutputMaximum(255);
	floatRescaler->SetInput(vectorCastFilter->GetOutput());
	ImageType::Pointer inputImage = ImageType::New();
	inputImage =(floatRescaler->GetOutput());
	floatRescaler->Update();
	CharImageWriterType::Pointer writer = CharImageWriterType::New();
	writer->SetInput(floatRescaler->GetOutput());
	// Write out S component
	vectorCastFilter->SetIndex(1);
	writer->SetFileName(imageFileName + std::string("_HSV_S.png"));
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Problem encountered while writing ";
		std::cerr << " image file : " << imageFileName+ std::string("_HSV_S.png") << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	float mean1 =Set_S_Image(inputImage);


	std::vector<float> userMeans;
	for( int k = 0; k < numberOfInitialClasses; k++ )
	{
		float userProvidedInitialMean = mean1+(pow(-1.0,k)*6*numberOfInitialClasses);
		userMeans.push_back(userProvidedInitialMean);
	}


	// Instantiate the ScalarImageKmeansImageFilter  
	typedef itk::ScalarImageKmeansImageFilter<ImageType > KMeansFilterType;

	KMeansFilterType::Pointer kmeansFilter = KMeansFilterType::New();

	kmeansFilter->SetInput( inputImage);

	// Make the output image intellegable by expanding the range of output image values, if desired

	kmeansFilter->SetUseNonContiguousLabels( useNonContiguousLabels );

	// initialize using the user input means

	for( int k = 0; k < numberOfInitialClasses; k++ )
	{
		kmeansFilter->AddClassWithInitialMean( userMeans[k] );
	}

	// Create and setup a writer

	typedef KMeansFilterType::OutputImageType  OutputImageType;
	kmeansFilter->Update();

	//typedef itk::ImageFileWriter< OutputImageType > WriterType;



	//writer->SetInput( kmeansFilter->GetOutput() );

	//writer->SetFileName( outputImageFileName );

	//// execute the pipeline
	//try
	//{
	//	writer->Update();
	//}
	//catch( itk::ExceptionObject & excp )
	//{
	//	std::cerr << "Problem encountered while writing ";
	//	std::cerr << " image file : " << outputImageFileName << std::endl;
	//	std::cerr << excp << std::endl;
	//	return EXIT_FAILURE;
	//}



	// inspect the means
	KMeansFilterType::ParametersType estimatedMeans = 
		kmeansFilter->GetFinalMeans();

	const unsigned int numberOfClasses = estimatedMeans.Size();

	for ( unsigned int i = 0 ; i < numberOfClasses ; ++i )
	{
		
		
		std::cout << "cluster[" << i << "] ";
		std::cout << "    estimated mean : " << estimatedMeans[i] << std::endl;
	}

	// //the following code extracts the bone part from the clustered image 

	ImageType::Pointer Bone_Image = ImageType::New();
	Bone_Image = GzetImage(kmeansFilter->GetOutput());
	kmeansFilter->Update();
	//writer->SetInput(  Bone_Image );
	//writer->SetFileName(imageFileName + std::string("_segmented.png") );
	//try
	//{
	//	writer->Update();
	//}
	//catch( itk::ExceptionObject & excp )
	//{
	//	std::cerr << "Problem encountered while writing ";
	//	std::cerr << " image file : " << imageFileName + std::string("_segmented.png")<< std::endl;
	//	std::cerr << excp << std::endl;
	//	return EXIT_FAILURE;
	//}


	// the ConnectedComponentImageFilter and RelabelComponentImageFilte give the labelled image from the extracted bone image.

	typedef itk::ConnectedComponentImageFilter< ImageType, ShortImageType > CCFilterType;
	typedef itk::RelabelComponentImageFilter<  ShortImageType,  ImageType > RelabelType;

	CCFilterType::Pointer ccfilter = CCFilterType::New();
	RelabelType::Pointer relabel = RelabelType::New();

	ccfilter->SetInput(	Bone_Image  );
	ccfilter->FullyConnectedOn();

	relabel->SetInput( ccfilter->GetOutput() );
	relabel->SetMinimumObjectSize(1);
	relabel->InPlaceOn();

	try
	{
		relabel->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Relabel: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return EXIT_FAILURE;
	}
	
	bone_label_image = relabel->GetOutput();

	writer->SetInput(   relabel->GetOutput() );
	writer->SetFileName(imageFileName + std::string("_label.tif") );
	std::cout << "Writing out file: " << imageFileName + std::string("_label.tif") << std::endl;
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Problem encountered while writing ";
		std::cerr << " image file : " << imageFileName + std::string("_label.png") << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}


	//// Write out H component
	//CharImageWriterType::Pointer writer = CharImageWriterType::New();
	//writer->SetInput(floatRescaler->GetOutput());
	//writer->SetFileName(argv[1]+ std::string("_HSV_H.png"));
	//writer->Update();
	//   

	//// Write out V component
	//vectorCastFilter->SetIndex(2);
	//writer->SetFileName(argv[1] + std::string("_HSV_V.png"));
	//writer->Update();

	return EXIT_SUCCESS;
}