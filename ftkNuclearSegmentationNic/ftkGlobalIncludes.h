// ############################################################################################################################################################################
#ifndef ftkGlobalIncludes_h_
#define ftkGlobalIncludes_h_
// ############################################################################################################################################################################


// STANDARD INCLUES
#include <vector>
#include <iostream>
#include <map>
#include <limits>

// FTK INCLUDES
#include <ftkObject.h>
#include <ftkImage/ftkImage.h>
#include <ftkCommon/ftkUtils.h>


// ITK INCLUDES
#include <itkImage.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>

// VTK INCLUDES



// #include <time.h>
// #include <iostream>
// #include <vector>
// #include <cmath>
// 
// #include "itkImage.h"
// #include "itkImageFileReader.h"
// #include "itkImageFileWriter.h"
// #include "itkRescaleIntensityImageFilter.h"
// 
// // #include "itkSobelOperator.h"
// // #include "itkNeighborhoodInnerProduct.h"
// #include "itkMinimumMaximumImageCalculator.h"
// 
// // #include "itkComposeImageFilter.h"
// // #include "itkVectorImage.h"
// 
// // #include "itkCannyEdgeDetectionImageFilter.h"
// 
// // #include "itkRecursiveGaussianImageFilter.h"
// #include "itkImageDuplicator.h"
// 
// // Connected components
// // #include "itkConnectedComponentImageFilter.h"


#ifdef _OPENMP
#include "omp.h"
#endif

// ############################################################################################################################################################################

namespace ftk{
	/**
	* namespace corresponding to the segmentation algorithm 
	* implemented by nicolas
	*/
 	namespace nucSecNic{
// 		template <typename T1, typename T2>
// 		int writeImage(typename T1::Pointer im, const char* filename)
// 		{
// 			typedef itk::RescaleIntensityImageFilter< T1, T2 > RescaleFilterType;
// 			typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
// 			rescaleFilter->SetInput( im );
// 			rescaleFilter->SetOutputMaximum( 255/*std::numeric_limits< T2::PixelType >::max()*/ );
// 			rescaleFilter->SetOutputMinimum( 0 ); 
// 			try{
// 				rescaleFilter->Update();}
// 			catch( itk::ExceptionObject & err ){
// 				std::cerr << "ExceptionObject caught!" << std::endl;
// 				std::cerr << err << std::endl;
// 				return 1;}
// 
// 			// Set Up the Writer and Write the result
// 			typedef itk::ImageFileWriter < T2 >  WriterType;
// 			typename WriterType::Pointer writer = WriterType::New();
// 			writer->SetInput( rescaleFilter->GetOutput() );
// 			writer->SetFileName( filename );
// 			try{
// 				writer->Update();}
// 			catch( itk::ExceptionObject & err ){
// 				std::cerr << "ExceptionObject caught!" << std::endl;
// 				std::cerr << err << std::endl;
// 				return 1;}
// 			return 0;
// 		};

		template <typename T1>
		int writeImageNoScaling(typename T1::Pointer im, const char* filename)
		{
			// Set Up the Writer and Write the result
			typedef itk::ImageFileWriter < T1 >  WriterType;
			typename WriterType::Pointer writer = WriterType::New();
			writer->SetInput( im );
			writer->SetFileName( filename );
			try{
				writer->Update();}
			catch( itk::ExceptionObject & err ){
				std::cerr << "ExceptionObject caught!" << std::endl;
				std::cerr << err << std::endl;
				return 1;}
			return 0;
		};

// 		template <typename T1, typename T2>
// 		typename T1::Pointer readImage( const char* filename )
// 		{
// 			// Set Up the Reader
// 			std::cout << std::endl << std::endl << "Reading ... " << filename;
// 			typedef itk::ImageFileReader < T1 >  ReaderType;
// 			typename ReaderType::Pointer reader = ReaderType::New();
// 			reader->SetFileName( filename );
// 			try{
// 				reader->Update();}
// 			catch( itk::ExceptionObject & err ){
// 				std::cerr << "ExceptionObject caught!" << std::endl;
// 				std::cerr << err << std::endl;
// 			}
// 			typename T1::Pointer inputImage = reader->GetOutput(); // Pass the readed image to a "regular image"
// 
// 			// RescaleIntensityImageFilter
// 			typedef itk::RescaleIntensityImageFilter< T1, T2 > RescaleFilterType;
// 			typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
// 			rescaleFilter->SetInput( inputImage );
// 			rescaleFilter->SetOutputMaximum( 1.0 );
// 			rescaleFilter->SetOutputMinimum( 0.0 ); 
// 			try{
// 				rescaleFilter->Update();}
// 			catch( itk::ExceptionObject & err ){
// 				std::cerr << "ExceptionObject caught!" << std::endl;
// 				std::cerr << err << std::endl;
// 			}
// 			inputImage = rescaleFilter->GetOutput();
// 			return inputImage;
// 		};
	};
};


void stopProgram( void );


#endif
// ############################################################################################################################################################################