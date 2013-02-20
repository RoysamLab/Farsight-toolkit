// ############################################################################################################################################################################
#ifndef FTKVOTINGGLOBAL_H
#define FTKVOTINGGLOBAL_H

#include <time.h>
#include<iostream>
#include<vector>
#include <cmath>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkSobelOperator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "itkComposeImageFilter.h"
#include "itkVectorImage.h"

#include "itkCannyEdgeDetectionImageFilter.h"

#include "itkRecursiveGaussianImageFilter.h"
#include "itkImageDuplicator.h"


#ifdef _OPENMP
#include "omp.h"
#endif

// ############################################################################################################################################################################

namespace nftkVot{

	// Input Image Type
	typedef double InputPixelType;
	const unsigned int Dimension = 2;
	typedef itk::Image< InputPixelType, Dimension > InputImageType;

	const unsigned int Dimension_3 = 3;
	typedef itk::Image< InputPixelType, Dimension_3 > InputImageType_3D;

	typedef unsigned short InputPixelType_16;
	typedef itk::Image< InputPixelType_16, Dimension_3 > InputImageType_3D_16;

	typedef unsigned char InputPixelType_8;
	typedef itk::Image< InputPixelType_8, Dimension_3 > InputImageType_3D_8;


	template <typename T1, typename T2>
	int writeImage(typename T1::Pointer im, const char* filename)
	{
		typedef itk::RescaleIntensityImageFilter< T1, T2 > RescaleFilterType;
		typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetInput( im );
		rescaleFilter->SetOutputMaximum( 255 );
		rescaleFilter->SetOutputMinimum( 0 ); 
		try{
			rescaleFilter->Update();}
		catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught!" << std::endl;
			std::cerr << err << std::endl;
			return 1;}

		// Set Up the Writer and Write the result
		typedef itk::ImageFileWriter < T2 >  WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput( rescaleFilter->GetOutput() );
		writer->SetFileName( filename );
		try{
			writer->Update();}
		catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught!" << std::endl;
			std::cerr << err << std::endl;
			return 1;}
		return 0;
	};

	template <typename T1>
	int writeImage_mhdDouble(typename T1::Pointer im, const char* filename)
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

	template <typename T1>
	int writeImage_mhd256(typename T1::Pointer im, const char* filename)
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


	template <typename T1, typename T2>
	typename T1::Pointer readImage( const char* filename )
	{
		// Set Up the Reader
		std::cout << std::endl << std::endl << "Reading ... " << filename;
		typedef itk::ImageFileReader < T1 >  ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename );
		try{
			reader->Update();}
		catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught!" << std::endl;
			std::cerr << err << std::endl;
		}
		typename T1::Pointer inputImage = reader->GetOutput(); // Pass the readed image to a "regular image"

		// RescaleIntensityImageFilter
		typedef itk::RescaleIntensityImageFilter< T1, T2 > RescaleFilterType;
		typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetInput( inputImage );
		rescaleFilter->SetOutputMaximum( 1.0 );
		rescaleFilter->SetOutputMinimum( 0.0 ); 
		try{
			rescaleFilter->Update();}
		catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught!" << std::endl;
			std::cerr << err << std::endl;
		}
		inputImage = rescaleFilter->GetOutput();
		return inputImage;
	};

	template <typename T1, typename T2>
	typename T2::Pointer readImage_3D( const char* filename )
	{
		// Set Up the Reader
		std::cout << std::endl << std::endl << "Reading ... " << filename;
		typedef itk::ImageFileReader < T1 >  ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename );
		try{
			reader->Update();}
		catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught!" << std::endl;
			std::cerr << err << std::endl;
		}


		int nxx = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
		int nyx = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
		int nzx = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[2];
		int npix = nxx*nyx*nzx;

		std::cout << std::endl << "The Size of the Image is: "<<nxx<<" "<<nyx<<" "<<nzx;

		// This is temporary
		// RescaleIntensityImageFilter
		typedef itk::RescaleIntensityImageFilter< T1, T2 > RescaleFilterType;
		typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetInput( reader->GetOutput() );
		rescaleFilter->SetOutputMaximum( 1.0 );
		rescaleFilter->SetOutputMinimum( 0.0 ); 
		try{
			rescaleFilter->Update();}
		catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught!" << std::endl;
			std::cerr << err << std::endl;
		}
		return rescaleFilter->GetOutput();
	};

	template <typename T1, typename T2>
	int writeImage_3D(typename T1::Pointer im, const char* filename)
	{
		typedef itk::RescaleIntensityImageFilter< T1, T2 > RescaleFilterType;
		typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetInput( im );
		//rescaleFilter->SetOutputMaximum( 255 ); //65,535
		rescaleFilter->SetOutputMaximum( 65535 ); //65,535
		rescaleFilter->SetOutputMinimum( 0 ); 
		try{
			rescaleFilter->Update();}
		catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught!" << std::endl;
			std::cerr << err << std::endl;
			return 1;}

		// Set Up the Writer and Write the result
		typedef itk::ImageFileWriter < T2 >  WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput( rescaleFilter->GetOutput() );
		writer->SetFileName( filename );
		try{
			writer->Update();}
		catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught!" << std::endl;
			std::cerr << err << std::endl;
			return 1;}
		return 0;
	};

	double diffclock(clock_t clock1,clock_t clock2);

	int round_double( double x );

	int round_double2( double x );

	void stopProgram( void );

};

#endif
// ############################################################################################################################################################################
