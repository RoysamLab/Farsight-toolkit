// ############################################################################################################################################################################
#ifndef FTKVOTINGGLOBAL_H
#define FTKVOTINGGLOBAL_H

#include <time.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

// ############################################################################################################################################################################

namespace nftkVotingGlobal{

	// Input Image Type
	typedef double InputPixelType;
	const unsigned int Dimension = 2;
	typedef itk::Image< InputPixelType, Dimension > InputImageType;


template <typename T1, typename T2>
int writeImage(typename T1::Pointer im, const char* filename)
{
	typedef itk::RescaleIntensityImageFilter< T1, T2 > RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput( im );
	rescaleFilter->SetOutputMaximum( 255 );
	rescaleFilter->SetOutputMinimum( 0 ); 
	rescaleFilter->Update();

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
	printf("Reading %s ... ",filename);
	typedef itk::ImageFileReader < T1 >  ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( filename );
	try{
		reader->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	T1::Pointer inputImage = reader->GetOutput(); // Pass the readed image to a "regular image"

	// RescaleIntensityImageFilter
	typedef itk::RescaleIntensityImageFilter< T1, T2 > RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput( inputImage );
	rescaleFilter->SetOutputMaximum( 1.0 );
	rescaleFilter->SetOutputMinimum( 0.0 ); 
	rescaleFilter->Update();
	inputImage = rescaleFilter->GetOutput();
	return inputImage;
};

double diffclock(clock_t clock1,clock_t clock2);

int round_double( double x );

int round_double2( double x );

};

#endif
// ############################################################################################################################################################################