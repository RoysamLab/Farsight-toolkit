#ifndef GENERALIO_H
#define GENERALIO_H

#include <time.h>
#include<iostream>
#include<vector>
#include <cmath>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#ifdef _OPENMP
#include "omp.h"
#endif

template <typename T1>
typename T1::Pointer readImage( const char* filename )
	{
		typedef itk::ImageFileReader < T1 >  ReaderType;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( filename );
		try{
			reader->Update();}
		catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught in reading " <<filename<<" !" << std::endl;
			std::cerr << err << std::endl;
		}
		typename T1::Pointer inputImage = reader->GetOutput(); 

		return inputImage;
	};

template <typename T1, typename T2>
int writeImage(typename T1::Pointer im, const char* filename)
{
		typedef itk::ImageFileWriter < T2 >  WriterType;
		typename WriterType::Pointer writer = WriterType::New();

		typedef itk::RescaleIntensityImageFilter< T1, T2 > RescaleFilterType;
		typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
		rescaleFilter->SetInput( im );
		rescaleFilter->SetOutputMaximum( 255 );
		rescaleFilter->SetOutputMinimum( 0 ); 
		try
		{
			rescaleFilter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught in rescaling "<<filename<< " !" << std::endl;
			std::cerr << err << std::endl;
			return 1;
		}
		writer->SetInput( rescaleFilter->GetOutput () );
		writer->SetFileName( filename );
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught in writting "<<filename<< " !" << std::endl;
			std::cerr << err << std::endl;
			return 1;
		}

		return 0;
};

template <typename T1 >
int writeImage(typename T1::Pointer im, const char* filename)
{
		typedef itk::ImageFileWriter < T1 >  WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput( im );
		writer->SetFileName( filename );
		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught in writting "<<filename<< " !" << std::endl;
			std::cerr << err << std::endl;
			return 1;
		}

		return 0;
};

#endif