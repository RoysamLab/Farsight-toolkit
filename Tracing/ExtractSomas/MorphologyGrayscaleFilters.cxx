/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: MathematicalMorphologyGrayscaleFilters.cxx,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

//  Software Guide : BeginCommandLineArgs
//    INPUTS: {BrainProtonDensitySlice.png}
//    OUTPUTS: {MathematicalMorphologyGrayscaleErosionOutput.png}
//    OUTPUTS: {MathematicalMorphologyGrayscaleDilationOutput.png}
//    150 180
//  Software Guide : EndCommandLineArgs


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 


int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImageFile  ";
    std::cerr << " outputImageFileErosion  outputImageFileDilation radius" << std::endl;
    return EXIT_FAILURE;
    }



  const unsigned int Dimension = 3;
  
  typedef unsigned char   InputPixelType;
  typedef unsigned char   OutputPixelType;

  typedef itk::Image< InputPixelType,  Dimension >   InputImageType;
  typedef itk::Image< OutputPixelType, Dimension >   OutputImageType;

  typedef itk::ImageFileReader< InputImageType  >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;



  typedef itk::BinaryBallStructuringElement< 
                      InputPixelType,
                      Dimension  >             StructuringElementType;



  typedef itk::GrayscaleErodeImageFilter<
                            InputImageType, 
                            OutputImageType,
                            StructuringElementType >  ErodeFilterType;

  typedef itk::GrayscaleDilateImageFilter<
                            InputImageType, 
                            OutputImageType, 
                            StructuringElementType >  DilateFilterType;


  // Creation of Reader and Writer filters
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writerDilation = WriterType::New();
  WriterType::Pointer writerErosion  = WriterType::New();



  ErodeFilterType::Pointer  grayscaleErode  = ErodeFilterType::New();
  DilateFilterType::Pointer grayscaleDilate = DilateFilterType::New();


  StructuringElementType  structuringElement;

  structuringElement.SetRadius( atoi( argv[4] ) );  // 3x3 structuring element

  structuringElement.CreateStructuringElement();

  grayscaleErode->SetKernel(  structuringElement );
  grayscaleDilate->SetKernel( structuringElement );


  reader->SetFileName( argv[1] );
 
  writerErosion->SetFileName(  argv[2] );
  writerDilation->SetFileName( argv[3] );
  



  grayscaleErode->SetInput(  reader->GetOutput() );
//  grayscaleDilate->SetInput( reader->GetOutput() );
  grayscaleDilate->SetInput( grayscaleErode->GetOutput() );




  writerDilation->SetInput( grayscaleDilate->GetOutput() );
  writerDilation->Update();

  writerErosion->SetInput( grayscaleErode->GetOutput() );
  writerErosion->Update();



  return EXIT_SUCCESS;
}

