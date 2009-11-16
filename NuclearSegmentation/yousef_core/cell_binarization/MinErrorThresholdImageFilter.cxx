/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*=========================================================================
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: MinErrorThresholdImageFilter.cxx,v $
  Language:  C++
  Date:      $Date: 2006/05/21 13:00:23 $
  Version:   $Revision: 1.10 $

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
//    OUTPUTS: {MinErrorThresholdImageFilterOutput.png}
//    255 0
//  Software Guide : EndCommandLineArgs

// Software Guide : BeginLatex
//
// This example illustrates how to use the \doxygen{MinErrorThresholdImageFilter}.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkMinErrorThresholdImageFilter.h"
// Software Guide : EndCodeSnippet

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImageFile outputImageFile ";  
    std::cerr << " insideValue    outsideValue   "  << std::endl;  
    return EXIT_FAILURE;
    }
  
  //  Software Guide : BeginLatex
  //
  //  The next step is to decide which pixel types to use for the input and output
  //  images.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef  unsigned char  InputPixelType;
  typedef  unsigned char  OutputPixelType;
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  The input and output image types are now defined using their respective
  //  pixel types and dimensions.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::Image< InputPixelType,  2 >   InputImageType;
  typedef itk::Image< OutputPixelType, 2 >   OutputImageType;
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  The filter type can be instantiated using the input and output image
  //  types defined above.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::MinErrorThresholdImageFilter<
               InputImageType, OutputImageType >  FilterType;
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  An \doxygen{ImageFileReader} class is also instantiated in order to read
  //  image data from a file. (See Section \ref{sec:IO} on page 
  //  \pageref{sec:IO} for more information about reading
  //  and writing data.) 
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //  
  // An \doxygen{ImageFileWriter} is instantiated in order to write the output
  // image to a file.
  //
  //  Software Guide : EndLatex 


  // Software Guide : BeginCodeSnippet
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  Both the filter and the reader are created by invoking their \code{New()}
  //  methods and assigning the result to \doxygen{SmartPointer}s.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  ReaderType::Pointer reader = ReaderType::New();
  FilterType::Pointer filter = FilterType::New();
  // Software Guide : EndCodeSnippet

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  reader->SetFileName( argv[1] );


  //  Software Guide : BeginLatex
  //  
  //  The image obtained with the reader is passed as input to the
  //  MinErrorThresholdImageFilter.
  //
  //  \index{itk::MinError\-Threshold\-Image\-Filter!SetInput()}
  //  \index{itk::FileImageReader!GetOutput()}
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  filter->SetInput( reader->GetOutput() );
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //  
  //  The method \code{SetOutsideValue()} defines the intensity value to be
  //  assigned to those pixels whose intensities are outside the range defined
  //  by the lower and upper thresholds. The method \code{SetInsideValue()}
  //  defines the intensity value to be assigned to pixels with intensities
  //  falling inside the threshold range.
  //  
  //  \index{itk::MinError\-Threshold\-Image\-Filter!SetOutsideValue()}
  //  \index{itk::MinError\-Threshold\-Image\-Filter!SetInsideValue()}
  //  \index{SetOutsideValue()!itk::MinError\-Threshold\-Image\-Filter}
  //  \index{SetInsideValue()!itk::MinError\-Threshold\-Image\-Filter}
  //
  //  Software Guide : EndLatex 

  const OutputPixelType outsideValue = atoi( argv[3] );
  const OutputPixelType insideValue  = atoi( argv[4] );

  // Software Guide : BeginCodeSnippet
  filter->SetOutsideValue( outsideValue );
  filter->SetInsideValue(  insideValue  );
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //  
  //  The method \code{SetNumberOfHistogramBins()} defines the number of bins
  //  to be used for computing the histogram. This histogram will be used
  //  internally in order to compute the MinError threshold.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  filter->SetNumberOfHistogramBins( 128 );
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //  
  //  The execution of the filter is triggered by invoking the \code{Update()}
  //  method.   If the filter's output has been passed as input to subsequent
  //  filters, the \code{Update()} call on any posterior filters in the
  //  pipeline will indirectly trigger the update of this filter.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  filter->Update();
  // Software Guide : EndCodeSnippet



  //  Software Guide : BeginLatex
  //  
  //  We print out here the Threshold value that was computed internally by the
  //  filter. For this we invoke the \code{GetThreshold} method.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  int threshold = filter->GetThreshold();
  std::cout << "Threshold = " << threshold << std::endl;
  // Software Guide : EndCodeSnippet



  //  Software Guide : BeginLatex
  //  
  // \begin{figure}
  // \center
  // \includegraphics[width=0.44\textwidth]{BrainProtonDensitySlice.eps}
  // \includegraphics[width=0.44\textwidth]{MinErrorThresholdImageFilterOutput.eps}
  // \itkcaption[MinErrorThresholdImageFilter output]{Effect of the MinErrorThresholdImageFilter on a slice from a MRI
  // proton density image  of the brain.}
  // \label{fig:MinErrorThresholdImageFilterInputOutput}
  // \end{figure}
  //
  //  Figure \ref{fig:MinErrorThresholdImageFilterInputOutput} illustrates the
  //  effect of this filter on a MRI proton density image of the brain. This
  //  figure shows the limitations of this filter for performing segmentation
  //  by itself. These limitations are particularly noticeable in noisy images
  //  and in images lacking spatial uniformity as is the case with MRI due to
  //  field bias.
  //
  //  \relatedClasses
  //  \begin{itemize}
  //  \item \doxygen{ThresholdImageFilter}
  //  \end{itemize}
  //
  //  Software Guide : EndLatex 

  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS;
}

