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
  Module:    $RCSfile: itkLabelBorderImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2007/10/05 10:31:58 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkLabelBorderImageFilter_h
#define __itkLabelBorderImageFilter_h

#include "itkInPlaceImageFilter.h"
#include "itkImage.h"
#include "itkConceptChecking.h"
#include <vector>
#include <map>
#include "itkProgressReporter.h"
#include "itkBarrier.h"

namespace itk
{

/**
 * \class LabelBorderImageFilter
 * \brief Give the pixels on the border of an object.
 *
 * LabelBorderImageFilter labels the pixels on the borders
 * of the objects in a binary image.
 *
 * \sa ImageToImageFilter 
 */

template <class TInputImage, class TOutputImage>
class ITK_EXPORT LabelBorderImageFilter : 
    public InPlaceImageFilter< TInputImage, TOutputImage > 
{
public:
  /**
   * Standard "Self" & Superclass typedef.
   */
  typedef LabelBorderImageFilter                   Self;
  typedef InPlaceImageFilter< TInputImage, TOutputImage > Superclass;

  /**
   * Types from the Superclass
   */
  typedef typename Superclass::InputImagePointer InputImagePointer;

  /**
   * Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same.
   */
  typedef typename TOutputImage::PixelType         OutputPixelType;
  typedef typename TOutputImage::InternalPixelType OutputInternalPixelType;
  typedef typename TInputImage::PixelType          InputPixelType;
  typedef typename TInputImage::InternalPixelType  InputInternalPixelType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  
  /**
   * Image typedef support
   */
  typedef TInputImage                       InputImageType;
  typedef typename TInputImage::IndexType   IndexType;
  typedef typename TInputImage::SizeType    SizeType;
  typedef typename TInputImage::OffsetType  OffsetType;
  typedef typename TInputImage::PixelType   InputImagePixelType;

  typedef TOutputImage                      OutputImageType;
  typedef typename TOutputImage::RegionType RegionType;
  typedef typename TOutputImage::IndexType  OutputIndexType;
  typedef typename TOutputImage::SizeType   OutputSizeType;
  typedef typename TOutputImage::OffsetType OutputOffsetType;
  typedef typename TOutputImage::PixelType  OutputImagePixelType;

  typedef std::list<IndexType>              ListType;

  /** 
   * Smart pointer typedef support 
   */
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /**
   * Run-time type information (and related methods)
   */
  itkTypeMacro(LabelBorderImageFilter, ImageToImageFilter);
  
  /**
   * Method for creation through the object factory.
   */
  itkNewMacro(Self);

  /**
   * Set/Get whether the connected components are defined strictly by
   * face connectivity or by face+edge+vertex connectivity.  Default is
   * FullyConnectedOff.  For objects that are 1 pixel wide, use
   * FullyConnectedOn.
   */
  itkSetMacro(FullyConnected, bool);
  itkGetConstReferenceMacro(FullyConnected, bool);
  itkBooleanMacro(FullyConnected);
  
  // Concept checking -- input and output dimensions must be the same
  itkConceptMacro(SameDimension,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),
       itkGetStaticConstMacro(OutputImageDimension)>));

  /**
   */
  itkSetMacro(BackgroundValue, OutputImagePixelType);
  itkGetMacro(BackgroundValue, OutputImagePixelType);

protected:
  LabelBorderImageFilter() 
    {
    m_FullyConnected = false;
    m_BackgroundValue = NumericTraits< OutputImagePixelType >::Zero;
    }
  virtual ~LabelBorderImageFilter() {}
  LabelBorderImageFilter(const Self&) {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /**
   * Standard pipeline methods.
   */
  void BeforeThreadedGenerateData ();
  void AfterThreadedGenerateData ();
  void ThreadedGenerateData (const RegionType& outputRegionForThread, int threadId) ;

  /** LabelBorderImageFilter needs the entire input. Therefore
   * it must provide an implementation GenerateInputRequestedRegion().
   * \sa ProcessObject::GenerateInputRequestedRegion(). */
  void GenerateInputRequestedRegion();

  /** LabelBorderImageFilter will produce all of the output.
   * Therefore it must provide an implementation of
   * EnlargeOutputRequestedRegion().
   * \sa ProcessObject::EnlargeOutputRequestedRegion() */
  void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));

  bool m_FullyConnected;
  
private:
  OutputImagePixelType m_BackgroundValue;

  // some additional types
  typedef typename TOutputImage::RegionType::SizeType OutSizeType;

  // types to support the run length encoding of lines
  class runLength
    {
    public:
    // run length information - may be a more type safe way of doing this
    long int length;
    typename InputImageType::IndexType where; // Index of the start of the run
    InputImagePixelType label;
    };

  typedef std::vector<runLength> lineEncoding;

  // the map storing lines
  typedef std::vector<lineEncoding> LineMapType;
  
  typedef std::vector<long> OffsetVec;

  // the types to support union-find operations
  typedef std::vector<unsigned long int> UnionFindType;

  bool CheckNeighbors(const OutputIndexType &A, 
                      const OutputIndexType &B);

  void CompareLines(lineEncoding &current, const lineEncoding &Neighbour);


  void SetupLineOffsets(OffsetVec &LineOffsets);

  void Wait()
    {
    long nbOfThreads = this->GetNumberOfThreads();
    if( itk::MultiThreader::GetGlobalMaximumNumberOfThreads() != 0 )
      {
      nbOfThreads = std::min( this->GetNumberOfThreads(), itk::MultiThreader::GetGlobalMaximumNumberOfThreads() );
      }
    if( nbOfThreads > 1 )
      {
      m_Barrier->Wait();
      }
    }

  typename std::vector< long > m_NumberOfLabels;
  typename std::vector< long > m_FirstLineIdToJoin;
  typename Barrier::Pointer m_Barrier;
  LineMapType m_LineMap;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabelBorderImageFilter.txx"
#endif

#endif

