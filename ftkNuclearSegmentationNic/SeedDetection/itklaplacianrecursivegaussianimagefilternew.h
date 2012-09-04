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
  Module:    $RCSfile: itkLaplacianRecursiveGaussianImageFilterNew.h,v $
  Language:  C++
  Date:      $Date: 2007/09/27 11:36:40 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLaplacianRecursiveGaussianImageFilterNew_h
#define __itkLaplacianRecursiveGaussianImageFilterNew_h

#include "itkRecursiveGaussianImageFilter.h"
#include "itkNthElementImageAdaptor.h"
#include "itkImage.h"
#include "itkPixelTraits.h"
#include "itkCommand.h"


namespace itk
{

/** \class LaplacianRecursiveGaussianImageFilterNew
 * \brief Computes the Laplacian of an image by convolution with the second derivative of a Gaussian.
 * 
 * This filter is implemented using the recursive gaussian filters
 *
 * 
 * \ingroup GradientFilters   
 * \ingroup Singlethreaded
 */
// NOTE that the ITK_TYPENAME macro has to be used here in lieu 
// of "typename" because VC++ doesn't like the typename keyword 
// on the defaults of template parameters
template <typename TInputImage, 
          typename TOutputImage= TInputImage >
class ITK_EXPORT LaplacianRecursiveGaussianImageFilterNew:
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef LaplacianRecursiveGaussianImageFilterNew  Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                   Pointer;
  typedef SmartPointer<const Self>        ConstPointer;
  
  
  /** Pixel Type of the input image */
  typedef TInputImage                                    InputImageType;
  typedef typename InputImageType::PixelType             PixelType;


  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  typedef typename NumericTraits<PixelType>::RealType      RealType;

  /** Define the image type for internal computations 
      RealType is usually 'double' in NumericTraits. 
      Here we prefer float in order to save memory.  */

  typedef float                                            InternalRealType;
  typedef Image<InternalRealType, 
                itkGetStaticConstMacro(ImageDimension) >   RealImageType;

  /**  Smoothing filter type */
  typedef RecursiveGaussianImageFilter<
    RealImageType,
    RealImageType
    >    GaussianFilterType;

  /**  Derivative filter type, it will be the first in the pipeline  */
  typedef RecursiveGaussianImageFilter<
    InputImageType,
    RealImageType
    >    DerivativeFilterType;

  /**  Pointer to a gaussian filter.  */
  typedef typename GaussianFilterType::Pointer     GaussianFilterPointer;

  /**  Pointer to a derivative filter.  */
  typedef typename DerivativeFilterType::Pointer   DerivativeFilterPointer;

  /**  Pointer to the Output Image */
  typedef typename TOutputImage::Pointer           OutputImagePointer;

  /** Type of the output Image */
  typedef TOutputImage      OutputImageType;
  typedef typename          OutputImageType::PixelType      OutputPixelType;

  /**  Auxiliary image for holding the values of the squared gradient components */
  typedef Image< InternalRealType, 
                 itkGetStaticConstMacro(ImageDimension) >      CumulativeImageType;
  typedef typename CumulativeImageType::Pointer                CumulativeImagePointer;

  /**  Command for observing progress of internal pipeline filters */
  typedef          MemberCommand< Self >     CommandType;  
  typedef typename CommandType::Pointer      CommandPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(LaplacianRecursiveGaussianImageFilterNew, 
               ImageToImageFilter);

  /** Set Sigma value. Sigma is measured in the units of image spacing. */
  void SetSigma( RealType sigma );

  /** by Yousef: use this for the anisotropic case */
  void SetSigmaAnis(RealType* sigma);
  

  /** Define which normalization factor will be used for the Gaussian */
  void SetNormalizeAcrossScale( bool normalizeInScaleSpace );
  itkGetMacro( NormalizeAcrossScale, bool );

  /** LaplacianRecursiveGaussianImageFilterNew needs all of the input to produce an
   * output. Therefore, LaplacianRecursiveGaussianImageFilterNew needs to provide
   * an implementation for GenerateInputRequestedRegion in order to inform
   * the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);

protected:
  LaplacianRecursiveGaussianImageFilterNew();
  virtual ~LaplacianRecursiveGaussianImageFilterNew() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion(DataObject *output);

  /** Compute progress by weighting the contributions of the internal filters */
  void ReportProgress(const Object * object, const EventObject & event );

private:
  LaplacianRecursiveGaussianImageFilterNew(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  GaussianFilterPointer         m_SmoothingFilters[ImageDimension-1];
  DerivativeFilterPointer       m_DerivativeFilter;

  CumulativeImagePointer        m_CumulativeImage;

  CommandPointer                m_ProgressCommand;
  float                         m_Progress;

  /** Normalize the image across scale space */
  bool m_NormalizeAcrossScale; 


};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itklaplacianrecursivegaussianimagefilternew.txx"
#endif

#endif




