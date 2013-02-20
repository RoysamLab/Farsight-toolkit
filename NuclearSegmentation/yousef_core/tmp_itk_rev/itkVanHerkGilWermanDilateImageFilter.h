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
  Module:    $RCSfile: itkVanHerkGilWermanDilateImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-02-24 19:03:15 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkVanHerkGilWermanDilateImageFilter_h
#define __itkVanHerkGilWermanDilateImageFilter_h

#include "itkVanHerkGilWermanErodeDilateImageFilter.h"
#include "vnl/vnl_math.h"

namespace itk {
template <class TPixel>
class MaxFunctor
{
public:
  MaxFunctor(){}
  ~MaxFunctor(){}
  inline TPixel operator()(const TPixel &A, const TPixel &B) const
    {
    return vnl_math_max(A, B);
    }
};


template<class TImage, class TKernel>
class  ITK_EXPORT VanHerkGilWermanDilateImageFilter :
    public VanHerkGilWermanErodeDilateImageFilter<TImage, TKernel, MaxFunctor<typename TImage::PixelType> >

{
public:
  typedef VanHerkGilWermanDilateImageFilter Self;
  typedef VanHerkGilWermanErodeDilateImageFilter<TImage, TKernel, MaxFunctor<typename TImage::PixelType> > Superclass;

  /** Runtime information support. */
  itkTypeMacro(VanHerkGilWermanDilateImageFilter, 
               VanHerkGilWermanErodeDilateImageFilter);

  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  virtual ~VanHerkGilWermanDilateImageFilter() {}

protected:

  typedef typename TImage::PixelType  PixelType;

  VanHerkGilWermanDilateImageFilter()
    {
    this->m_Boundary = NumericTraits< PixelType >::NonpositiveMin();
    }

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    os << indent << "VanHerkGilWerman erosion: " << std::endl;
    }

private:
  
  VanHerkGilWermanDilateImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


} // namespace itk

#endif
