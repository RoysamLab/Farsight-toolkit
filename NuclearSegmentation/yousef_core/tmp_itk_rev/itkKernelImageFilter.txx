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
  Module:    $RCSfile: itkKernelImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008-08-12 12:55:28 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkKernelImageFilter_txx
#define __itkKernelImageFilter_txx

#include "itkKernelImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkFlatStructuringElement.h"

namespace itk {

template <class TInputImage, class TOutputImage, class TKernel>
KernelImageFilter<TInputImage, TOutputImage, TKernel>
::KernelImageFilter()
{
  this->SetRadius( 1 );
}


template <class TInputImage, class TOutputImage, class TKernel>
void
KernelImageFilter<TInputImage, TOutputImage, TKernel>
::SetRadius( const RadiusType & radius )
{
  // Try to use a FlatStructuringElement if possible, because it is
  // much efficient with van Herk / Gil Werman filters
  // check that the type is exactly FlatKernelType, to be sure to
  // not replace m_Kernel data with the data of a FlatKernelType
  // if KernelType is a subclass of FlatKernelType.
  if( typeid(KernelType) == typeid(FlatKernelType) )
    {
    // SetKernel() must be called, because it can be overloaded in a subclass
    // - MovingHistogramImageFilterBase for example.
    FlatKernelType flatKernel = FlatKernelType::Box( radius );
    KernelType * kernel = reinterpret_cast< KernelType* >( & flatKernel );
    this->SetKernel( *kernel );
    }
  else
    {
    // Superclass::SetRadius() is called in SetKenel() - no need to call it here
    // Superclass::SetRadius( radius );
    KernelType kernel;
    kernel.SetRadius( radius );
    for( typename KernelType::Iterator kit=kernel.Begin(); kit != kernel.End(); kit++ )
      {
      *kit = 1;
      }
    this->SetKernel( kernel );
    }
}


template <class TInputImage, class TOutputImage, class TKernel>
void
KernelImageFilter<TInputImage, TOutputImage, TKernel>
::SetKernel( const KernelType & kernel )
{
  if( m_Kernel != kernel )
    {
    m_Kernel = kernel;
    this->Modified();
    }
  // set the radius of the super class to be the same than the kernel one
  Superclass::SetRadius( kernel.GetRadius() );
}


template <class TInputImage, class TOutputImage, class TKernel>
void
KernelImageFilter<TInputImage, TOutputImage, TKernel>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Kernel: " << m_Kernel << std::endl;
}

}


#endif
