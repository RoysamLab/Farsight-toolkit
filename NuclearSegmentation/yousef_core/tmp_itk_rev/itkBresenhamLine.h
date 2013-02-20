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
  Module:    $RCSfile: itkBresenhamLine.h,v $
  Language:  C++
  Date:      $Date: 2008-08-06 16:49:13 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkBresenhamLine_h
#define __itkBresenhamLine_h

#include "itkVector.h"
#include "itkOffset.h"
#include "itkIndex.h"
#include <vector>

namespace itk {

/* a simple class that will return an array of indexes that are
* offsets along the line. The line will be described by a vector and a
* length */

template <unsigned int VDimension>
class ITK_EXPORT BresenhamLine
{
public:
  typedef BresenhamLine                      Self;
  // This defines the line direction
  typedef Vector<float, VDimension>          LType;
  typedef Offset<VDimension>                 OffsetType;
  typedef Index<VDimension>                  IndexType;
  typedef std::vector<OffsetType>            OffsetArray;

  typedef typename IndexType::IndexValueType IndexValueType;

  // constructurs
  BresenhamLine(){}
  ~BresenhamLine(){}

  OffsetArray buildLine(LType Direction, unsigned int length);

};


}


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBresenhamLine.txx"
#endif


#endif
