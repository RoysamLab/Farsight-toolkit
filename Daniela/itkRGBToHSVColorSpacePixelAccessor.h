/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRGBToHSVColorSpacePixelAccessor.h,v $
  Language:  C++
  Date:      $Date: 2009-03-03 15:08:46 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRGBToHSVColorSpacePixelAccessor_h
#define __itkRGBToHSVColorSpacePixelAccessor_h


#include "itkRGBPixel.h"
#include "itkVector.h"
#include "vnl/vnl_math.h"
#include "itkNumericTraits.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageAdaptor.h"
#include "itkRedPixelAccessor.h"
#include "itkGreenPixelAccessor.h"
#include "itkBluePixelAccessor.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
namespace Accessor
{
/**
 * \class RGBToHSVColorSpacePixelAccessor
 * \brief Give access to a RGBPixel as if it were in HSI Color Space as a Vector type.
 *
 * This class is intended to be used as parameter of 
 * an ImageAdaptor to make an RGBPixel image appear as being
 * an image of Vector pixel type in HSI Color Space.
 *
 * \sa ImageAdaptor
 * \ingroup ImageAdaptors
 *
 */

template <class TInput, class TOutput>
class ITK_EXPORT RGBToHSVColorSpacePixelAccessor
{
public:
  /** Standard class typedefs. */
  typedef   RGBToHSVColorSpacePixelAccessor        Self;

 /** External typedef. It defines the external aspect
   * that this class will exhibit */
  typedef  Vector<TOutput,3>     ExternalType;

  /** Internal typedef. It defines the internal real
   * representation of data */
  typedef   RGBPixel<TInput>    InternalType;

  /** Write access to the RGBToHSVColorSpace component */
  inline void Set( InternalType & output, const ExternalType & input ) const
    { 
    // Normalize RGB values.
    double r = (double)input[0] / (double)NumericTraits<TInput>::max();
    double g = (double)input[1] / (double)NumericTraits<TInput>::max();
    double b = (double)input[2] / (double)NumericTraits<TInput>::max();
       
    double min = vnl_math_min(r, g);
    min = vnl_math_min(min, b);
    double max = vnl_math_max(r, g);
    max = vnl_math_max(max, b);
    double delta = max - min;
    
    double H = 0;
    double S = 0;
    double V = max;
    
    if(delta == 0)
      {
	  H = 0;
	  S = 0;
	  }
    else
      {
	  S = delta / max;
	  
	  double r2 = (((max - r) / 6.0) + (max / 2.0)) / max;
	  double g2 = (((max - g) / 6.0) + (max / 2.0)) / max;
	  double b2 = (((max - b) / 6.0) + (max / 2.0)) / max;
	  
	  if(r == max)
	    H = b2 - g2;
	  else if(g == max)
	    H = (1.0 / 3.0) + r2 - b2;
	  else if (b == max)
	    H = (2.0 / 3.0) + g2 -r2;
	    
	  if(H < 0)
	    H += 1;
	  if(H > 1)
	    H -= 1;
	  }

    output[0] = static_cast<TInput>(H); // H
    output[1] = static_cast<TInput>(S); // S
    output[2] = static_cast<TInput>(V); // V
    
    return output;
    }

  /** Read access to the RGBToHSVColorSpace component */
  inline ExternalType Get( const InternalType & input ) const
    {
    // Normalize RGB values.
    double r = (double)input[0] / (double)NumericTraits<TInput>::max();
    double g = (double)input[1] / (double)NumericTraits<TInput>::max();
    double b = (double)input[2] / (double)NumericTraits<TInput>::max();
       
    double min = vnl_math_min(r, g);
    min = vnl_math_min(min, b);
    double max = vnl_math_max(r, g);
    max = vnl_math_max(max, b);
    double delta = max - min;
    
    double H = 0;
    double S = 0;
    double V = max;
    
    if(delta == 0)
      {
	  H = 0;
	  S = 0;
	  }
    else
      {
	  S = delta / max;
	  
	  double r2 = (((max - r) / 6.0) + (max / 2.0)) / max;
	  double g2 = (((max - g) / 6.0) + (max / 2.0)) / max;
	  double b2 = (((max - b) / 6.0) + (max / 2.0)) / max;
	  
	  if(r == max)
	    H = b2 - g2;
	  else if(g == max)
	    H = (1.0 / 3.0) + r2 - b2;
	  else if (b == max)
	    H = (2.0 / 3.0) + g2 -r2;
	    
	  if(H < 0)
	    H += 1;
	  if(H > 1)
	    H -= 1;
	  }

    ExternalType output;
    output[0] = static_cast<TOutput>(H); // H
    output[1] = static_cast<TOutput>(S); // S
    output[2] = static_cast<TOutput>(V); // V
    
    return output;
    }

private:
};
  
}  // end namespace Accessor
}  // end namespace itk

#endif
