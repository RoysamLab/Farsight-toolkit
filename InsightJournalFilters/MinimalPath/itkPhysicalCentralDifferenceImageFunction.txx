/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkPhysicalCentralDifferenceImageFunction.txx,v $
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPhysicalCentralDifferenceImageFunction_txx
#define __itkPhysicalCentralDifferenceImageFunction_txx

#include "itkPhysicalCentralDifferenceImageFunction.h"

namespace itk
{


	/**
	* Constructor
	*/
	template <class TInputImage, class TCoordRep>
	PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>
		::PhysicalCentralDifferenceImageFunction()
	{
		m_Interpolator = InterpolateImageFunctionType::New();
	}


	/**
	*
	*/
	template <class TInputImage, class TCoordRep>
	void
		PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>
		::PrintSelf(std::ostream& os, Indent indent) const
	{
		this->Superclass::PrintSelf(os,indent);
	}


	/**
	*
	*/
	template <class TInputImage, class TCoordRep>
	typename PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>::OutputType
		PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>
		::Evaluate( const PointType& point ) const
	{
		OutputType derivative;
		derivative.Fill( 0.0 );

		const	typename TInputImage::RegionType image_region = Superclass::m_Image->GetBufferedRegion();
		ContinuousIndexType index_of_point;

		for(unsigned int dim=0; dim<TInputImage::ImageDimension; dim++)
		{
			double denominator = 0.0;
			// Get the left neighbor
			PointType pointLeft( point );
			pointLeft[dim] += -1 * Superclass::m_Image->GetSpacing()[dim];
			Superclass::m_Image->TransformPhysicalPointToContinuousIndex(pointLeft, index_of_point);
			if (image_region.IsInside(index_of_point))
			{
				denominator += 1.0;
			}
			else
			{
				pointLeft[dim] = point[dim];
			}
			TCoordRep valueLeft = m_Interpolator->Evaluate( pointLeft );


			// Get the right neighbor
			PointType pointRight( point );
			pointRight[dim] += 1 * Superclass::m_Image->GetSpacing()[dim];
			Superclass::m_Image->TransformPhysicalPointToContinuousIndex(pointRight, index_of_point);
			if (image_region.IsInside(index_of_point))
			{ 
				denominator += 1.0;
			}
			else
			{
				pointRight[dim] = point[dim];
			}
			TCoordRep valueRight = m_Interpolator->Evaluate( pointRight );

			if (denominator == 0.0)
			{
				derivative[dim] = 0.0;	//one-pixel wide case
			}
			else
			{	
				// Compute derivative
				derivative[dim] = (valueRight - valueLeft) / 
					(denominator * Superclass::m_Image->GetSpacing()[dim]);
			}			
		}
		return ( derivative );
	}

} // end namespace itk

#endif
