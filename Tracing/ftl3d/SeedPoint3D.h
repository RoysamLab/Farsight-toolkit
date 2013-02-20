/**
 \brief Class for seed points in a 3D volume. 
 \author $ Author: Amit Mukherjee $
 \version $Revision 1.0 $
 \date May 05 2009
*/
// Copyright Farsight Toolkit , Rensselaer Polytechnic institute Troy NY 12180.

#ifndef SEEDPOINT3D_H
#define SEEDPOINT3D_H
#include "itkImage.h"
#include <vnl/vnl_vector_fixed.h>

typedef vnl_vector_fixed<double,3> Vect3;
typedef itk::Image<float, 3> ImageType3D;
typedef itk::Image<float, 2> ImageType2D;

class SeedPoint3D	{
	public :
	Vect3 getPosition() {return(pos);}
	float getScale() {return(Scale);}
	float getIntensity() {return(Intensity);}
	inline ImageType3D::RegionType::IndexType getIndex()	{
		ImageType3D::RegionType::IndexType ndx;
		ndx[0] = int(pos[0]);
		ndx[1] = int(pos[1]);
		ndx[2] = int(pos[2]);
		return ndx;
	}

	inline void setXYPosition(ImageType2D::RegionType::IndexType ndx) {
		pos[0] = static_cast<double>(ndx[0]);
		pos[1] = static_cast<double>(ndx[1]);
	}

	void setSeed3D(ImageType3D::IndexType ndx, ImageType3D::PixelType val, long sc) {
		pos[0] = static_cast<double>(ndx[0]);
		pos[1] = static_cast<double>(ndx[1]);
		pos[2] = static_cast<double>(ndx[2]);
		Intensity = static_cast<float>(val);
		Scale = static_cast<float>(sc)/3.0;
	}

	inline void setScale(float i) {Scale = i;}
	inline void setIntensity(float i) {Intensity = i;}
	inline void setZPosition(unsigned int i) {pos[2] = static_cast<double>(i);}
	void PrintSelf();

	private:
	Vect3 pos;
    float Intensity;
    float Scale;
 };

 #endif
