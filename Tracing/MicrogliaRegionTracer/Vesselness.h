#ifndef Vesselness_H
#define Vesselness_H


#include "Cell.h"
#include "Hessian.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkMinimumMaximumImageCalculator.h"

class Vesselness
{
private:
	typedef Cell::ImageType				ImageType;
	typedef Cell::VesselnessImageType	VesselnessImageType;
	typedef Cell::HessianImageType		HessianImageType;
	typedef	Cell::HessianTensorType		HessianTensorType;

public:
	Vesselness();
	~Vesselness();

	VesselnessImageType::Pointer	RunVesselness(ImageType::Pointer image, double scale);
	double							ComputeVesselness( HessianTensorType::EigenValuesArrayType::ValueType ev1, HessianTensorType::EigenValuesArrayType::ValueType ev2, HessianTensorType::EigenValuesArrayType::ValueType ev3, double maximum_intensity );
	void							RunMultiscaleVesselness(ImageType::Pointer image);

	void Vesselness::WriteVesselnessImage(std::string filename, VesselnessImageType::Pointer image);
};

#endif
