#ifndef _Common_h
#define _Common_h

#include <iostream>
#include <fstream>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkMedianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "QuickView.h"
#include "itkImageToVTKImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkDivideImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkStatisticsImageFilter.h"

#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkActor.h"
#include "vtkImageActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkVolume.h"
#include "vtkPiecewiseFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkFixedPointVolumeRayCastMapper.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "vtkGPUVolumeRayCastMapper.h"
#include "vtkOpenGLGPUVolumeRayCastMapper.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkIndent.h"
#include "vtkCamera.h"
#include "vtkImageCast.h"
#include "vtkImageViewer2.h"

typedef float PixelType;
typedef itk::Image<PixelType, 3> ImageType3D;
typedef unsigned char RenderPixelType;
typedef itk::Image<RenderPixelType, 3> RenderImageType3D;

typedef itk::ImageToVTKImageFilter<RenderImageType3D> ITKToVTKConnectorType;
typedef itk::ShiftScaleImageFilter<ImageType3D,ImageType3D> ShiftScaleFilterType;
typedef itk::RegionOfInterestImageFilter<ImageType3D, ImageType3D> VolumeOfInterestFilterType;
typedef itk::MinimumMaximumImageCalculator<ImageType3D> MinMaxCalculatorType;
typedef itk::DivideImageFilter<ImageType3D, ImageType3D, ImageType3D> DivideImageFilterType;
typedef itk::StatisticsImageFilter<ImageType3D> StatisticsFilterType;
typedef itk::ImageDuplicator<ImageType3D> DuplicatorType;

namespace Common{
	
	/**
	 * Read a 3D TIFF image to disk.
	 * (data path with extension, empty data pointer)
	 */
	void ReadImage3D(const std::string&, ImageType3D::Pointer&);

	/**
	 * Wirte a 3D TIFF image to disk. 
	 * (data path with extension, data pointer)
	 */
	void WriteTIFFImage3D(const std::string&, ImageType3D::Pointer&);

	/**
	 * Wirte a 3D image as mhd file to disk. 
	 * (data path with extension, data pointer)
	 */
	void WriteImage3D(const std::string&, ImageType3D::Pointer&);

	/**
	 * Curvature anisotropic diffusion filter used by Amit.
	 * (N_iterations, conductance, data pointer)
	 */
	void CurvatureAnisotropicDiffusion(int&, int&, ImageType3D::Pointer&); 

	/**
	 * Median filtering.
	 * (radius, data_pointet)
	 */
	void MedianFilter(const int, ImageType3D::Pointer&);

	/**
	 * Calculate the gradient vector field and write the result in mhd files.
	 * (sigma for smoothing, path for writing files, data pointer)
	 */
	void GVFDiffusion(float&, const std::string&, ImageType3D::Pointer&);
	
	/** Clip the given data using the epsilon value.
	 * (epsilon, float data)
	 */
	float inline EpsilonClip(double, float);

	/** Clip the given data using the epsilon value.
	 * (epsilon, double data)
	 */
	double inline EpsilonClip(double, double);

	/** Render 3D ITK image using VTK.
	 * (ITK image ptr)
	 */
	void RenderImage3D(RenderImageType3D::Pointer);

	/** Rescale data for rendering. Rescaling since input data is 16 bit.
	 * (ITK image ptr, ITK image ptr for rendering)
	 */
	void RescaleDataForRendering(ImageType3D::Pointer, RenderImageType3D::Pointer&);
	
	/** Normalize data using its max value, return tge maax value.
	 * (input data, data to normalize)
	 */
	PixelType NormalizeData(ImageType3D::Pointer, ImageType3D::Pointer&);

	// Cannot add new functions to namespace!!
	/** Return the sign of the value
	 * (value)
	 */
	//int inline GetSign(double);
} // end namespace Common

#endif // _Common_h
