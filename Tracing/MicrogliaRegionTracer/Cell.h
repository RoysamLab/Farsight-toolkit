#ifndef Cell_H
#define Cell_H

#include <deque>

#include "itkImage.h"
#include "itkIntTypes.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"

#include "fregl/fregl_roi.h"

#include "itkMaskNegatedImageFilter.h"
#include "itkMaximumEntropyThresholdImageFilter.h"
#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryThinningImageFilter3D.h"
#include "itkShiftScaleImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"

#include <fstream>
#include <cstring>

#include "ftkTimeStampOverflowSafeUpdate.h"


class Cell
{
public:
	typedef	unsigned char								InputPixelType;
	typedef fregl_roi< InputPixelType >::ImageType		ImageType;

	typedef itk::Image< float, 3 >						LoGImageType;
	typedef itk::Image< float, 3 >						VesselnessImageType;
	typedef itk::Image< float, 3 >						DistanceImageType;
	typedef itk::Image< unsigned char, 3 >				MaskImageType;
	typedef itk::Image< unsigned char, 3 >				SomaImageType;
	typedef itk::Image< unsigned short, 3 >				LabelImageType;
	typedef itk::ImageFileReader < SomaImageType >		SomaReaderType;

public:
	explicit Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z);
	
	itk::uint64_t getX() const;  
	itk::uint64_t getY() const;
	itk::uint64_t getZ() const;

	void SetSize(ImageType::SizeType roi_size);
	ImageType::SizeType GetSize();

	void SetOrigin(ImageType::PointType roi_origin);
	ImageType::PointType GetOrigin();

	void setRequestedSize(ImageType::SizeType cell_requested_size);
	ImageType::SizeType getRequestedSize();

	void setShiftIndex(ImageType::IndexType shift_index);
	ImageType::IndexType getShiftIndex();

	void ComputeCriticalPointsVector(ImageType::Pointer critical_points_image);

	//Various methods to perform filters on image
	void GetMask(std::string soma_filename);
	void ComputeMaskedImage();

	//Writes various images
	void WriteImage(std::string filename, itk::Image< unsigned char, 3>::Pointer image);
	void WriteImage(std::string filename, itk::Image< unsigned short, 3>::Pointer image);
	void WriteImage(std::string filename, itk::Image< float , 3 >::Pointer image);
	

public:
	std::deque< ImageType::IndexType > critical_points_queue; //deque is double-ended queue
	
	ImageType::Pointer				image;
	ImageType::Pointer				critical_point_image;
	LoGImageType::Pointer			multiscale_LoG_image;
	LoGImageType::Pointer			ridge_image;
	VesselnessImageType::Pointer	vesselness_image;
	ImageType::Pointer				masked_image;
	MaskImageType::Pointer			mask;
	LabelImageType::Pointer			soma_label_image;
	DistanceImageType::Pointer		speed_image;
	ImageType::Pointer				isometric_image;

	itk::int64_t next_available_ID;

private:
	itk::uint64_t cell_x;
	itk::uint64_t cell_y;
	itk::uint64_t cell_z;

	ImageType::PointType roi_origin;
	ImageType::SizeType roi_size;

	ImageType::SizeType cell_requested_size;
	ImageType::IndexType shift_index;

	
};

#endif
