#ifndef Cell_H
#define Cell_H

#include <deque>

#include "itkImage.h"
#include "itkIntTypes.h"

#include "fregl/fregl_roi.h"


class Cell
{
public:
	typedef	unsigned char							InputPixelType;
	typedef fregl_roi< InputPixelType >::ImageType	ImageType;
	
	//typedef ImageType::Pointer					ImageTypePointer;
	//typedef ImageType::SizeType					ImageSizeType;
	//typedef ImageType::IndexType					ImageIndexType;

	typedef itk::Image<float, 3>					LoGImageType;
	typedef itk::Image<float, 3>					VesselnessImageType;

public:
	Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z);
	
	itk::uint64_t getX() const;
	itk::uint64_t getY() const;
	itk::uint64_t getZ() const;

	void setSize(ImageType::SizeType size);
	ImageType::SizeType getSize();

	void setRequestedSize(ImageType::SizeType cell_requested_size);
	ImageType::SizeType getRequestedSize();

	void setShiftIndex(ImageType::IndexType shift_index);
	ImageType::IndexType getShiftIndex();

public:
	std::deque< ImageType::IndexType > critical_points_queue; //deque is double-ended queue
	
	ImageType::Pointer image;
	ImageType::Pointer critical_point_image;
	LoGImageType::Pointer multiscale_LoG_image;
	VesselnessImageType::Pointer vesselness_image;
	ImageType::Pointer masked_image;

private:
	itk::uint64_t cell_x;
	itk::uint64_t cell_y;
	itk::uint64_t cell_z;

	ImageType::SizeType cell_size;
	ImageType::SizeType cell_requested_size;
	ImageType::IndexType shift_index;
};

#endif
