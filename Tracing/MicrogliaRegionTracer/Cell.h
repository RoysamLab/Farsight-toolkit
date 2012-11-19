#ifndef Cell_H
#define Cell_H

#include <deque>

#include "itkImage.h"
#include "itkIntTypes.h"
#include "itkImageFileReader.h"

#include "fregl/fregl_roi.h"

#include <fstream>
#include <cstring>

#include "ftkTimeStampOverflowSafeUpdate.h"

#include "Cell.h"
#include "Tree.h" 

class Cell
{
public:
	typedef	unsigned char								InputPixelType;
	typedef fregl_roi< InputPixelType >::ImageType		ImageType;
    
	typedef itk::Image< float, 3 >						LoGImageType;
	typedef itk::Image< float, 3 >						VesselnessImageType;
	typedef itk::Image< float, 3 >						DistanceImageType;
	typedef itk::Image< unsigned char, 3 >				SomaImageType;
	typedef SomaImageType								MaskImageType;
	typedef itk::Image< unsigned short, 3 >				LabelImageType;
	typedef itk::ImageFileReader < SomaImageType >		SomaReaderType;
    
public:
	explicit Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z);
	
	itk::uint64_t           getX() const;
	itk::uint64_t           getY() const;
	itk::uint64_t           getZ() const;
    
	void                    SetSize(const ImageType::SizeType & roi_size);
	ImageType::SizeType     GetSize() const;
    
	void                    SetOrigin(const ImageType::PointType & roi_origin);
	ImageType::PointType    GetOrigin() const;
    
	void                    SetRequestedSize(const ImageType::SizeType & cell_requested_size);
	
    ImageType::SizeType     GetRequestedSize() const;
    
	void                    SetShiftIndex(ImageType::IndexType shift_index);
	ImageType::IndexType    GetShiftIndex() const;
    
	void                    ComputeCriticalPointsVector(const ImageType::Pointer & critical_points_image);
    
	//Various methods to perform filters on image
	void                    GetMask(const std::string & soma_filename);
	void                    ComputeMaskedImage();
    
	//Writes various images
	void                    WriteImage(const std::string & filename, const itk::Image< unsigned char, 3>::Pointer & image) const;
	void                    WriteImage(const std::string & filename, const itk::Image< unsigned short, 3>::Pointer & image) const;
	void                    WriteImage(const std::string & filename, const itk::Image< float , 3 >::Pointer & image) const;
    
    void                    WriteTreeToSWCFile(Tree* tree, std::string filename, std::string filename_local);
	void                    WriteLinkToParent(Node* node, itk::uint64_t tree_depth, std::ofstream &traceFile, std::ofstream &traceFile_local);
	
    
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
