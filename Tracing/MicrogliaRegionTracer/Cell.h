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
#include "AspectRatioResampler.h"

class Cell
{
private:
	typedef itk::CovariantVector< double, 3 >					GradientVectorType;

public:
	typedef	unsigned char										InputPixelType;
	typedef fregl_roi< InputPixelType >::ImageType				ImageType;
    
	typedef itk::Image< double, 3 >								LoGImageType;
	typedef itk::Image< double, 3 >								VesselnessImageType;
	typedef itk::Image< double, 3 >								DistanceImageType;
	typedef itk::Image< unsigned char, 3 >						SomaImageType;
	typedef SomaImageType										MaskImageType;
	typedef itk::Image< unsigned short, 3 >						LabelImageType;
	typedef itk::ImageFileReader < SomaImageType >				SomaReaderType;
	typedef	itk::Image< GradientVectorType, 3 >					GVFImageType;
	typedef GVFImageType										GradientImageType;
    
public:
	explicit Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z, double aspect_ratio);
	
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
	
	void                    CreateIsotropicImage();
    void					CreateLoGImage();
	void					CreateVesselnessImage();
	void					CreateSpeedImage();

	//Writes various images
	static void             WriteImage(const std::string & filename, const itk::Image< unsigned char, 3>::Pointer & image);
	static void             WriteImage(const std::string & filename, const itk::Image< unsigned short, 3>::Pointer & image);
	static void             WriteImage(const std::string & filename, const itk::Image< float , 3 >::Pointer & image);
	static void             WriteImage(const std::string & filename, const itk::Image< double , 3 >::Pointer & image);
    
	void					WriteTreeToSWCFile(Tree* tree, std::string filename, std::string filename_local);
    void                    WriteTreeToSWCFileDepthFirst(Tree* tree, std::ofstream &traceFile, std::ofstream &traceFile_local);
	void					WriteTreeToSWCFileBreadthFirst(Tree* tree, std::ofstream &traceFile, std::ofstream &traceFile_local);
	void                    WriteLinkToParent(Node* node, itk::uint64_t tree_depth, std::ofstream &traceFile, std::ofstream &traceFile_local);

private:    
    void					CreateGVFVesselnessImage(float noise_level, int num_iterations);
	void					CreateGVFImage(float noise_level, int num_iterations);
	double					GetVesselnessValue(const GradientVectorType & grad_Dx_vector, const GradientVectorType & grad_Dy_vector, const GradientVectorType & grad_Dz_vector);
    
	void					CreateHessianVesselnessImage();

	void					SplitITKCovariantVectorImage(const itk::Image< itk::CovariantVector< double, 3 >, 3 >::Pointer & covar_image, itk::Image< double, 3>::Pointer & x_image, itk::Image< double, 3>::Pointer & y_image, itk::Image< double, 3>::Pointer & z_image);

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
	ImageType::Pointer				isotropic_image;
	GVFImageType::Pointer			gvf_image;
	VesselnessImageType::Pointer	adjusted_seed_img;
	ImageType::Pointer				threshold_label_mask;
    
	itk::int64_t next_available_ID;
    
private:
	itk::uint64_t cell_x;
	itk::uint64_t cell_y;
	itk::uint64_t cell_z;
    
	double aspect_ratio;

	ImageType::PointType roi_origin;
	ImageType::SizeType roi_size;
    
	ImageType::SizeType cell_requested_size;
	ImageType::IndexType shift_index;

	AspectRatioResampler::SamplingType sampling_type;
};

#endif
