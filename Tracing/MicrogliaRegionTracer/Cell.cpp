#include "Cell.h"

#include "itkMaskNegatedImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkPowImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientVectorFlowImageFilter.h"
#include "itkVectorImageToImageAdaptor.h"
#include "AspectRatioResampler.h"

#include "EigenAnalysis.h"

#include "LoG.h"

Cell::Cell(itk::uint64_t cell_x, itk::uint64_t cell_y, itk::uint64_t cell_z, double aspect_ratio) :
	next_available_ID(1),
	cell_x(cell_x),
	cell_y(cell_y),
	cell_z(cell_z),
	aspect_ratio(aspect_ratio),
    sampling_type(AspectRatioResampler::DownSample)
{
}

itk::uint64_t Cell::getX() const
{
	return cell_x;
}

itk::uint64_t Cell::getY() const
{
	return cell_y;
}

itk::uint64_t Cell::getZ() const
{
	return cell_z;
}

void Cell::SetSize(const ImageType::SizeType & roi_size)
{
	this->roi_size = roi_size;
}

Cell::ImageType::SizeType Cell::GetSize() const
{
	return roi_size;
}

void Cell::SetOrigin(const ImageType::PointType & roi_origin)
{
	this->roi_origin = roi_origin;
}

Cell::ImageType::PointType Cell::GetOrigin() const
{
	return roi_origin;
}

void Cell::SetRequestedSize(const ImageType::SizeType & cell_requested_size)
{
	this->cell_requested_size = cell_requested_size;
}

Cell::ImageType::SizeType Cell::GetRequestedSize() const
{
	return cell_requested_size;
}

void Cell::SetShiftIndex(ImageType::IndexType shift_index)
{
	this->shift_index = shift_index;
}

Cell::ImageType::IndexType Cell::GetShiftIndex() const
{
	return shift_index;
}

void Cell::ComputeCriticalPointsVector(const ImageType::Pointer & critical_points_image)
{
	critical_points_queue.clear();

	itk::ImageRegionIterator< ImageType > critical_points_img_iter(critical_points_image, critical_points_image->GetLargestPossibleRegion());
	critical_points_img_iter.GoToBegin();
	while(!critical_points_img_iter.IsAtEnd())
	{
		if (critical_points_img_iter.Get() != 0)
			critical_points_queue.push_back(critical_points_img_iter.GetIndex());

		++critical_points_img_iter;
	}
}

void Cell::GetMask(const std::string & soma_filename)
{
	typedef itk::ImageFileReader< SomaImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(soma_filename);

	//**WARNING** PURPOSELY COMMENTED OUT SO YOU DO NOT TRY TO WRITE THIS CODE, THIS KILLS PERFORMANCE BECAUSE IT ATTEMPTS TO READ THE ENTIRE IMAGE FOR EACH CELL INSTEAD OF JUST THE ROI
	//try
	//{
	//	reader->Update();
	//}
	//catch (itk::ExceptionObject &err)
	//{
	//	std::cerr << "reader Exception: " << err << std::endl;
	//}

	typedef itk::RegionOfInterestImageFilter< SomaImageType, MaskImageType > ROIFilterType;
	ROIFilterType::Pointer roi_filter = ROIFilterType::New();

	ImageType::IndexType start;
	start[0] = this->roi_origin[0];
	start[1] = this->roi_origin[1];
	start[2] = this->roi_origin[2];

	ImageType::SizeType size = this->roi_size;

	ImageType::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);

	roi_filter->SetRegionOfInterest(desiredRegion);
	roi_filter->SetInput(reader->GetOutput());

	try
	{
		//roi_filter->Update();
		ftk::TimeStampOverflowSafeUpdate(roi_filter.GetPointer());
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "roi_filter Exception: " << err << std::endl;
	}

	this->mask = roi_filter->GetOutput();
	this->mask->DisconnectPipeline();	//Disconnect pipeline so we don't propagate...

	ImageType::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	this->mask->SetOrigin(origin);

	//Make the file name of the mask image
	std::stringstream mask_filename_stream;
	mask_filename_stream << cell_x << "_" << cell_y << "_" << cell_z << "_mask.TIF";	//X_Y_Z_masked.TIF

	//Write the masked cell image
	//WriteImage(mask_filename_stream.str(), this->mask);

	//Get the label image from the binary image
	typedef itk::BinaryImageToLabelMapFilter< MaskImageType > BinaryToLabelFilterType;
	BinaryToLabelFilterType::Pointer labelMapFilter = BinaryToLabelFilterType::New();
	labelMapFilter->SetInput(this->mask);

	try
	{
		ftk::TimeStampOverflowSafeUpdate( labelMapFilter.GetPointer() );
		//labelMapFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "labelMapFilter exception: " << err << std::endl;
		std::cerr << "Mask image: " << this->mask << std::endl;
		std::cerr << this->mask << std::endl;
		std::cerr << labelMapFilter << std::endl;
	}

	BinaryToLabelFilterType::OutputImageType::Pointer label_map_image = labelMapFilter->GetOutput();
	label_map_image->DisconnectPipeline();

	typedef itk::LabelMapToLabelImageFilter< BinaryToLabelFilterType::OutputImageType, LabelImageType > LabelMapToLabelImageFilterType;
	LabelMapToLabelImageFilterType::Pointer labelImageFilter = LabelMapToLabelImageFilterType::New();
	labelImageFilter->SetInput(label_map_image);
	try
	{
		ftk::TimeStampOverflowSafeUpdate( labelImageFilter.GetPointer() );
		//labelImageFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "labelImageFilter exception: " << err << std::endl;
	}

	this->soma_label_image = labelImageFilter->GetOutput();
	this->soma_label_image->DisconnectPipeline();	//Disconnect pipeline so we don't propagate...
}

void Cell::ComputeMaskedImage()
{
	typedef itk::MaskNegatedImageFilter< ImageType, MaskImageType, ImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	maskFilter->SetMaskImage(this->mask);
	maskFilter->SetInput(this->image);
	try
	{
		//ftk::TimeStampOverflowSafeUpdate( maskFilter.GetPointer() );
		maskFilter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "maskFilter Exception: " << err << std::endl;
	}

	this->masked_image = maskFilter->GetOutput();
	this->masked_image->DisconnectPipeline();		//Disconnect pipeline so we don't propagate...

	//Make the file name of the masked cell image
	std::stringstream masked_cell_filename_stream;
	masked_cell_filename_stream << cell_x << "_" << cell_y << "_" << cell_z << "_masked.TIF";	//X_Y_Z_masked.TIF

	//Write the masked cell image
	//WriteImage(masked_cell_filename_stream.str(), this->masked_image);
}

void Cell::WriteImage(const std::string & filename, const itk::Image< unsigned char, 3>::Pointer & image)
{
	typedef itk::ImageFileWriter< itk::Image< unsigned char, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		writer->Update();
		//ftk::TimeStampOverflowSafeUpdate( writer.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
		std::cerr << writer << std::endl;
	}
}

void Cell::WriteImage(const std::string & filename, const itk::Image< unsigned short, 3>::Pointer & image)
{
	typedef itk::ImageFileWriter< itk::Image< unsigned short, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		writer->Update();
		//ftk::TimeStampOverflowSafeUpdate( writer.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
		std::cerr << writer << std::endl;
	}
}

void Cell::WriteImage(const std::string & filename, const itk::Image< float , 3 >::Pointer & image)
{
	typedef itk::ImageFileWriter< itk::Image< float, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		writer->Update();
		//ftk::TimeStampOverflowSafeUpdate( writer.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
		std::cerr << writer << std::endl;
	}
}

void Cell::WriteImage(const std::string & filename, const itk::Image< itk::CovariantVector < float, 3 >, 3 >::Pointer & image)
{
	typedef itk::ImageFileWriter< itk::Image< itk::CovariantVector < float, 3 >, 3 > > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(image);	//image is from function parameters!
	writer->SetFileName(filename);
	try
	{
		writer->Update();
		//ftk::TimeStampOverflowSafeUpdate( writer.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "writer Exception: " << err << std::endl;
		std::cerr << writer << std::endl;
	}
}

/*	This function will take a tree structure and write out the "local" and "global" SWC file corresponding to that tree structure.
*	Local here means an SWC file that will fit into the region of interest we are interested in and global means that it will fit into the montage */
void Cell::WriteTreeToSWCFile(Tree * tree, std::string filename, std::string filename_local)
{
	std::cout << "Entering WriteTreeToSWCFile" << std::endl;
	std::ofstream traceFile, traceFile_local;

	std::cout << "Opening " << filename << std::endl;
	traceFile.open(filename.c_str());
	std::cout << "Opening " << filename_local << std::endl;
	traceFile_local.open(filename_local.c_str());

	Node* root = tree->GetRoot();

	itk::uint64_t tree_depth = 0; //root node is defined as tree depth 0
	WriteLinkToParent(root, tree_depth, traceFile, traceFile_local);	//Recursive function that does the actual tree traversal and writing the SWC lines

	traceFile.close();
	traceFile_local.close();
}

//This function does a depth-first traversal of the tree, maybe it makes sense to do a breadth-first traversal for some uses?
//CAREFUL: This function may overflow the stack due to recursion pushing parameters onto the stack and will crash if you have a tree deep enough. Rewrite iteratively or increase stack size if this is the case...
void Cell::WriteLinkToParent(Node* node, itk::uint64_t tree_depth, std::ofstream &traceFile, std::ofstream &traceFile_local)
{
	//Calculate some node indices
	ImageType::PointType node_index, node_index_local;
	node_index_local[0] = node->x;
	node_index_local[1] = node->y;
	node_index_local[2] = node->z;
	node_index[0] = node_index_local[0] + this->getX() - this->GetRequestedSize()[0]/2 - this->GetShiftIndex()[0];
	node_index[1] = node_index_local[1] + this->getY() - this->GetRequestedSize()[1]/2 - this->GetShiftIndex()[1];
	node_index[2] = node_index_local[2] + this->getZ() - this->GetRequestedSize()[2]/2 - this->GetShiftIndex()[2];

	itk::int64_t parent_node_id;
	if (node->GetParent() == NULL)
		parent_node_id = -1; //This node has no parent, so this is the root node, so parent ID is -1
	else
		parent_node_id = node->GetParent()->getID();

	std::vector< Node* > children = node->GetChildren();

	if (tree_depth == 1 && children.size() == 0)
		return;														//BASE CASE: Don't write out a trace if we are at depth one and have no children because we are a trace to the edge of the soma

	//Write out the SWC lines
	traceFile		<< node->getID() << " 3 " << node_index[0]			<< " " << node_index[1]			<< " "  << node_index[2]		<< " 1 " << parent_node_id << std::endl;
	traceFile_local << node->getID() << " 3 " << node_index_local[0]	<< " " << node_index_local[1]	<< " "  << node_index_local[2]	<< " 1 " << parent_node_id << std::endl;

	if (children.size() == 0)
		return;														//BASE CASE: No children, so don't visit them

	std::vector< Node* >::iterator children_iter;
	for (children_iter = children.begin(); children_iter != children.end(); ++children_iter)
	{
		Node* child_node = *children_iter;
		WriteLinkToParent(child_node, tree_depth+1, traceFile, traceFile_local);
	}

	return;															//BASE CASE: Finished visiting the entire subtree that is rooted at this node
}

/* This function resamples the image and creates an isotropic image according to the ratio set in the constructor */
void Cell::CreateIsotropicImage()
{
    std::cout << "Creating Isotropic Image" << std::endl;
	this->isotropic_image = AspectRatioResampler::ResampleImage< ImageType >(this->image, this->aspect_ratio, this->sampling_type);
    
    //Make the file name of the raw cell image
	std::stringstream isotropic_image_filename_stream;
	isotropic_image_filename_stream << this->getX() << "_" << this->getY() << "_" << this->getZ() << "_isotropic.nrrd";
	WriteImage(isotropic_image_filename_stream.str(), this->isotropic_image);
}

Cell::GVFImageType::Pointer Cell::CreateGVFImage(float noise_level, int num_iterations)
{
	//Gradient filter
	typedef itk::GradientImageFilter< ImageType > GradientImageFilterType;
	GradientImageFilterType::Pointer gradient_filter = GradientImageFilterType::New();

	gradient_filter->SetInput(this->isotropic_image);

	try
	{
		gradient_filter->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << "Exception caught: " << err << std::endl;
	}

	//Gradient Vector Flow filter
	typedef itk::GradientVectorFlowImageFilter< GradientImageFilterType::OutputImageType, GradientImageFilterType::OutputImageType > GradientVectorFlowFilterType;
	GradientVectorFlowFilterType::Pointer gvf_filter = GradientVectorFlowFilterType::New();

	gvf_filter->SetInput(gradient_filter->GetOutput());
	gvf_filter->SetNoiseLevel(noise_level);
	gvf_filter->SetIterationNum(num_iterations);

	try
	{
		gvf_filter->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << "Exception caught: " << err << std::endl;
	}

	GVFImageType::Pointer gvf_image = gvf_filter->GetOutput();

	//Make the file name of the GVF image
	std::stringstream gvf_image_filename_stream;
	gvf_image_filename_stream << this->getX() << "_" << this->getY() << "_" << this->getZ() << "_gvf.mhd";
	//WriteImage(gvf_image_filename_stream.str(), gvf_image);

	return gvf_image;
}

void Cell::CreateGVFVesselnessImage(float noise_level, int num_iterations)
{
	std::cout << "Creating Vesselness image by Gradient Vector Flow" << std::endl;
	GVFImageType::Pointer gvf_image = CreateGVFImage(noise_level, num_iterations);

	//Separate the GVF image into Partial Derivatives
	typedef itk::Image< float, 3 > PartialDerivativeImageType;

	PartialDerivativeImageType::Pointer Dx = PartialDerivativeImageType::New();
	PartialDerivativeImageType::Pointer Dy = PartialDerivativeImageType::New();
	PartialDerivativeImageType::Pointer Dz = PartialDerivativeImageType::New();

	PartialDerivativeImageType::SizeType size = gvf_image->GetLargestPossibleRegion().GetSize();
	PartialDerivativeImageType::IndexType start;
	start.Fill(0);
	PartialDerivativeImageType::RegionType region(start, size);

	Dx->SetRegions(region);
	Dy->SetRegions(region);
	Dz->SetRegions(region);

	Dx->Allocate();
	Dy->Allocate();
	Dz->Allocate();

	itk::ImageRegionConstIterator< GVFImageType > gvf_image_iter(gvf_image, gvf_image->GetLargestPossibleRegion());
	itk::ImageRegionIterator< PartialDerivativeImageType > Dx_iter(Dx, Dx->GetLargestPossibleRegion());
	itk::ImageRegionIterator< PartialDerivativeImageType > Dy_iter(Dy, Dy->GetLargestPossibleRegion());
	itk::ImageRegionIterator< PartialDerivativeImageType > Dz_iter(Dz, Dz->GetLargestPossibleRegion());

	while (!gvf_image_iter.IsAtEnd())
	{
		GVFImageType::PixelType vector = gvf_image_iter.Get();

		Dx_iter.Set(vector[0]);
		Dy_iter.Set(vector[1]);
		Dz_iter.Set(vector[2]);

		++gvf_image_iter;
		++Dx_iter;
		++Dy_iter;
		++Dz_iter;
	}
	
	std::ostringstream Dx_filename_stream;
	Dx_filename_stream << this->getX() << "_" << this->getY() << "_" << this->getZ() << "_Dx.nrrd";
	//WriteImage(Dx_filename_stream.str(), Dx);

	std::ostringstream Dy_filename_stream;
	Dy_filename_stream << this->getX() << "_" << this->getY() << "_" << this->getZ() << "_Dy.nrrd";
	//WriteImage(Dy_filename_stream.str(), Dy);

	std::ostringstream Dz_filename_stream;
	Dz_filename_stream << this->getX() << "_" << this->getY() << "_" << this->getZ() << "_Dz.nrrd";
	//WriteImage(Dz_filename_stream.str(), Dz);

	//Calculate gradient of the GVF image
	//Gradient of Dx
	typedef itk::GradientImageFilter< PartialDerivativeImageType > GradientImageFilterType;
	GradientImageFilterType::Pointer gradient_filter = GradientImageFilterType::New();
	gradient_filter->SetInput(Dx);
	try
	{
		gradient_filter->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << "Exception caught: " << err << std::endl;
	}
	GradientImageType::Pointer grad_Dx = gradient_filter->GetOutput();
	grad_Dx->DisconnectPipeline();

	//Gradient of Dy
	gradient_filter->SetInput(Dy);
	try
	{
		gradient_filter->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << "Exception caught: " << err << std::endl;
	}
	GradientImageType::Pointer grad_Dy = gradient_filter->GetOutput();
	grad_Dy->DisconnectPipeline();

	//Gradient of Dz
	gradient_filter->SetInput(Dz);
	try
	{
		gradient_filter->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		std::cerr << "Exception caught: " << err << std::endl;
	}
	GradientImageType::Pointer grad_Dz = gradient_filter->GetOutput();
	grad_Dz->DisconnectPipeline();

	//Calculate the vesselness value
	itk::ImageRegionConstIterator< GradientImageType > grad_Dx_iter(grad_Dx, grad_Dx->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator< GradientImageType > grad_Dy_iter(grad_Dy, grad_Dy->GetLargestPossibleRegion());
	itk::ImageRegionConstIterator< GradientImageType > grad_Dz_iter(grad_Dz, grad_Dz->GetLargestPossibleRegion());

	//GradientImageType::RegionType region = grad_Dx->GetLargestPossibleRegion();

	VesselnessImageType::Pointer isotropic_vesselness_image = VesselnessImageType::New();
	isotropic_vesselness_image->SetRegions(region);
	isotropic_vesselness_image->Allocate();

	itk::ImageRegionIterator< VesselnessImageType > vesselness_img_iter(isotropic_vesselness_image, isotropic_vesselness_image->GetLargestPossibleRegion());

	while (!vesselness_img_iter.IsAtEnd())
	{
		GradientVectorType grad_Dx_vector = grad_Dx_iter.Get();
		GradientVectorType grad_Dy_vector = grad_Dy_iter.Get();
		GradientVectorType grad_Dz_vector = grad_Dz_iter.Get();

		vesselness_img_iter.Set(GetVesselnessValue(grad_Dx_vector, grad_Dy_vector, grad_Dz_vector));

		++vesselness_img_iter;
		++grad_Dx_iter;
		++grad_Dy_iter;
		++grad_Dz_iter;
	}

	std::ostringstream isotropic_vesselness_filename_stream;
	isotropic_vesselness_filename_stream << this->getX() << "_" << this->getY() << "_" << this->getZ() << "_isotropic_vesselness.nrrd";
	//WriteImage(isotropic_vesselness_filename_stream.str(), isotropic_vesselness_image);

	this->vesselness_image = AspectRatioResampler::UnsampleImage< VesselnessImageType >(isotropic_vesselness_image, this->aspect_ratio, this->sampling_type);
}

float Cell::GetVesselnessValue(const GradientVectorType & grad_Dx_vector, const GradientVectorType & grad_Dy_vector, const GradientVectorType & grad_Dz_vector)
{
	float Dxx = grad_Dx_vector[0];
	float Dxy = grad_Dx_vector[1];
	float Dxz = grad_Dx_vector[2];
	float Dyx = grad_Dy_vector[0];
	float Dyy = grad_Dy_vector[1];
	float Dyz = grad_Dy_vector[2];
	float Dzx = grad_Dz_vector[0];
	float Dzy = grad_Dz_vector[1];
	float Dzz = grad_Dz_vector[2];

	double grad_GVF_matrix[3][3];
	grad_GVF_matrix[0][0] = Dxx;
	grad_GVF_matrix[0][1] = Dxy;
	grad_GVF_matrix[0][2] = Dxz;
	grad_GVF_matrix[1][0] = Dyx;
	grad_GVF_matrix[1][1] = Dyy;
	grad_GVF_matrix[1][2] = Dyz;
	grad_GVF_matrix[2][0] = Dzx;
	grad_GVF_matrix[2][1] = Dzy;
	grad_GVF_matrix[2][2] = Dzz;

	double eigenvalues[3];
	double eigenvectors[3][3];

	EigenAnalysis::eigen_decomposition(grad_GVF_matrix, eigenvectors, eigenvalues);

	double Lambda1 = eigenvalues[0];
	double Lambda2 = eigenvalues[1];
	double Lambda3 = eigenvalues[2];

	if ( Lambda2 >= 0.0 || Lambda3 > 0.0 )
	{
		return 0.0;
	}
	else
	{
		static const double FrangiAlpha = 0.5;
		static const double FrangiBeta = 0.5;
		static const double FrangiC = 10;

		const double A = 2 * pow(FrangiAlpha,2);
		const double B = 2 * pow(FrangiBeta,2);
		const double C = 2 * pow(FrangiC,2);

		double Ra  = Lambda2 / Lambda3; 
		double Rb  = Lambda1 / vcl_sqrt ( vnl_math_abs( Lambda2 * Lambda3 )); 
		double S  = vcl_sqrt( pow(Lambda1,2) + pow(Lambda2,2) + pow(Lambda3,2));

		double vesMeasure_1  = ( 1 - vcl_exp(-1.0*(( vnl_math_sqr( Ra ) ) / ( A ))) );
		double vesMeasure_2  = vcl_exp ( -1.0 * ((vnl_math_sqr( Rb )) /  ( B )));
		double vesMeasure_3  = ( 1 - vcl_exp( -1.0 * (( vnl_math_sqr( S )) / ( C ))) );

		float V_Saliency = vesMeasure_1 * vesMeasure_2 * vesMeasure_3;

		return V_Saliency;
	}
}

void Cell::CreateLoGImage()
{
	//Calculate the LoG on multiple scales and store into an image
	std::cout << "Calculating Multiscale LoG" << std::endl;
	LoGImageType::Pointer resampled_multiscale_LoG_image = LoG::RunMultiScaleLoG(*this);

	//Make the file name of the isotropic LoG image and write it out
	std::stringstream multiscaled_LoG_image_filename_stream;
	multiscaled_LoG_image_filename_stream << this->getX() << "_" << this->getY() << "_" << this->getZ() << "_LoG.nrrd";
	WriteImage(multiscaled_LoG_image_filename_stream.str(), resampled_multiscale_LoG_image);

	ImageType::SizeType outputSize = this->image->GetLargestPossibleRegion().GetSize();

	ImageType::SpacingType outputSpacing;
	outputSpacing[0] = 1.0;
	outputSpacing[1] = 1.0;
	outputSpacing[2] = 1/this->aspect_ratio;
	//outputSpacing[2] = 1.0;

	typedef itk::IdentityTransform< double, 3 > TransformType;
	typedef itk::ResampleImageFilter< LoGImageType, LoGImageType > ResampleImageFilterType;
	ResampleImageFilterType::Pointer resample_filter = ResampleImageFilterType::New();
	resample_filter->SetInput(resampled_multiscale_LoG_image);
	resample_filter->SetSize(outputSize);
	resample_filter->SetOutputSpacing(outputSpacing);
	resample_filter->SetTransform(TransformType::New());

	try
	{
		resample_filter->UpdateLargestPossibleRegion();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "CreateLoGImage() resample_filter exception: " << err << std::endl;
		return;
	}

	this->multiscale_LoG_image = resample_filter->GetOutput();

	ImageType::SpacingType spacing;
	spacing.Fill(1.0);
	this->multiscale_LoG_image->SetSpacing(spacing);	//Convert back to the original spacing
}

void Cell::CreateVesselnessImage()
{
	CreateGVFVesselnessImage(2000.0, 10);
	//CreateHessianVesselnessImage();

	std::ostringstream vesselness_filename_stream;
	vesselness_filename_stream << this->getX() << "_" << this->getY() << "_" << this->getZ() << "_vesselness.nrrd";
	WriteImage(vesselness_filename_stream.str(), this->vesselness_image);
}

void Cell::CreateHessianVesselnessImage()
{
    std::cout << "Creating Vesselness Image by Multiscale Hessian method" << std::endl;
    
    typedef itk::Hessian3DToVesselnessMeasureImageFilter< float > VesselnessFilterType;
	VesselnessFilterType::Pointer vesselness_filter = VesselnessFilterType::New();
	vesselness_filter->SetAlpha1(0.5);
	vesselness_filter->SetAlpha2(0.5);

	typedef itk::SymmetricSecondRankTensor< double, 3 >	HessianTensorType;
	typedef itk::Image< HessianTensorType, 3 >			HessianImageType;

	typedef itk::MultiScaleHessianBasedMeasureImageFilter< ImageType, HessianImageType, VesselnessImageType > MultiscaleHessianFilterType;

	MultiscaleHessianFilterType::Pointer multiscale_hessian_filter = MultiscaleHessianFilterType::New();
	//multiscale_hessian_filter->SetInput(this->image);
	multiscale_hessian_filter->SetInput(this->isotropic_image);
	multiscale_hessian_filter->SetSigmaStepMethodToEquispaced();
	multiscale_hessian_filter->SetSigmaMinimum(1.0);
	multiscale_hessian_filter->SetSigmaMaximum(10.0);
	multiscale_hessian_filter->SetNumberOfSigmaSteps(10);
	multiscale_hessian_filter->SetNonNegativeHessianBasedMeasure(true);
	multiscale_hessian_filter->SetHessianToMeasureFilter(vesselness_filter);

	try
	{
		ftk::TimeStampOverflowSafeUpdate( multiscale_hessian_filter.GetPointer() );
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "multiscale_hessian_filter exception: " << err << std::endl;
	}

	this->vesselness_image = AspectRatioResampler::UnsampleImage< VesselnessImageType >(multiscale_hessian_filter->GetOutput(), this->aspect_ratio, this->sampling_type);
}

void Cell::CreateSpeedImage()
{
	//Normalize vesselness image
	typedef itk::RescaleIntensityImageFilter< DistanceImageType > RescaleIntensityFilterType;
	RescaleIntensityFilterType::Pointer rescale_filter = RescaleIntensityFilterType::New();
	rescale_filter->SetOutputMinimum(0.0);
	rescale_filter->SetOutputMaximum(1.0);
	rescale_filter->SetInput(this->vesselness_image);

	try
	{
		rescale_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "rescale_filter exception: " << err << std::endl;
		return;
	}

	//this->speed_image = rescale_filter->GetOutput();

	//Raise normalized vesselness_image to third power
	typedef itk::PowImageFilter< DistanceImageType > PowImageFilterType;
	PowImageFilterType::Pointer pow_image_filter = PowImageFilterType::New();
	pow_image_filter->SetInput1(rescale_filter->GetOutput());
	pow_image_filter->SetConstant2( 0.33 );

	try
	{
		pow_image_filter->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cerr << "pow_image_filter exception: " << err << std::endl;
		return;
	}

	this->speed_image = pow_image_filter->GetOutput();

	std::stringstream speed_image_filename_stream;
	speed_image_filename_stream << this->getX() << "_" << this->getY() << "_" << this->getZ() << "_speed_image.nrrd";
	WriteImage(speed_image_filename_stream.str(), this->speed_image);
}