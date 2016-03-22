#ifndef SYNAPSE_SEG_txx
#define SYNAPSE_SEG_txx

#include "Synapse_Seg.h"

// Constructor
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::Synapse_Seg()
{
	nucImage = NULL;
	assImage = NULL;
	synImage = NULL;
	labelImage = NULL;

	this->centroids.clear();
	this->areas.clear();
	this->associate_intensity.clear();
	
	// Defaults Parameters:
	this->paras.nuc_mask_threshold = 0.0;
	this->paras.association_intensity_threshold = 50.0;
	this->paras.association_distance_threshold = 100;
	this->paras.mask_size = 70;
	this->paras.mask_sigma = 0;
	this->paras.syn_threshold = 90.0;
}

// Set Inpuut Images and Label Image
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
bool Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::SetInput(IntensityImagePointer nucImgIn, IntensityImagePointer assImgIn, IntensityImagePointer synImgIn,LabelImagePointer lblImgIn)
{
	if( !(nucImgIn && assImgIn && synImgIn && lblImgIn))
	{std::cout<<"no images to process !"<<std::endl; return false;}

	typename IntensityImageType::RegionType intRegion = nucImgIn->GetLargestPossibleRegion();
	typename LabelImageType::RegionType lblRegion = lblImgIn->GetLargestPossibleRegion();
	if( lblRegion != intRegion )
	{std::cout<<"mismatch size of inensity image and label image process !"<<std::endl; return false;}

	this->nucImage = nucImgIn;
	this->assImage = assImgIn;
	this->synImage = synImgIn;
	this->labelImage = lblImgIn;

	this->num_row = intRegion.GetSize()[0];
	this->num_colum = intRegion.GetSize()[1];

	return true;
}

// Set Parameters from xml file;
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::SetParameters(std::vector<Param > parameters)
{
	this->paras.nuc_mask_threshold = parameters[0].value;
	this->paras.association_intensity_threshold = parameters[1].value;
	this->paras.association_distance_threshold = parameters[2].value;
	this->paras.mask_size = parameters[3].value;
	this->paras.mask_sigma = parameters[4].value;
	this->paras.syn_threshold = parameters[5].value;
	
	return;
}

// Update
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::Update()
{
	std::cout<<"Segmenting!!"<<std::endl;
	this->ComputeLabelStatistics();
	this ->CreatNucMask();
	std::cout<<"Associating!!"<<std::endl;
	this->ComputeMaskedAssociateImage();
	this->ComputeAssociateIntensity();
	std::cout<<"Classifing!!"<<std::endl;
	this->ClassifyCells();
	this->ComputeDistanceTarNonTar();
	this->GetTarNeighbour();
	std::cout<<"Creating Mask!!"<<std::endl;
	this->CreatGlobalMask();
	std::cout<<"Thresholding!!"<<std::endl;
	this->ThresholdSynImage();
	std::cout<<"Masking!!"<<std::endl;
	this->Masking();
	std::cout<<"Finished!!"<<std::endl;
}

// Compute label statistic including area and centroids
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::ComputeLabelStatistics()
{
	typedef itk::LabelGeometryImageFilter< LabelImageType> LabelGeometryImageFilterType;
	LabelGeometryImageFilterType::Pointer labGeometryFilter = LabelGeometryImageFilterType::New();
	labGeometryFilter->SetInput( this->labelImage);
	labGeometryFilter->Update();
	LabelGeometryImageFilterType::LabelsType allLabels = labGeometryFilter->GetLabels();
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt = allLabels.begin();
	for( ++allLabelsIt; allLabelsIt != allLabels.end(); allLabelsIt++ )
	{
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		std::vector< int> centroid;
		centroid.push_back(labGeometryFilter->GetCentroid(labelValue)[0]);
		centroid.push_back(labGeometryFilter->GetCentroid(labelValue)[1]);
		this->centroids.push_back(centroid);
		this->areas.push_back(labGeometryFilter->GetVolume(labelValue));
		this->num_label = labelValue;
	}
}

// Creat nucImage mask to masking associate image
template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::CreatNucMask()
{
	typedef itk::BinaryThresholdImageFilter< LabelImageType, IntensityImageType> BinaryThresholdImageFilterType;
	BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
	thresholdFilter->SetInput(this->labelImage);
	//thresholdFilter->SetLowerThreshold(0);
	thresholdFilter->SetUpperThreshold(0);
	thresholdFilter->SetInsideValue(0);
	thresholdFilter->SetOutsideValue(1);
	this->nuc_mask = thresholdFilter->GetOutput();
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::ComputeMaskedAssociateImage()
{
	typedef itk::MultiplyImageFilter <IntensityImageType, IntensityImageType >MultiplyImageFilterType;
	MultiplyImageFilterType::Pointer multiplyFilter = MultiplyImageFilterType::New ();
	multiplyFilter->SetInput1(this->assImage);
	multiplyFilter->SetInput2(this->nuc_mask);
	this->assmaskedImage = multiplyFilter->GetOutput();
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::ComputeAssociateIntensity()
{
	labIteratorType labIt( labelImage,labelImage->GetLargestPossibleRegion());
	//intIteratorType masked_assIt(assmaskedImage, assmaskedImage->GetRequestedRegion()); 
	intIteratorType masked_assIt(assImage, assImage->GetLargestPossibleRegion()); 

	this->associate_intensity.resize(num_label);
	std::fill (associate_intensity.begin(),associate_intensity.end(), 0);
	for (long int i = 0; i < num_colum * num_row; i++)
	{
		TLPixel label = labIt.Get();
		if (label > 0)
		{ 
			this->associate_intensity[label - 1] += masked_assIt.Get(); 
		}
		labIt++;
		masked_assIt++;
	}

	for(int i = 0; i < num_label; i++)
	{
		associate_intensity[i] = associate_intensity[i] / areas[i];
	}
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::ClassifyCells()
{
	this->labels.resize(num_label);
	this->num_tar = 0;
	this->num_nontar = 0;

	for(int i = 0; i < num_label; i++)
	{
		if(associate_intensity[i] > paras.association_intensity_threshold)
		{
			labels[i] = 1;
			tar_label.push_back(i);
			num_tar++;
		}
		else
		{
			labels[i] = 0;
			non_tar_label.push_back(i);
			num_nontar++;
		}
	}
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::ComputeDistanceTarNonTar()
{
	this->distance_mat.resize(num_tar);
	for(int i = 0; i < num_tar; i++)
	{
		int index1 = tar_label[i];
		distance_mat[i].resize(num_nontar);
		for(int j = 0; j < num_nontar; j++)
		{
			int index2 = non_tar_label[j];
			double temp = (centroids[index1][0] - centroids[index2][0]) * (centroids[index1][0] - centroids[index2][0]) 
									+ (centroids[index1][1] - centroids[index2][1]) * (centroids[index1][1] - centroids[index2][1]);
			distance_mat[i][j] = sqrt(temp);
		}
	}
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::GetTarNeighbour()
{
	this->syn_pairs.clear();
	for(int i = 0; i < num_tar; i++)
	{	
		int min_distance = *std::min_element(distance_mat[i].begin(),distance_mat[i].end());
		if (min_distance < paras.association_distance_threshold)
		{
			int min_index = std::min_element(distance_mat[i].begin(),distance_mat[i].end()) - distance_mat[i].begin();
			syn_pairs.push_back(std::make_pair(tar_label[i],non_tar_label[min_index]));
		}
	}

	num_region = this->syn_pairs.size();
	this->centers.resize(num_region);
	for (int k = 0; k< num_region; k++)
	{
		centers[k].resize(VImageDimension);
		int index1 = syn_pairs[k].first;
		int index2 = syn_pairs[k].second;
		centers[k][0] = (centroids[index1][0] + centroids[index2][0]) / 2;
		centers[k][1] = (centroids[index1][1] + centroids[index2][1]) / 2;
	}
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::CreatGlobalMask()
{
	IntensityImageType::RegionType region = nucImage->GetLargestPossibleRegion();
	mask = IntensityImageType::New();
	mask->SetRegions(region);
	mask->Allocate();
	mask->Update();
 
	IntensityImageType::SizeType regionSize = region.GetSize();
	intIteratorType maskIterator(mask,region);
 
	while(!maskIterator.IsAtEnd())
	{
		int flag = 0;
		for (int  k = 0; k < num_region; k ++)
		{
			double temp_dis = (maskIterator.GetIndex()[0] - centers[k][0]) * (maskIterator.GetIndex()[0] - centers[k][0])
				                     +  (maskIterator.GetIndex()[1] - centers[k][1]) * (maskIterator.GetIndex()[1] - centers[k][1]);
			if( sqrt(temp_dis) < paras.mask_size/2)
			{
				maskIterator.Set(1);
				flag = 1;
				break;
			}
		}
		if(flag == 0)
			maskIterator.Set(0);
		++maskIterator;
	}

	//IntensityImageWriter::Pointer image_writer = IntensityImageWriter::New();
	//image_writer->SetFileName("thisismymask.tif");
	//image_writer->SetInput(mask);
	//image_writer->Update();
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::ThresholdSynImage()
{
	typedef itk::ThresholdImageFilter <IntensityImageType> ThresholdImageFilterType;
 
	ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
	thresholdFilter->SetInput(this->synImage);
	thresholdFilter->ThresholdBelow(paras.syn_threshold);
	thresholdFilter->SetOutsideValue(0);
	thresholdFilter->Update();
	this->masked_syn = thresholdFilter->GetOutput();
	

	//IntensityDuplicatorType::Pointer duplicator = IntensityDuplicatorType::New();

	//IntensityImagePointer temp =  thresholdFilter->GetOutput() ;
	//IntensityImageWriter::Pointer image_writer = IntensityImageWriter::New();
	//image_writer->SetFileName("pleasepleasagain.tif");
	//image_writer->SetInput(masked_syn);
	//image_writer->Update();
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::Masking()
{
	typedef itk::MaskImageFilter< IntensityImageType, IntensityImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	maskFilter->SetInput(masked_syn);
	maskFilter->SetMaskImage(mask);	
	maskFilter->Update();
	this->SynOutImage = maskFilter->GetOutput();
	

	//IntensityImageWriter::Pointer image_writer = IntensityImageWriter::New();
	//image_writer->SetFileName("whyisthishappenning.tif");
	//image_writer->SetInput(SynOutImage);
	//image_writer->Update();
	
}

template< typename TIPixel, typename TLPixel, unsigned int VImageDimension > 
void Synapse_Seg< TIPixel, TLPixel, VImageDimension>
::GetOutput()
{
	IntensityImageWriter::Pointer image_writer = IntensityImageWriter::New();
	image_writer->SetFileName("synapse.tif");
	image_writer->SetInput(SynOutImage);
	image_writer->Update();
	//return this->SynOutImage;
}


#endif // SYNAPSE_SEG_txx