#ifndef SYNAPSE_SEG_H
#define SYNAPSE_SEG_H

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <utility>

#include <QtGui/QDialog>
#include <QtGui/QFileDialog>
#include <QString>

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkImageDuplicator.h"


#include "ftkCommon/ftkUtils.h"

struct Params 
{
	double nuc_mask_threshold;
	double association_intensity_threshold;
	double association_distance_threshold;
	int mask_size;
	double mask_sigma;
	double syn_threshold;

};

struct Param{ std::string name; double value; };

template< typename TIPixel = unsigned char, typename TLPixel = unsigned short, unsigned int VImageDimension = 2> 
class Synapse_Seg
{
public:

	Synapse_Seg();
	~Synapse_Seg(){};

	typedef TIPixel IntensityPixelType;
	typedef TLPixel LabelPixelType;

	typedef itk::Image< IntensityPixelType, VImageDimension > IntensityImageType;
	typedef itk::Image< LabelPixelType, VImageDimension > LabelImageType;

	typedef itk::ImageRegionIterator< IntensityImageType> intIteratorType;
	typedef itk::ImageRegionIterator< LabelImageType> labIteratorType;

	typedef typename IntensityImageType::Pointer IntensityImagePointer;
	typedef typename LabelImageType::Pointer LabelImagePointer;

	typedef itk::ImageFileWriter< IntensityImageType > IntensityImageWriter;
	typedef itk::ImageFileWriter< LabelImageType > LabelImageWriter;

	typedef itk::ImageDuplicator< IntensityImageType > IntensityDuplicatorType;

	bool SetInput( IntensityImagePointer nucImgIn, IntensityImagePointer assImgIn, IntensityImagePointer synImgIn,LabelImagePointer lblImgIn);
	void SetParameters(std::vector<Param > parameters);
	void Update();
	void GetOutput();

protected:

private:

	IntensityImagePointer nucImage;	    //Input nuclei image;
	IntensityImagePointer assImage;	//Input associate image;
	IntensityImagePointer synImage;	    //Input synapse image;
	LabelImagePointer labelImage;			//Input label image;
	IntensityImagePointer SynOutImage;	    //Outputput synapse image;

	long int num_colum;
	long int num_row;
	Params paras;

	IntensityImagePointer nuc_mask;	    // nuclei mask;
	IntensityImagePointer assmaskedImage;	//masked associated image;
	IntensityImagePointer masked_syn;	
	IntensityImagePointer mask;	    // mask used to masking synapse with local regions
	int num_label;
	int num_tar;
	int num_nontar;
	int num_region;

	std::vector< std::vector< int> > centroids;
	std::vector<double > areas;
	std::vector<double > associate_intensity;
	std::vector<int > tar_label;
	std::vector<int > non_tar_label;
	std::vector<int > labels;
	std::vector<std::vector< double> > distance_mat;
	std::vector<std::vector< double> > centers;
	std::vector<std::pair<int,int > > syn_pairs;

	void ComputeLabelStatistics();
	void CreatNucMask();
	void ComputeMaskedAssociateImage();
	void ComputeAssociateIntensity();
	void ClassifyCells();
	void ComputeDistanceTarNonTar();
	void GetTarNeighbour();
	void CreatGlobalMask();
	void ThresholdSynImage();
	void Masking();

	QString lastPath;	

};     // end class Synapse_Seg

#include "Synapse_Seg.txx"

#endif // end SYNAPSE_SEG_H


