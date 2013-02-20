/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef _MODEL_NUCLEUS_SEG_H_
#define _MODEL_NUCLEUS_SEG_H_

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <itkSignedDanielssonDistanceMapImageFilter.h>


#include <tinyxml/tinyxml.h>
#include <ftkCommon/ftkUtils.h>
#include <ftkFeatures/ftkObjectAssociation.h>

#include "yousef_core/yousef_seg.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBoundingBox.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkNeighborhood.h"
#include "itkImageDuplicator.h"
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_matrix.h>
#include <vcl_vector.h>
#include <vcl_iostream.h>
#include <vcl_fstream.h>
#include <ftkFeatures/ftkLabelImageToFeatures.h>
#include <ftkCommon/ftkgnt.h>
#include "ftkIntrinsicFeatures.h"

#include <PatternAnalysis/embrex/kpls.h>
#include <PatternAnalysis/libsvm/svm.h>
#include <float.h>
#include "itkPasteImageFilter.h"

#include <iostream>
#include <cmath>
//#include "itkListSample2.h"
#include "itkFixedArray.h"

//#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_determinant.h>
#include "itkVariableSizeMatrix.h"
#include <tinyxml/tinyxml.h>
#include <itkSmartPointer.h>

#include <limits.h>
#include "NuclearSegmentation/Nuclear_Association/ftkNuclearAssociationRules.h"
#include "CytoplasmSegmentation/whole_cell.h"
#include "itkTileImageFilter.h"

#include <ftkObject.h>
#include "ftkImage.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRegionOfInterestImageFilter.h"


#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <vtkTable.h>
#include <itkMedianImageFilter.h>
#include <PatternAnalysis/agf/agf.h>
#include "itkOtsuMultipleThresholdsCalculator.h"


#include <map>

using namespace std;
using namespace boost;


#define MM_PI		3.14159265358979323846
#define MM_PI_2		1.57079632679489661923

class model_nucleus_seg 
{
public:
  
  static const int myDimension = 3;
  typedef itk::Image< unsigned char, myDimension>InputImageType;
  typedef itk::Image< unsigned short, myDimension > OutputImageType;
  typedef itk::Image< double, myDimension > DistImageType;
  typedef itk::Image< unsigned short int, 2 > UShortImageType;
  
  typedef OutputImageType::RegionType RegionType;
  typedef InputImageType::RegionType IRegionType;


  typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelComponentType;
  typedef itk::LabelStatisticsImageFilter< InputImageType, OutputImageType > LabelStatsFilterType;  
  typedef ftk::LabelImageToFeatures< unsigned char,  unsigned short, myDimension > FeatureCalcType;
  typedef itk::PasteImageFilter<OutputImageType> PasteFilterType;
  typedef itk::SignedDanielssonDistanceMapImageFilter<DistImageType, DistImageType > DTFilter ;
  typedef itk::ImageDuplicator< OutputImageType > DuplicatorType;

  typedef itk::VariableSizeMatrix<double>  CovarianceType;
  typedef itk::ImageRegionIterator< InputImageType > IIteratorType;
  typedef itk::ImageRegionIterator< OutputImageType> IteratorType;
  typedef itk::ImageRegionIterator< DistImageType> DIteratorType;
  typedef itk::ConstNeighborhoodIterator< OutputImageType > NeighborhoodIteratorType;

  typedef itk::MedianImageFilter< InputImageType, InputImageType > MedianFilterType;


  typedef itk::BoundingBox<unsigned short, myDimension, double> BB; 
  typedef ftk::IntrinsicFeatures FeaturesType;	
  typedef struct { int number; std::string name; std::string type; } Channel;	
	
  InputImageType::Pointer inputImage;
  OutputImageType::Pointer bImage;	
  OutputImageType::Pointer bImage2;

   int featureNumb; 
	//: constructor
	model_nucleus_seg();
 	//: destructor
	~model_nucleus_seg();



	vnl_matrix <double> normTrSet(char* train_fname,int inputdDimensions);
	vnl_matrix <double> getTrSet(char* train_fname,int inputdDimensions);

	void SetRawImage(const char* fname);	
	void SetImages(InputImageType::Pointer rawImage, OutputImageType::Pointer labelImage);
	void SetSplitImage(OutputImageType::Pointer img1);
	OutputImageType::Pointer SplitImage(std::vector<unsigned short> outliers,InputImageType::Pointer rimage, OutputImageType::Pointer limage);
	void SetTrainingFile(char* fname, int inputDimensions);				
	std::vector<unsigned short> getListofIds();
	std::vector<int> getNucIdsShattered();
	inline std::vector< unsigned short > getLabelIndex() {return labelIndex;};
	inline std::vector<FeaturesType> getFeats() {return allFeat;};
	inline int getNumFeat() {return NUM_FEAT;};
	inline vnl_matrix<double> getnormFeats() { return normFeats;};
	inline vnl_matrix<double> gettrainFeats() { return Feats;};
	inline InputImageType::Pointer getRawImage() { return inputImage;};
	inline OutputImageType::Pointer getLabImage() { return bImage;};

	std::vector<OutputImageType> getIndLabelImages(){ };
	std::vector<InputImageType> getIndRawImages(){ };
	void GetFeatsnImages();
	//void GetScoresfromKPLS(ftkgnt::MTreeType mTree,KPLS *kpls);
	void SetAssocFeatNumber(int nFeat);
	void PerformPCA();
	void Tset2EigenSpace(vnl_matrix<double> nFeats,vnl_matrix<double> Feats);
	inline int getasscofeatnumb() { return NUM_ASSOC_FEAT;};
		

	void Associations(char* xmlImage,char* projDef);
	void Associations_NucEd(ftk::Image::Pointer xmlImage, std::vector<ftk::AssociationRule> assocRules);
	std::vector < std::vector<double> > ComputeAssociations(void);
	void splitAssociations();
	OutputImageType::Pointer Shatter_Nuclei(InputImageType::Pointer roiImage,unsigned long * maxLabel);
	
	void SelectHypothesis();
	double CalcGenModelScore(double** train, long ntrain,long nvar,long ntest,double **test);
	OutputImageType::Pointer SplitImage(std::vector<int> outliers,InputImageType::Pointer rimage, OutputImageType::Pointer limage);
	std::vector<unsigned short> runSVM(vnl_matrix<double> feats,model_nucleus_seg::FeatureCalcType* filter);
	std::vector<unsigned short> Detect_undersegmented_cells();
	bool LoadAssoc(std::string filename);
	FeaturesType get_merged_features(set<int> currRPS);
	void getFeatureVectorsFarsight(OutputImageType::Pointer im, InputImageType::Pointer in_image, std::vector<FeaturesType> & feature_vector);
	OutputImageType::Pointer extract_label_image(int label, float bbox[6],OutputImageType::Pointer l);
	InputImageType::Pointer extract_raw_image(float bbox[6],InputImageType::Pointer r);
	void getFeatureVectorsFarsight_merge_test(OutputImageType::Pointer im, InputImageType::Pointer in_image, std::vector<FeaturesType> & feature_vector);
	std::vector<double> Convert2EigenSpace(std::vector<double> feat);
	void labelChange(unsigned short id1,unsigned short id2) ;
	std::vector< std::vector<double> > GetFeaturesOriginal();
	std::vector<double>  get_merged_assoc_features(set<int> currRPS,float* bbox);
	void PerformMerges(const char* x);
	OutputImageType::Pointer PerformMerges_NucEd();
	double GetScoreforId(std::vector<double> filtered_features);
	std::vector<double> GetOriginalScores();
	void UpdateBoolMerge();
	void DispScores(std::vector<double> scores);
	


	std::vector<ftk::AssociationRule> ReadAssociationRules(TiXmlElement * inputElement);
	std::vector<double> ComputeOneAssocMeasurement(itk::SmartPointer<OutputImageType> trgIm, int ruleID, std::vector<int>objID,float* b1);
	bool sorthelp (double i,double j); 
	bool uniquehelp (unsigned short i, unsigned short j);
	void norm_vec_std(double **mat1, long m, long n1, double **mat2, long n2);
	int Num_Lines2(char *root, int mode);
	FILE* FDeclare2(char *root, char *extension, char key);
	void calc_norm(double **mat, long D,long n, double *ave,double *stdev);
	std::string convert2string(unsigned long id);
	std::vector<double **> PrepareDataforScores();
	void SetTrainingFile(char* fname);
	double ComputeTotal(std::vector<int> LST);
	double ComputeAverage(std::vector<int> LST);
	//ftkgnt::MTreeType BuildMergeTreeDcon(ftkgnt::RAGraph,unsigned short id,std::vector< std::set<int> > hypothesis);
	void GetScoresfromKPLS(ftkgnt::MTreeType mTree);
	void getFeatureNames(char* fname);
	std::vector<double> populateVector(ftk::IntrinsicFeatures *f);
	//void calcvolLimits(vnl_matrix<double> feats);
	std::vector<double> calcStats(std::vector<double> classVol);
	unsigned short returnthresh( OutputImageType::Pointer input_image, int num_bin_levs, int num_in_fg );
	bool LoadSegParams(std::string filename);
	int GetYousefSeg();
	void RemoveSmallComponents(int minObjSize,const char* x);	

	double MAX_VOL, MIN_VOL;
	int MAX_DEPTH;
	double WC_DEFAULT;
	int SPLIT;

	int splitflag;
	
private:	

	std::vector<Channel> inputs;	
	int globalcount;
	std::vector <double> originalScores;
	std::vector <double> mhmRows; // multiple hypothesis matrix rows
	std::vector <double> mhmBRows; // multiple hypothesis binary matrix rows
	std::vector< std::vector<double> > hypoMatrix; // hypothesis matrix
	std::vector< std::vector<double> > hypoBMatrix; // hypothesis matrix
	std::vector <double> mergeindices;	
	std::vector< std::set<int> > hypothesis;
	unsigned short numLabels;
	std::vector<double **> train4class;
	std::vector<double> prior;
	double **PCAdata;	// Stores the training set after reducing the number of dimensions
	double *myKnownClass;
	 

	std::vector<OutputImageType::Pointer> li;
	std::vector<InputImageType::Pointer> ri;
	std::vector<FeaturesType> allFeat;
	std::vector< FeatureCalcType::LabelPixelType > labels;
	std::vector<ftk::AssociationRule> associationRules;	
	std::vector< std::vector<int> > classIndex;
	std::vector<std::string> featureNames;
	vnl_matrix<double> Feats;
	vnl_matrix<double> normFeats;

	int NUM_ASSOC_FEAT;
	int NUM_COMP;
	int NUM_FEAT;
	int NUM_CLASS;
	int SPLIT_SENSITIVITY;

	
	
	int minScale;
	int maxScale;

	unsigned short maxlabel;
	double medianScore;
	double *std1, *ave;
	double cFlag;
	FeatureCalcType::Pointer labFilter;
    FeatureCalcType::Pointer labFilter_b;
	int classNumb;
    std::vector< double > volLimit;
	std::vector<double> volumesys;
	std::vector<unsigned short> boolMergeTree;
	std::vector< vnl_vector<double> > eigenvecs; // eigenvector matrix
	int numSplit;
    int preparedatacounter;
	ftk::NuclearAssociationRules *assoc;
	std::vector< std::vector< double> > AssocFeat; 
	std::vector< std::vector< double > > objfeaturesNN; 
	
	unsigned short maxL;
	std::vector< unsigned short > labelIndex;

	double* trMean;
	double* trSTD;
	double totalCount;
	
	map<unsigned short,unsigned short> mymap;
	unsigned short thresh;
	int x_Size;
	int y_Size;
	int z_Size;
	
	char* myfilenamenorm;
	char* myfilenamePCA;
	

  DuplicatorType::Pointer duplicator; 
  //Declare the images used
  OutputImageType::Pointer inpImage ;		
  OutputImageType::Pointer bImagetemp;  
  OutputImageType::Pointer dmImage;
  DistImageType::Pointer GMimage  ;
  OutputImageType::IndexType pixelIndex;
  OutputImageType::Pointer clonedImage;
  OutputImageType::Pointer finalImage;
  ftk::Image::Pointer myFtkImage;
  char * myfilename;
  std::vector<unsigned short> imgindex;

};

#endif

