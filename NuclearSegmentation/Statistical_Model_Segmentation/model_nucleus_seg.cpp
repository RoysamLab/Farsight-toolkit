
#include "model_nucleus_seg.h"
#include "glpk.h"
#include <stdlib.h>
#include <time.h>
#include <float.h>
//#include <conio.h>


using namespace std;


template <typename T>
typename T::Pointer returnROI(typename T::Pointer image,typename T::RegionType region1)
{

	typedef itk::RegionOfInterestImageFilter< T,T> ROIFilterType;
	typename ROIFilterType::Pointer roifilter = ROIFilterType::New();
	roifilter->SetInput(image);
	roifilter->SetRegionOfInterest(region1); 	

	try
	{
		roifilter->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught!" <<std::endl;
		std::cout << err << std::endl;
		//return EXIT_FAILURE;
	}
	printf("Done.\n");
	return roifilter->GetOutput();

}

template <typename T>
typename T::Pointer readImage(const char *filename)
{
	printf("Reading %s ... ",filename);
	typedef typename itk::ImageFileReader<T> ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();

	ReaderType::GlobalWarningDisplayOff();
	reader->SetFileName(filename);
	try
	{
		reader->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught!" <<std::endl;
		std::cout << err << std::endl;
		//return EXIT_FAILURE;
	}
	printf("Done.\n");
	return reader->GetOutput();
}



template <typename T>
int writeImage(typename T::Pointer im, const char* filename)
{
	printf("Writing %s ... ",filename);
	typedef typename itk::ImageFileWriter<T> WriterType;

	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(im);
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught!" <<std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}
	printf("Done.\n");
	return EXIT_SUCCESS;
}



model_nucleus_seg::model_nucleus_seg()
{

	inputImage = NULL;
	bImage  = NULL;
	bImage2 = NULL;
	duplicator = NULL; 
	bImagetemp= NULL;  
	clonedImage = NULL;
	featureNumb = NULL;
	PCAdata = NULL;

	myfilenamenorm ="";
	myfilenamePCA ="PCAfile.txt";

	objfeaturesNN.clear();
	classIndex.clear();
	myFtkImage =NULL;
	inputs.clear();
	hypothesis.clear();
	eigenvecs.clear();
	imgindex.clear();
	AssocFeat.clear();
	Feats.clear();
	normFeats.clear();
	
	MAX_VOL = 1e5;
	MAX_DEPTH = 6;
	WC_DEFAULT = 5.; //WC_DEFAULT is used in calc scores: libagf parameter
	SPLIT = 0;
	globalcount = 0;
	splitflag = 0;
	maxL = 0;
}


///////////////////////
// SET THE RAW IMAGE
///////////////////////
void model_nucleus_seg::SetRawImage(const char* fname)
{
	inputImage = readImage<InputImageType>(fname);

	////Get gradient image from cytoplasm image
	MedianFilterType::Pointer  mfilter = MedianFilterType::New();
	InputImageType::SizeType radius; 
	radius[0] = 2; // radius along x 
	radius[1] = 2; // radius along y 
	radius[2] = 1; // radius along y 
	mfilter->SetRadius(radius);
	mfilter->SetInput( inputImage );
	mfilter->Update();
	inputImage = mfilter->GetOutput();
}		

void model_nucleus_seg::SetAssocFeatNumber(int nFeat)
{
	NUM_ASSOC_FEAT = nFeat;
	featureNumb = this->NUM_ASSOC_FEAT + this->NUM_FEAT;
}	


/////////////////////////
// SET IMAGES
/////////////////////////
void model_nucleus_seg::SetImages(InputImageType::Pointer rawImage, OutputImageType::Pointer labelImage)
{
	inputImage = rawImage;
	bImage = labelImage;
}


/////////////////////////
// SET THE TRAINING FILE
/////////////////////////
void model_nucleus_seg::SetTrainingFile(char* fname)
{	

	int inputDimensions = this->NUM_FEAT+this->NUM_ASSOC_FEAT;

	Feats = this->getTrSet(fname,inputDimensions); //Samples * Feats+2 matrix
	
	// Calculate the (Avg Volume - 1*Std of vol ) for each class
	//This will be useful for constructing trees.
	//calcvolLimits(Feats);
	
	//Normalize the features
	normFeats = this->normTrSet(fname,inputDimensions); 
	
	//Compute the sample covariance matrix
	vnl_matrix<double> normFeats_tr = normFeats.transpose();// Feats * Samples matrix
	vnl_matrix<double> covMat = normFeats_tr*normFeats/(normFeats.rows()); // Feats * Feats matrix
	
	//Compute the eigen vectors of the covariance matrix
	vnl_symmetric_eigensystem<double> eig(covMat);

	std::vector<double> currvec;
	currvec.resize(inputDimensions);
	eigenvecs.resize(inputDimensions);


	for(int i = 0; i < inputDimensions; i++)
	{
		eigenvecs[i] = eig.get_eigenvector(i);
		//std::cout<< eigenvecs[i] << std::endl;
		//std::cout<< eig.get_eigenvalue(i) << std::endl;
	}
	Tset2EigenSpace(normFeats,Feats);
	getFeatureNames(fname); //get the feature names used in this model
}		


//void model_nucleus_seg::calcvolLimits(vnl_matrix<double> feats)
//{
//
//	std::vector< std::vector<double> > vols;
//	
//	//Get the volumes of each class into a vector
//	for(int classval = 1; classval<NUM_CLASS+1 ; classval++) //NUM_CLASS
//	{
//		std::vector<double> classVol;
//		for(unsigned int r=0; r<feats.rows(); ++r)
//		{	
//			if(feats(r,NUM_FEAT+NUM_ASSOC_FEAT)==classval)
//				classVol.push_back(feats(r,0));
//		}
//		vols.push_back(classVol);
//	}
//
//	
//	int ctr=0;
//	
//	//extract data from the model and get min/max values:
//	for(unsigned int r=0; r<vols.size(); ++r)
//	{	
//		std::vector<double> classVol = vols[r];
//		std::vector<double> stats = calcStats(classVol);
//		volLimit[r] = stats[0]-(1*stats[1]); 
//
//		if(volLimit[r] < 0 )
//			volLimit[r] = stats[0]; 
//	}
//}



// Calulates the mean and standard deviation of the elements of a vector
std::vector<double> model_nucleus_seg::calcStats(std::vector<double> classVol)
{
	double sum =0;
	std::vector<double> stats;
	stats.resize(NUM_CLASS);//NUM_CLASS
	for(int j =0;j<classVol.size();++j)
	{
		sum = sum+classVol[j];
	}

	double mu = sum/classVol.size();
	sum = 0;

	for(int j =0;j<classVol.size();++j)
	{
		double x = classVol[j];
		sum = sum + (x - mu)*(x -mu);
	}
	
	double std = sum/classVol.size();
	stats[0] = mu;
	stats[1] = std;
	return stats;
}

// Get the Feature names which are used in the model
void model_nucleus_seg::getFeatureNames(char* fname)
{
	vtkSmartPointer<vtkTable> model_table;
	model_table = ftk::LoadTable(fname);

	for(int d=0; d<(int)model_table->GetNumberOfColumns(); ++d)
	{
		std::string model_column = model_table->GetColumnName(d);
		this->featureNames.push_back(model_column);
	}
}


//Compute the assoications
void model_nucleus_seg::Associations(char* xmlImage,char* projDef )
{	
	allFeat.clear();
	labelIndex.clear();
	
	labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( inputImage, bImage );
	labFilter->SetLevel(3);
	labFilter->ComputeHistogramOn();
	labFilter->Update();
	labels = labFilter->GetLabels();
	
	//Any Id greater than maxL has been split !
	maxL = labels.size()-1;

	getFeatureVectorsFarsight(bImage,inputImage,allFeat);
	labelIndex.resize(allFeat.size());	
	
	for(int counter=0; counter < allFeat.size(); counter++)
	{
		labelIndex[counter] =  allFeat[counter].num;
	}	

	myFtkImage = ftk::LoadXMLImage(xmlImage);
	LoadAssoc(projDef);
	AssocFeat = ComputeAssociations();
}

void model_nucleus_seg::Associations_NucEd(ftk::Image::Pointer xmlImage, std::vector<ftk::AssociationRule> assocRules )
{	
	allFeat.clear();
	labelIndex.clear();
	
	labFilter = FeatureCalcType::New();
	labFilter->SetImageInputs( inputImage, bImage );
	labFilter->SetLevel(3);
	labFilter->ComputeHistogramOn();
	labFilter->Update();
	labels = labFilter->GetLabels();
	
	//Any Id greater than maxL has been split !
	maxL = labels.size()-1;

	getFeatureVectorsFarsight(bImage,inputImage,allFeat);
	labelIndex.resize(allFeat.size());	
	
	for(int counter=0; counter < allFeat.size(); counter++)
	{
		labelIndex[counter] =  allFeat[counter].num;
	}	

	//myFtkImage = ftk::LoadXMLImage(xmlImage);
	//LoadAssoc(projDef);
	myFtkImage = xmlImage;
	associationRules = assocRules;
	AssocFeat = ComputeAssociations();
}


// COmpute associations if the splitting has been performed
void model_nucleus_seg::splitAssociations()
{	
	AssocFeat.clear();
	AssocFeat = ComputeAssociations();
}


///////////////////////
// SET SPLIT IMAGE
///////////////////////
void model_nucleus_seg::SetSplitImage(OutputImageType::Pointer img1)
{	
	this->bImage  = NULL;
	this->bImage = img1;
}	


std::vector<unsigned short> model_nucleus_seg::Detect_undersegmented_cells()
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////	
	//RELABEL THE LABELS IN THE IMAGE TO MAKE THEM SERIALLY ORDERED
	//////////////////////////////////////////////////////////////////////////////////////////////////////	
	RelabelComponentType::Pointer relabel = RelabelComponentType::New();
	relabel->SetInput( bImage );
	relabel->Update();
	bImage = relabel->GetOutput();
	//////////////////////////////////////////////////////////////////////////////////////////////////////	
	

	///////////////////////////////////////////////////////////////////////////////////////////////////////	
	// Run SVM : DETECT THE OBJECTS THAT ARE UNDERSEGMENTED
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<unsigned short> outliers = runSVM(this->gettrainFeats(),labFilter);
	std::cout<< outliers.size() <<" ---- cells are undersegmented " <<std::endl;
	
	return outliers;
}	

std::vector<unsigned short> model_nucleus_seg::getListofIds()
{	
	return boolMergeTree;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE FEATURES OF THE MERGED OBJECT BASED ON THE RPS
// USES BOUNDING BOX INFORMATION
///////////////////////////////////////////////////////////////////////////////////////////////////////

model_nucleus_seg::FeaturesType model_nucleus_seg::get_merged_features(set<int> currRPS)
{
	set<int>::iterator RPSIterator;
	RPSIterator = currRPS.begin();	
	std::vector<model_nucleus_seg::OutputImageType::Pointer> sli;
	std::vector<model_nucleus_seg::InputImageType::Pointer> sri;
	std::vector< std::vector<int> > bounds; 	 

	sli.resize(currRPS.size());
	sri.resize(currRPS.size());
	bounds.resize(currRPS.size());

	int lbounds[6] = { 100000 , -100000, 100000 , -100000 , 100000 , -100000 };
	int ctr = 0;



///////////////////////////////////////////////////////////////////////////////////////////////////////
// Get the Bounding boxes

	for(; RPSIterator != currRPS.end(); ++RPSIterator)
	{
		vector<unsigned short>::iterator posn1 = find(labelIndex.begin(), labelIndex.end(), *RPSIterator);
		sli[ctr] = li[posn1 - labelIndex.begin()];
		sri[ctr] = ri[posn1 - labelIndex.begin()];	

		if(allFeat[posn1 - labelIndex.begin()].BoundingBox[0] < lbounds[0] )
		{
			lbounds[0] = allFeat[posn1 - labelIndex.begin()].BoundingBox[0];
		}


		if(allFeat[posn1 - labelIndex.begin()].BoundingBox[2] < lbounds[2] )
		{
			lbounds[2] = allFeat[posn1 - labelIndex.begin()].BoundingBox[2];
		}

		if(allFeat[posn1 - labelIndex.begin()].BoundingBox[4] < lbounds[4] )
		{
			lbounds[4] = allFeat[posn1 - labelIndex.begin()].BoundingBox[4];
		}

		if(allFeat[posn1 - labelIndex.begin()].BoundingBox[1] > lbounds[1] )
		{
			lbounds[1] = allFeat[posn1 - labelIndex.begin()].BoundingBox[1];
		}

		if(allFeat[posn1 - labelIndex.begin()].BoundingBox[3] > lbounds[3] )
		{
			lbounds[3] = allFeat[posn1 - labelIndex.begin()].BoundingBox[3];
		}

		if(allFeat[posn1 - labelIndex.begin()].BoundingBox[5] > lbounds[5] )
		{
			lbounds[5] = allFeat[posn1 - labelIndex.begin()].BoundingBox[5];
		}

		++ctr;
	}	  
///////////////////////////////////////////////////////////////////////////////////////////////////////

	OutputImageType::SizeType ls;
	ls[0] = lbounds[1]-lbounds[0]+1;
	ls[1] = lbounds[3]-lbounds[2]+1;
	ls[2] = lbounds[5]-lbounds[4]+1;

	OutputImageType::Pointer p = OutputImageType::New();
	InputImageType::Pointer r = InputImageType::New();
	OutputImageType::IndexType lindex;
	lindex.Fill(0);
	OutputImageType::RegionType lregion;
	lregion.SetIndex(lindex);
	lregion.SetSize(ls);

	p->SetRegions(lregion);
	p->Allocate();
	p->FillBuffer(0);
	r->SetRegions(lregion);
	r->Allocate();
	r->FillBuffer(0);

	RPSIterator = currRPS.begin();

	int counter = 0;
	int ctr1 = 0;
	std::vector<int> currlabs;
	currlabs.resize(currRPS.size());

	//std::cout<<"Merged ---> " << std::endl;
	RPSIterator = currRPS.begin();
	for(; RPSIterator != currRPS.end(); ++RPSIterator)
	{
		vector<unsigned short>::iterator posn1 = find(labelIndex.begin(), labelIndex.end(), *RPSIterator);		
		currlabs[ctr1] = allFeat[posn1 - labelIndex.begin()].num;
		ctr1++;
	}


	RPSIterator = currRPS.begin();

	for(; RPSIterator != currRPS.end(); ++RPSIterator)
	{

		OutputImageType::Pointer lbl;
		InputImageType::Pointer raw;		
		std::vector<unsigned short>::iterator posn1 = find(labelIndex.begin(), labelIndex.end(), *RPSIterator);		
		
		lbl = sli[counter];
		raw = sri[counter];

		IteratorType liter1(lbl,lbl->GetLargestPossibleRegion());
		IIteratorType riter1(raw,raw->GetLargestPossibleRegion());

		lindex[0] = allFeat[posn1 - labelIndex.begin()].BoundingBox[0]-lbounds[0];
		lindex[1] = allFeat[posn1 - labelIndex.begin()].BoundingBox[2]-lbounds[2];
		lindex[2] = allFeat[posn1 - labelIndex.begin()].BoundingBox[4]-lbounds[4];		

		lregion.SetSize(lbl->GetLargestPossibleRegion().GetSize());
		lregion.SetIndex(lindex);

		IteratorType liter(p,lregion);
		IIteratorType riter(r,lregion);		
		for(liter1.GoToBegin(),riter1.GoToBegin(),liter.GoToBegin(),riter.GoToBegin();!liter1.IsAtEnd(); ++liter1,++riter1,++liter,++riter)
		{
			if(currlabs.end() != find(currlabs.begin(), currlabs.end(), liter1.Get()))	
			{
				liter.Set(255);	
			}
			else 
			{
				liter.Set(liter1.Get());			
			}	
			riter.Set(riter1.Get());
		}		
		++counter;	
	}	  

	std::vector<FeaturesType> f1;
	model_nucleus_seg::getFeatureVectorsFarsight_merge_test(p,r,f1);
	
	FeaturesType f = f1[0];
	f.Centroid[0]+=lbounds[0];
	f.Centroid[1]+=lbounds[2];
	f.Centroid[2]+=lbounds[4];
	f.WeightedCentroid[0]+=lbounds[0];
	f.WeightedCentroid[2]+=lbounds[2];
	f.WeightedCentroid[4]+=lbounds[4];
	f.BoundingBox[0]+=lbounds[0];
	f.BoundingBox[1]+=lbounds[0];
	f.BoundingBox[2]+=lbounds[2];
	f.BoundingBox[3]+=lbounds[2];
	f.BoundingBox[4]+=lbounds[4];
	f.BoundingBox[5]+=lbounds[5];


	IteratorType liter(p,p->GetLargestPossibleRegion());

	for(liter.GoToBegin();!liter.IsAtEnd(); ++liter)
	{
		if(liter.Get()!=255)	
		{
			liter.Set(0);
		}
	}		
	return f;
}		

///////////////////////////////////////////////////////////////////////////////////////////////////////
// STORE ALL INDIVIDUAL SEGMENTATIONS AND RAW FILES IN SEPARATE IMAGES 
// ALSO GET THE FEATURES OF OBJECTS IN "allFeat" vector
///////////////////////////////////////////////////////////////////////////////////////////////////////
void model_nucleus_seg::GetFeatsnImages()
{
	
	allFeat.clear();
	labelIndex.clear();
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// NEED TO CALCULATE THE FEATURES OF THE NEW SPLIT IMAGE WITH A NEW FILTER 
	///////////////////////////////////////////////////////////////////////////////////////////////////////	
	labFilter = FeatureCalcType::New();
	labFilter->SetCompleteImageInputs( inputImage, bImage);
	labFilter->SetLevel(3);
	labFilter->ComputeHistogramOn();
	labFilter->Update();
	labels = labFilter->GetLabels();
	
	getFeatureVectorsFarsight(bImage,inputImage,allFeat);
	labelIndex.resize(allFeat.size());	
	
	for(int counter=0; counter < allFeat.size(); counter++)
	{	
		li.push_back(model_nucleus_seg::extract_label_image(allFeat[counter].num,allFeat[counter].BoundingBox,bImage));
		ri.push_back(model_nucleus_seg::extract_raw_image(allFeat[counter].BoundingBox,inputImage));		 
		labelIndex[counter] =  allFeat[counter].num;
	}	

	//Allocate memory for the Multiple Hypothesis Matrix rows
	// depends on the number of fragments : allFeat
	mhmRows.resize(allFeat.size());
	mhmBRows.resize(allFeat.size());
	totalCount = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE FEATURES AND INSERT IN "allFeat"
///////////////////////////////////////////////////////////////////////////////////////////////////////

void model_nucleus_seg::getFeatureVectorsFarsight(OutputImageType::Pointer im, InputImageType::Pointer in_image, std::vector<FeaturesType> & feature_vector)
{
	for(unsigned int counter=0; counter< labels.size(); counter++)
	{
		if(labels[counter]==0)
			continue;
		feature_vector.push_back(*(labFilter->GetFeatures(labels[counter])));
		feature_vector.back().num = labels[counter];			
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE FEATURES OF THE MERGED OBJECT
///////////////////////////////////////////////////////////////////////////////////////////////////////
void model_nucleus_seg::getFeatureVectorsFarsight_merge_test(OutputImageType::Pointer im, InputImageType::Pointer in_image, std::vector<FeaturesType> & feature_vector)
{
	typedef ftk::LabelImageToFeatures<InputImageType::PixelType,OutputImageType::PixelType,myDimension> FeatureCalculatorType;
	FeatureCalculatorType::Pointer fc = FeatureCalculatorType::New();
	fc->SetImageInputs(in_image,im);
	fc->SetLevel(3);
	//fc->ComputeTexturesOn();
	fc->ComputeHistogramOn();
	fc->Update();

	std::vector<OutputImageType::PixelType> labels = fc->GetLabels();
	for(unsigned int counter=0; counter< labels.size(); counter++)
	{
		if(labels[counter]!=255)
			continue;
		feature_vector.push_back(*(fc->GetFeatures(labels[counter])));
		feature_vector.back().num = labels[counter];			
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// EXTRACT LABEL IMAGES AND STORE THEM INDIVIDUALLY: SPEEDS UP MERGING CALCULATIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////
model_nucleus_seg::OutputImageType::Pointer model_nucleus_seg::extract_label_image(int label, float bbox[6],OutputImageType::Pointer l)
{
	int bb[6];
	for(int co = 0; co < 6; co++)
	{
		bb[co] = int(bbox[co]+0.5);
	}
	OutputImageType::Pointer lp = OutputImageType::New();
	OutputImageType::IndexType lindex;
	OutputImageType::SizeType lsize;
	OutputImageType::RegionType lregion;
	lindex.Fill(0);
	lsize[0] = bb[1]-bb[0]+1;
	lsize[1] = bb[3]-bb[2]+1;
	lsize[2] = bb[5]-bb[4]+1;
	lregion.SetIndex(lindex);
	lregion.SetSize(lsize);
	lp->SetRegions(lregion);
	lp->Allocate();
	lp->FillBuffer(0);
	IteratorType lpiter(lp,lp->GetLargestPossibleRegion());
	lindex[0] = bb[0];
	lindex[1] = bb[2];
	lindex[2] = bb[4];
	lregion.SetIndex(lindex);
	IteratorType liter(l,lregion);
	for(;!lpiter.IsAtEnd(); ++lpiter,++liter)
	{
		//if(liter.Get()==label)
		lpiter.Set(liter.Get());
	}

	return lp;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// EXTRACT RAW IMAGES AND STORE THEM INDIVIDUALLY : SPEEDS UP MERGING CALCULATIONS
///////////////////////////////////////////////////////////////////////////////////////////////////////
model_nucleus_seg::InputImageType::Pointer model_nucleus_seg::extract_raw_image(float bbox[6],InputImageType::Pointer r)
{
	int bb[6];
	for(int co = 0; co < 6; co++)
	{
		bb[co] = int(bbox[co]+0.5);
	}
	InputImageType::Pointer lp = InputImageType::New();
	InputImageType::IndexType lindex;
	InputImageType::SizeType lsize;
	InputImageType::RegionType lregion;
	lindex.Fill(0);
	lsize[0] = bb[1]-bb[0]+1;
	lsize[1] = bb[3]-bb[2]+1;
	lsize[2] = bb[5]-bb[4]+1;
	lregion.SetIndex(lindex);
	lregion.SetSize(lsize);
	lp->SetRegions(lregion);
	lp->Allocate();
	lp->FillBuffer(0);
	IIteratorType lpiter(lp,lp->GetLargestPossibleRegion());
	lindex[0] = bb[0];
	lindex[1] = bb[2];
	lindex[2] = bb[4];
	lregion.SetIndex(lindex);
	IIteratorType liter(r,lregion);
	for(;!lpiter.IsAtEnd(); ++lpiter,++liter)
	{
		lpiter.Set(liter.Get());
	}
	return lp;
}





///////////////////////////////////////////////////////////////////////////////////////////////////////
//GET THE TRAINING SET INTO A MATRIX
///////////////////////////////////////////////////////////////////////////////////////////////////////

vnl_matrix <double> model_nucleus_seg::getTrSet(char* train_fname,int inputdDimensions)
{
	//Read Training File:
	//FILE * fpin = FDeclare2(train_fname, "", 'r');
	FILE *fpin;
	fpin = fopen(train_fname,"r");

	double val;
	char str [80];

	int myrows = Num_Lines2(train_fname, 0) - 1 ;   // number of rows            
	int mycols = inputdDimensions ;  // number of columns

	//First column contains string names
	for(int c=0; c< mycols+2; ++c)
	{
		//Edit -- read from file
		fscanf(fpin, "%s",str);
	}

	vnl_matrix <double> FeatsMatrix(myrows,mycols+2);
	double **trVecs = (double **) malloc(myrows*sizeof(double *));

	//extract data from the model and get min/max values:
	for(unsigned int r=0; r < myrows ; ++r)
	{
		trVecs[r] = (double *)malloc((NUM_FEAT+NUM_ASSOC_FEAT) * sizeof(double));
		for(int c=0; c< mycols; ++c)
		{
			//Edit -- read from file
				fscanf(fpin, "%lf", &val);
			trVecs[r][c] = val;
			FeatsMatrix.put(r,c,val);
		}
		fscanf(fpin, "%lf", &val);
		FeatsMatrix.put(r,inputdDimensions,val);
		fscanf(fpin, "%lf", &val);
		FeatsMatrix.put(r,inputdDimensions+1,val);
	}

	fclose(fpin);

	return FeatsMatrix;
}



vnl_matrix <double> model_nucleus_seg::normTrSet(char* train_fname,int inputdDimensions)
{
	//Read Training File:
	//FILE * fpin = FDeclare2(train_fname, "", 'r');
	FILE *fpin;
	fpin = fopen(train_fname,"r");

	double val;
	char str [80];

	int myrows = Num_Lines2(train_fname, 0) - 1 ;   // number of rows            
	int mycols = inputdDimensions ;  // number of columns
	double **trVecs = (double **) malloc(myrows*sizeof(double *));

	//First column contains string names
	for(int c=0; c< mycols+2; ++c)
	{
		//Edit -- read from file
		fscanf(fpin, "%s",str);
	}

	//extract data from the file
	for(unsigned int r=0; r < myrows ; ++r)
	{
		trVecs[r] = (double *)malloc((inputdDimensions) * sizeof(double));
		for(int c=0; c< mycols; ++c)
		{
			//Edit -- read from file
			fscanf(fpin, "%lf", &val);
			trVecs[r][c] = val;
			//std::cout<<val<<"   ";
		}
		//std::cout<<""<<std::endl;
		fscanf(fpin, "%lf", &val);
		fscanf(fpin, "%lf", &val);
	}

	fclose(fpin);


	trMean=new double[inputdDimensions];
	trSTD =new double[inputdDimensions];
	int ctr =0;	
	//Get the Mean and StD for each feature
	calc_norm(trVecs,inputdDimensions,myrows,trMean, trSTD);

	std::vector< std::vector<double> > newFeat;
	newFeat.resize(myrows);
	std::vector<double> newRow;	
	newRow.resize(mycols);
	
	
	vnl_matrix <double> normMatrix(myrows,mycols);

	for(unsigned int r=0; r < myrows ; ++r)
	{
		for(int c=0; c< mycols; ++c)
		{
			normMatrix.put(r,c,(trVecs[r][c]-trMean[c])/(trSTD[c]));
		}		
	}

	return normMatrix;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// TRANSFORM TRAINING FILE FEATURES FROM ORIGINAL SPACE TO EIGEN SPACE/PCA SPACE
///////////////////////////////////////////////////////////////////////////////////////////////////////
void model_nucleus_seg::Tset2EigenSpace(vnl_matrix<double> nFeats,vnl_matrix<double> Feats)
{
	std::vector<double> measure;
	measure.resize(NUM_FEAT+NUM_ASSOC_FEAT);


	FILE * fpin1 = FDeclare2(this->myfilenamePCA, "", 'w');
	const unsigned int numberOfrows =  nFeats.rows();		
	
	PCAdata = (double **) malloc(numberOfrows*sizeof(double *));
	myKnownClass  = (double *)malloc(numberOfrows * sizeof(double));
	for (int i = 0; i < numberOfrows; i++) 
	{	
		
		PCAdata[i] = (double *)malloc(NUM_COMP * sizeof(double));

		for (int j = 0; j < nFeats.columns(); j++) 
		{
			measure[j] = nFeats(i,j);	
		}

		std::vector<double> eigenMeasure = Convert2EigenSpace(measure);

		for(int c=0; c< this->NUM_COMP; ++c)
		{
			fprintf(fpin1, "%lf    ", eigenMeasure[c]);
			PCAdata[i][c] = eigenMeasure[c];
		}		

		//Then get class and Id:
		fprintf(fpin1, "%lf    ", Feats(i,Feats.columns()-2));
		myKnownClass[i] = Feats(i,Feats.columns()-2);
		fprintf(fpin1, "%lf    ", Feats(i,Feats.columns()-1));
		fprintf(fpin1,"\n");
	}
	fclose(fpin1);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// CONVERT FEATURES TO EIGENSPACE/PCA SPACE
///////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> model_nucleus_seg::Convert2EigenSpace(std::vector<double> feat)
{
	std::vector<double> eigenFeature;
	vnl_vector<double> currEig;
	eigenFeature.resize(NUM_COMP);

	for(int i = 0; i < NUM_COMP; i++)
	{
		// since eigenvecs stores the eigen vectors in ascending order 
		// of eigen values and since we need the last NUM_COMP(5) eigen vectors 
		currEig = eigenvecs[NUM_FEAT+NUM_ASSOC_FEAT-1-i];
		eigenFeature[i] = 0;

		for(int j = 0; j < NUM_FEAT+NUM_ASSOC_FEAT; j++)
		{	
			//std::cout<< currEig[j] << "-" << feat[j] <<std::endl;
			eigenFeature[i] += currEig[j]*feat[j]; 
		}	
		currEig.clear();
	}
	return eigenFeature;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE FEATURES OF THE INITIAL OBJECTS. USED IN GETORIGINALSCORES() METHOD
///////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector< std::vector<double> > model_nucleus_seg::GetFeaturesOriginal()
{
	std::vector< std::vector<double> > allrelFeat;

	allrelFeat.resize(allFeat.size());

	std::vector<double> in_features;	
	in_features.resize(NUM_FEAT+NUM_ASSOC_FEAT);
	int ctr2;
	
	for(int ctr =0 ; ctr < allFeat.size() ; ++ctr ) 
	{
		ftk::IntrinsicFeatures  features = allFeat[ctr];					
		std::vector<double> in_features1 = populateVector(&features);

		for(int i=0; i<NUM_FEAT;++i)	
		{
			in_features[i] = in_features1[i]; 
		}

		ctr2 = 0;

		for(int c=NUM_FEAT; c< NUM_FEAT+NUM_ASSOC_FEAT; ++c)
		{
			std::vector<double> aM =	AssocFeat[ctr];
			in_features[c] = aM[ctr2];//ctr2 corresponds to the Association rule... ctr refers to the label
			++ctr2;	
		}		
		
		allrelFeat[ctr] = in_features;
		in_features.clear();	
		in_features.resize(NUM_FEAT+NUM_ASSOC_FEAT);
	}	
	return allrelFeat;
}



std::vector<double> model_nucleus_seg::populateVector(ftk::IntrinsicFeatures *f)
{
	std::vector<double> filtered_features;
	filtered_features.resize(this->NUM_FEAT);
	int ctr = 0;
	for(int k =0; k < this->NUM_FEAT; ++k)
	{
		for(int i=0; i<32;++i)	
		{
			//Converts upper case to lower case
			std::transform(f->Info[i].name.begin(), f->Info[i].name.end(),f->Info[i].name.begin(),::tolower);	
			if(this->featureNames[k] ==  f->Info[i].name)
			{
				filtered_features[ctr] = f->ScalarFeatures[i];
				if(this->featureNames[k] == "convexity")
				{
					//Make sure convexity is not INF ( happens for really small fragments)
					if(filtered_features[ctr] > 10000000)
					{
						filtered_features[ctr] = 1;
					}
				}
				++ctr;
			}
		}
	}
	return(filtered_features);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// USED TO CHANGE THE LABEL OF ONE SEGMENT TO ANOTHER SEGMENT
// USED TO PERFORM THE FINAL MERGES
///////////////////////////////////////////////////////////////////////////////////////////////////////
void model_nucleus_seg::labelChange(unsigned short id1,unsigned short id2)  
{

	//Define a labelstatisticsimagefilter so that we can get the region corresponding to
	//"id2" that needs to be changed. Only this region's value needs to be changed to 
	//"id1"
	RegionType region1;

	LabelStatsFilterType::Pointer lsfilterChange = LabelStatsFilterType::New();
	lsfilterChange->SetInput(inputImage);
	lsfilterChange->SetLabelInput(clonedImage); 
	lsfilterChange->Update();
	region1 = lsfilterChange->GetRegion(id2);
	IteratorType regioniterator(clonedImage,region1);


	IteratorType regioniterator2(bImage,region1);
	
	for(regioniterator.GoToBegin(),regioniterator2.GoToBegin(); !regioniterator.IsAtEnd(),!regioniterator2.IsAtEnd();++regioniterator,++regioniterator2)
	{
		if(regioniterator2.Value() == id2) 
		{
			regioniterator.Set(id1);	
		}
	} 	

}	


///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE SCORES FOR ALL THE POSSIBLE MERGES IN A MERGETREE
///////////////////////////////////////////////////////////////////////////////////////////////////////
void model_nucleus_seg::GetScoresfromKPLS(ftkgnt::MTreeType mTree)
{ 
	unsigned int counter = 0; 
	set<int> currRPS;	

	//std::cout<< num_vertices(mTree) <<std::endl;
	while (counter != num_vertices(mTree))
	{
		currRPS =  mTree[counter].RPS;
		set<int>::iterator RPSIterator;
		RPSIterator = currRPS.begin();
		if(currRPS.size()!=1)
		{
			ftk::IntrinsicFeatures f = get_merged_features(currRPS);

			std::vector<double> filtered_features;
			int ctr =0;
			double score ;

			if(!NUM_ASSOC_FEAT)
			{
				filtered_features.resize(NUM_FEAT);
				filtered_features = populateVector(&f);
				//Get the score for this ID
				score = GetScoreforId(filtered_features);
			}
			else
			{
				std::vector<double> as = get_merged_assoc_features(currRPS,f.BoundingBox);
				filtered_features.resize(NUM_FEAT+NUM_ASSOC_FEAT);

				std::vector<double> in_features1 = populateVector(&f);

				for(int i=0; i<NUM_FEAT;++i)	
				{
					filtered_features[i] = in_features1[i]; 
				}

				ctr =0;	

				for(int c=NUM_FEAT; c< NUM_FEAT+NUM_ASSOC_FEAT; ++c)
				{
					filtered_features[c] = as[ctr];
					ctr++;	
				}		

				//Get the score for this ID	
				score = GetScoreforId(filtered_features);
			}

			for(RPSIterator = currRPS.begin(); RPSIterator != currRPS.end(); ++RPSIterator)
			{
				vector<unsigned short>::iterator posn1 = find(labelIndex.begin(), labelIndex.end(), *RPSIterator);
				mhmRows[posn1 - labelIndex.begin()] = (score - originalScores[posn1 - labelIndex.begin()]) ; 
				mhmBRows[posn1 - labelIndex.begin()] = 1 ; 
			}	 

			//Insert the row in the hypothesis matrix
			hypoMatrix.push_back(mhmRows);	
			hypoBMatrix.push_back(mhmBRows);

			//Make a note of all the hypothesis
			hypothesis.push_back(currRPS);		

			//reset the rows	
			mhmRows.clear();	
			mhmBRows.clear();	
			mhmRows.resize(allFeat.size());		 
			mhmBRows.resize(allFeat.size());
		}
		counter = counter +1;
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE ASSOCIATIVE FEATURES FOR THE MERGED OBJECT
///////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double>  model_nucleus_seg::get_merged_assoc_features(set<int> currRPS,float* bbox)
{
	set<int>::iterator RPSIterator;
	RPSIterator = currRPS.begin();

	std::vector<double> assocMerge;
	std::vector<double> aMval;

	std::vector<int> channelnumb;

	std::vector<int> idvec;
	idvec.resize(currRPS.size());

	//	
	int pos =0;
	// Put all the labels in idvec before feeding it to ComputeOneAssoc 
	for(; RPSIterator != currRPS.end(); ++RPSIterator)
	{
		vector<unsigned short>::iterator posn1 = find(labelIndex.begin(), labelIndex.end(), *RPSIterator);		
		idvec[pos] = (allFeat[posn1 - labelIndex.begin()].num);
		pos++;
	}

	int count=0;
	//For each of the associations
	for(std::vector<ftk::AssociationRule>::iterator ascit=associationRules.begin(); ascit!=associationRules.end(); ++ascit )
	{
		int seg_channel_number=-1;
		int inp_channel_number=-1;

		if( strcmp( ascit->GetSegmentationFileName().c_str(), "NUCLEAR" ) == 0 ){
			seg_channel_number = 0;

		} 
		else{
			std::cout<<"Please check region type for associative feature computation\n";
		}

		for( int j=0; j<(int)myFtkImage->GetImageInfo()->channelNames.size(); ++j )
			if( strcmp( myFtkImage->GetImageInfo()->channelNames.at(j).c_str(), ascit->GetTargetFileNmae().c_str() ) == 0 ){
				inp_channel_number=j;
				//std::cout<<"Channel-->"<<inp_channel_number<<std::endl;
				break;
			}
			if( inp_channel_number == -1 ){
				std::cout<<"Unable to access grayscale image while computing associative feature: "<<ascit->GetRuleName()<<std::endl;
			}


		if(channelnumb.end() == find(channelnumb.begin(), channelnumb.end(), inp_channel_number))
		{
			OutputImageType::Pointer inp  = OutputImageType::New();	
			inp = myFtkImage->GetItkPtr<unsigned short>(0,inp_channel_number,ftk::Image::DEEP_COPY);
			thresh	 = mymap[inp_channel_number];

			aMval = ComputeOneAssocMeasurement(inp, 0, idvec,bbox);
			for(int l1 =0;l1<aMval.size();l1++)
				assocMerge.push_back(aMval[l1]);
			count++;
			channelnumb.push_back(inp_channel_number);
		}
	}
	
	return assocMerge;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// PERFORM THE FINAL MERGES 
///////////////////////////////////////////////////////////////////////////////////////////////////////

void model_nucleus_seg::PerformMerges(const char* x)
{
	////Create a clone of the original output. Will perform merges on this image	
	typedef itk::ImageDuplicator< OutputImageType > DuplicatorType; 
	DuplicatorType::Pointer duplicator1 = DuplicatorType::New(); 
	duplicator1->SetInputImage(bImage); 
	duplicator1->Update();
	clonedImage = duplicator1->GetOutput();	

	std::cout<<"No.of valid hypotheses:"<<mergeindices.size()<<std::endl;

	for(int i =0 ; i< mergeindices.size() ; ++i)
	{
		set<int> currRPS;
		set<int>::iterator RPSIterator;
		currRPS =  hypothesis[mergeindices[i]];
		RPSIterator = currRPS.begin();
		unsigned short root = *RPSIterator;

		std::cout<<"Merged"<<"--->"<<root<<"-";

		++RPSIterator;
		for(; RPSIterator != currRPS.end(); ++RPSIterator)
		{
			labelChange(root,*RPSIterator);
			std::cout<<*RPSIterator<<"-";
		}

		std::cout<<""<<std::endl;		
	}	
	RemoveSmallComponents(MIN_VOL, x);
	//writeImage<OutputImageType>(clonedImage,x);
}

model_nucleus_seg::OutputImageType::Pointer model_nucleus_seg::PerformMerges_NucEd()
{
	////Create a clone of the original output. Will perform merges on this image	
	typedef itk::ImageDuplicator< OutputImageType > DuplicatorType; 
	DuplicatorType::Pointer duplicator1 = DuplicatorType::New(); 
	duplicator1->SetInputImage(bImage); 
	duplicator1->Update();
	clonedImage = duplicator1->GetOutput();	

	std::cout<<"No.of valid hypotheses:"<<mergeindices.size()<<std::endl;

	for(int i =0 ; i< mergeindices.size() ; ++i)
	{
		set<int> currRPS;
		set<int>::iterator RPSIterator;
		currRPS =  hypothesis[mergeindices[i]];
		RPSIterator = currRPS.begin();
		unsigned short root = *RPSIterator;

		std::cout<<"Merged"<<"--->"<<root<<"-";

		++RPSIterator;
		for(; RPSIterator != currRPS.end(); ++RPSIterator)
		{
			labelChange(root,*RPSIterator);
			std::cout<<*RPSIterator<<"-";
		}

		std::cout<<""<<std::endl;		
	}	
	return clonedImage;
}



void model_nucleus_seg::RemoveSmallComponents(int minObjSize,const char* x)
{
	typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelFilterType;
	RelabelFilterType::Pointer relabel = RelabelFilterType::New();
	relabel->SetInput( this->clonedImage );
	relabel->SetMinimumObjectSize( minObjSize );
	relabel->InPlaceOn();
	
	try
    {
		relabel->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
		std::cerr << "Relabel: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return;
    }
	std::cout<<"Removing Small Components"<<std::endl;
	writeImage<OutputImageType>(relabel->GetOutput(),x);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// SHATTERS NUCLEI 
///////////////////////////////////////////////////////////////////////////////////////////////////////
int model_nucleus_seg::GetYousefSeg()
{
	RegionType region1;	
	InputImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
	//Allocate Memory for Images
	OutputImageType::IndexType start;
	start.Fill(0);
	region1.SetSize(size);
	region1.SetIndex( start );
	bImage = OutputImageType::New();
	bImage->SetRegions( region1 );
	bImage->Allocate();

	unsigned char *in_Image;
	in_Image = (unsigned char *) malloc (size[0]*size[1]*size[2]);
	//in_Image = (unsigned char *) malloc (size[0]*size[1]*size[2]);

	if( in_Image == NULL )
	{
		std::cout<<"Memory allocation for image failed\n";
		return 0;
	}
	memset(in_Image/*destination*/,0/*value*/,size[0]*size[1]*size[2]*sizeof(unsigned char)/*num bytes to move*/);

	IIteratorType pix_buf(inputImage,inputImage->GetRequestedRegion());
	int ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
		in_Image[ind]=(pix_buf.Get());

	yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->readParametersFromFile("");
	NucleusSeg->setDataImage(in_Image,size[0],size[1],size[2],"");

	// correspond to any existing ids. So, get the max label id in 
	// the segmented image and add this label to all the new fragments

	unsigned short *output_img;
	//segmentation steps
	//1-Binarization
	NucleusSeg->runBinarization();
	NucleusSeg->runSeedDetection();
	NucleusSeg->runClustering();
	NucleusSeg->runAlphaExpansion();
	output_img=NucleusSeg->getSegImage();


	 this->minScale = NucleusSeg->getScaleMin();
	 this->maxScale = NucleusSeg->getScaleMax();

	IteratorType biterator(bImage,bImage->GetRequestedRegion());
	for(unsigned int i=0; i<size[0]*size[1]*size[2]; i++)
	{		
		biterator.Set(output_img[i]);
		++biterator;	
	}
	return 1;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE SHAPE RATIO (NOT BEING USED CURRENTLY IN THE MODEL
///////////////////////////////////////////////////////////////////////////////////////////////////////

//double GetShaperatio(OutputImageType2D::Pointer im)
//{
//	
//	typedef itk::Vector<double,myDimension>  VectorType;
//	typedef itk::ImageMomentsCalculator<OutputImageType> CalculatorType1;
//	
//	CovarianceType inMat;//itkVariablesize matrix
//	inMat.SetSize(2,2);
//	CovarianceType meanMat;
//	meanMat.SetSize(2,2);
//	CovarianceType diffMat;
//	diffMat.SetSize(2,2);
//	
//	vnl_matrix_fixed<double,2,2> inMatvnl ;
//	double cellvol =0;
//	inMat[0][0] = 0;
//	inMat[1][1] = 0;
//	inMat[0][1] = 0;
//	inMat[1][0] = 0 ;										 
//	
//	IteratorType2D reseterator(im,im->GetRequestedRegion());
//	
//	for (reseterator.GoToBegin(); !reseterator.IsAtEnd();++reseterator)
//	{
//		if(reseterator.Get()>0)
//		{
//			OutputImageType2D::IndexType idx = reseterator.GetIndex();
//			inMat[0][0] += ((idx[1])*(idx[1]));
//			inMat[1][1] += ((idx[0])*(idx[0]));
//			inMat[0][1] += ((idx[0])*(idx[1]));
//			inMat[1][0] += ((idx[0])*(idx[1]));
//			cellvol++;				
//		}	
//	} 
//
//		/* Compute the moments */
//		CalculatorType1::Pointer moments = CalculatorType1::New();
//		moments->SetImage( im );
//		moments->Compute();
//		VectorType ccg = moments->GetCenterOfGravity();
//		meanMat[0][0] = ccg[1]*ccg[1]; 
//		meanMat[0][1] = ccg[1]*ccg[0];
//		meanMat[1][0] = ccg[1]*ccg[0];
//		meanMat[1][1] = ccg[0]*ccg[0]; 
//	
//	//DiffMat		
//		for(int i =0;i<2;++i)
//			{									 
//			  for(int j =0;j<2;++j)
//				 {
//					diffMat[i][j] = inMat[i][j]/cellvol - meanMat[i][j]; 
//				 }
//			}
//	
//			inMatvnl = diffMat.GetVnlMatrix();
//			double detval = vnl_determinant(inMatvnl);
//			double shaperatio = cellvol/(4*3.14*sqrt(detval)) ;							 
//											
//	return shaperatio;
//}				


///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE SCORE FOR EVERY ID. BASED ON THE LIBAGF/GENERATIVE MODEL
///////////////////////////////////////////////////////////////////////////////////////////////////////

double model_nucleus_seg::GetScoreforId(std::vector<double> filtered_features)
{
	//Compute region features:	
	double **test;
	//Load the feature vector of the current node after the merge
	double myscore = 0;
	
	for(int c=0; c< NUM_FEAT+NUM_ASSOC_FEAT; ++c)
	{
		filtered_features[c] = (filtered_features[c]-trMean[c])/(trSTD[c]);
	}		

	filtered_features = Convert2EigenSpace(filtered_features);

	test = (double **) malloc(1*sizeof(double *));
	test[0] = (double *)malloc(NUM_COMP * sizeof(double));

	for(int i=0; i<NUM_COMP;++i)
	{
		test[0][i] = filtered_features[i];
	}

	// need to get the data in the format that libagf requires
	train4class = PrepareDataforScores();
	
	
	for(int i=1;i<=NUM_CLASS ;++i)
	{
		std::vector<int> classes2 = classIndex[i-1];
		double** train = train4class[i-1];
		//	Calculate the generative model score ! 
		myscore =  myscore + CalcGenModelScore(train,classes2.size(),NUM_COMP,(long)1,test)*prior[i-1]; 
	}
		
		return myscore;
}


//void model_nucleus_seg::GetIndices()
//{
//	for(int k =0; k < this->NUM_FEAT; ++k)
//	{
//	//Converts upper case to lower case
//	std::transform(features->Info[i].name.begin(), features->Info[i].name.end(),features->Info[i].name.begin(),::tolower);	
//	if(this->featureNames[k] == "convexity")
//		this->convexityIndex = k
//	if(this->featureNames[k] == "volume")	
//		this->volIndex = k;
//	if(this->featureNames[k] == "eccentricity")	
//		this->eccIndex = k;
//	}
//}

std::vector<double> model_nucleus_seg::GetOriginalScores()
{
	std::vector<double> myscore;

	std::vector< std::vector< double > > objectFeatures; 
	objectFeatures.clear();
	objectFeatures.resize(allFeat.size());

	std::cout<< "Computing Original Features" <<std::endl;
	objectFeatures = GetFeaturesOriginal();	
	myscore.resize(objectFeatures.size());		

	int counter = 0;

	double **test;
	test = (double **) malloc(1*sizeof(double *));

	std::cout<< "Finished Computing Original Features" <<std::endl;


	for(int g=0; g < allFeat.size() ; g++)
	{
		std::vector<double> filtered_features2;
		filtered_features2.resize(NUM_FEAT+NUM_ASSOC_FEAT);
		filtered_features2 = objectFeatures.at(g);
		FeatureCalcType::LabelPixelType id = allFeat[g].num; 						
				
		
		double currVol = allFeat[g].ScalarFeatures[ftk::IntrinsicFeatures::VOLUME];
		double eccen = allFeat[g].ScalarFeatures[ftk::IntrinsicFeatures::ECCENTRICITY];
		double sb = allFeat[g].ScalarFeatures[ftk::IntrinsicFeatures::SHARED_BOUNDARY];


		if(!NUM_ASSOC_FEAT)
		{
			for(int c=0; c< NUM_FEAT; ++c)
			{
				filtered_features2[c] = (filtered_features2[c]-trMean[c])/(trSTD[c]);
			}		
		}
		else 
		{
			for(int c = 0; c< NUM_FEAT+NUM_ASSOC_FEAT; ++c)
			{
				filtered_features2[c] = (filtered_features2[c]-trMean[c])/(trSTD[c]);
			}		
		}

		filtered_features2 = Convert2EigenSpace(filtered_features2);
		test[0] = (double *)malloc(NUM_COMP * sizeof(double));

		for(int i=0; i<NUM_COMP;++i)
		{
			test[0][i] = filtered_features2[i];
		}
		
		// Insert ids into boolMerge tree if the shared bdary > 0 and 
		// for eccentricity values > 0.9 or if the fragment is really small
		// eccentricity is considered to account for oversegmented endothelials
		//if(sb>0 && (currVol < MIN_VOL || eccen>0.7))
		if(sb>0 && (currVol < MIN_VOL || eccen>0.9))
		{
			boolMergeTree.push_back(id);
		}

		//Format the data to input it to LibAGF	
		train4class = PrepareDataforScores();
		myscore[g] = 0;
		
		for(int classCount=1;classCount<=NUM_CLASS;classCount++)
		{
			double** train = train4class[classCount-1];
			std::vector<int> classes2 = classIndex[classCount-1];
			//	Calculate the generative model score ! 
			myscore[g]  =  myscore[g]+CalcGenModelScore(train,	classes2.size(),NUM_COMP,(long)1,test)*prior[classCount-1];	
		}
			counter = counter +1;	
	}

	// display the scores for reference on the command prompt.
	DispScores(myscore);

	counter = 0;
	return myscore;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// INSERT ALL THE FRAGMENTS WHICH HAVE BEEN SPLIT INTO BOOLMERGETREE
// BOOLMERGETREE CHECKS IF A MERGE TREE NEEDS TO BE FORMED FOR A PARTICULAR FRAGMENT
///////////////////////////////////////////////////////////////////////////////////////////////////////

void model_nucleus_seg::UpdateBoolMerge()
{
	for(unsigned short i = maxL+1; i<=allFeat[allFeat.size()-1].num ; ++i)
	{
		if(boolMergeTree.end() == find(boolMergeTree.begin(), boolMergeTree.end(), i))
		{
			boolMergeTree.push_back(i);
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// DISPLAY THE SCORE ON THE COMMAND PROMPT FOR REFERENCE 
///////////////////////////////////////////////////////////////////////////////////////////////////////
void model_nucleus_seg::DispScores(std::vector<double> scores)
{
	originalScores.resize(allFeat.size());
	for (unsigned short id = 0; id<allFeat.size();id++)
	{	
		unsigned short numval = allFeat[id].num; 
		originalScores[id] =  scores[id];
		//std::cout<<numval<<"-"<<originalScores[id]<<std::endl;
	}	
}


//****************************************************************************
// RunSVM - run support vector machine in one class mode to find outliers.
//****************************************************************************
std::vector<unsigned short> model_nucleus_seg::runSVM(vnl_matrix<double> feats,model_nucleus_seg::FeatureCalcType* filter)
{
	// Get all the samples of individual classes into a matrix and
	// store the matrices in a vector featsbyClass
	
	std::vector< vnl_matrix<double> > FeatsbyClass;
	//FeatsbyClass.resize(NUM_CLASS);
	

	for(unsigned int classnumb=0; classnumb < NUM_CLASS; ++classnumb)
	{
		vnl_matrix<double> temp_Matrix(10000,NUM_FEAT+NUM_ASSOC_FEAT);
		int rowcounter = 0;
		for(unsigned int r=0; r<feats.rows(); ++r)
			{
				vnl_vector<double> temp_row = feats.get_row(r);
				//(temp_row.get(temp_row.size()-2)) gives the class
				if(temp_row.get(temp_row.size()-2) == classnumb+1) // since classIndex starts from 1 and not 0
				{
					for(unsigned int featcount=0; featcount<NUM_FEAT+NUM_ASSOC_FEAT; ++featcount)
						temp_Matrix(rowcounter,featcount) = temp_row(featcount);  
					rowcounter++;
				}
			}
		FeatsbyClass.push_back(temp_Matrix.get_n_rows(0,rowcounter-1));
	}
	
	
	double nu = SPLIT_SENSITIVITY;
	double vol;

	int myrows = feats.rows();
	int mycols = feats.columns()-2;//Includes both associative and intrinsic features 

	//Setup the scaling values
	std::vector<double> f_max;                        //The maximum value of each feature
	f_max.assign(mycols,-(DBL_MAX));
	std::vector<double> f_min;                        //The minimum value of each feature
	f_min.assign(mycols,(DBL_MAX));

	std::vector< std::vector< double > > objfeatures;    //Will contain the normalized features of the objects in the image


	//Read Training File:
	double val;

	//Normalize/Scale the Features Data:
	double upper = 1;
	double lower = -1;
	double desired_range = upper-lower;

	
	//Get the volumes for each class and the max/min values of features
	// for normalization
	for(unsigned int r=0; r< myrows; ++r)
	{	
		for(int c=0; c < mycols ; ++c) // Edit
		{	
			val = feats(r,c);
			f_max.at(c) = val > f_max.at(c) ? val : f_max.at(c);
			f_min.at(c) = val < f_min.at(c) ? val : f_min.at(c);
		}
	}

	
	//Normalize the rows 
	for(unsigned int classnumb=0; classnumb < NUM_CLASS; ++classnumb)
	{
		vnl_matrix<double> temp_Matrix = FeatsbyClass[classnumb];
		for(unsigned int r=0; r<temp_Matrix.rows(); ++r)
		{
			for(int c=0; c < mycols ; ++c) // Edit
			{	
				double oldval = temp_Matrix(r,c);
				temp_Matrix(r,c) = lower + desired_range * (oldval - f_min.at(c)) / (f_max.at(c) - f_min.at(c));
			}
			FeatsbyClass[classnumb] = temp_Matrix;
		}
	}

	
	//Predict:

	std::vector<unsigned short> outliers;
	std::vector< FeatureCalcType::LabelPixelType > labels = filter->GetLabels();
	objfeatures.resize(labels.size()-1);
	objfeaturesNN.resize(labels.size()-1);


	//First get the Intrinsic features
	for (unsigned long r = 1; r<labels.size();++r)
	{
		FeatureCalcType::LabelPixelType id = labels.at(r); 
		ftk::IntrinsicFeatures* features = filter->GetFeatures(id);

		for(int k =0; k < this->NUM_FEAT; ++k)
		{
			for(int i=0; i<32;++i)	
			{
				//Converts upper case to lower case
				std::transform(features->Info[i].name.begin(), features->Info[i].name.end(),features->Info[i].name.begin(),::tolower);	
				if(this->featureNames[k] ==  features->Info[i].name)
				{
					val = features->ScalarFeatures[i];
					objfeatures.at(r-1).push_back(val);
					objfeaturesNN.at(r-1).push_back(val);
					//std::cout<<val<<" ";
				}
			}
		}
		//std::cout<<""<<std::endl;
	}


	//Then get the Associative Features here if present
	if(this->NUM_ASSOC_FEAT>0)
	{
		for (unsigned long r = 1; r<labels.size();++r)
		{
			std::vector<double> assoc = this->AssocFeat[r-1];
			for(int i=0; i< NUM_ASSOC_FEAT ;++i)	
				{
					val = assoc[i];
					objfeatures.at(r-1).push_back(val);
					objfeaturesNN.at(r-1).push_back(val);
				}
		}
	}

	//getch();

	for(unsigned int r=0; r< objfeatures.size(); ++r)
	{
		for(int c=0; c < mycols ; ++c) // Edit
		{
			double oldval = objfeatures.at(r).at(c);
			objfeatures.at(r).at(c) = lower + desired_range * (oldval - f_min.at(c)) / (f_max.at(c) - f_min.at(c));
		}
	}
	

	

		//Create the libSVM problem:
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

	// Will store 1 in the column if the id ( nucleus ) is an outlier corresponding to this model
	// rows of the matrix correspond to different ids
	vnl_matrix<int> id_outlier(objfeatures.size(),NUM_CLASS);
		
	for(unsigned int classnumb=0; classnumb < NUM_CLASS; ++classnumb)
	{	

		struct svm_problem prob;
		vnl_matrix<double> Feats_classnumb = FeatsbyClass[classnumb];
		prob.l = Feats_classnumb.rows();                    //Number of objects
		prob.y = Malloc(double,prob.l);                    //Array Containing target values (unknowns)
		prob.x = Malloc(struct svm_node *,prob.l);        //Array of Pointers to Nodes

		for(int r=0; r<prob.l; ++r)
		{
			prob.y[r] = 1;                                //This is the label (target) and it is unknown

			struct svm_node *x_space = Malloc(struct svm_node, mycols+1);    //Individual node -> Volume + Intensity + Class = 3 cols

			for(int c=0; c< mycols ; ++c)
			{
				x_space[c].index = c+1;
				x_space[c].value = Feats_classnumb(r,c);
			}
			x_space[mycols].index = -1;
			prob.x[r] = &x_space[0];    //Point to this new set of nodes.
		}

		//Set the Parameters
		struct svm_parameter param;
		param.svm_type = ONE_CLASS;
		param.kernel_type = RBF;
		param.degree = 3;
		param.gamma = 1.0/double(prob.l);    // 1/k
		//param.gamma = 1;
		param.coef0 = 0;
		//param.nu = 0.1;
		param.nu = nu;
		param.cache_size = 100;
		param.C = 1;
		param.eps = .001;
		param.p = 0.1;
		param.shrinking = 1;
		param.probability = 0;
		param.nr_weight = 0;
		param.weight_label = NULL;
		param.weight = NULL;

		//Now train
		struct svm_model *m_svm_model;
		m_svm_model = svm_train(&prob,&param);

		svm_destroy_param(&param);
		free(prob.y);
		free(prob.x);


		prob.l = (int)objfeatures.size();                    //Number of objects
		prob.y = Malloc(double,prob.l);                    //Array Containing target values (unknowns)
		prob.x = Malloc(struct svm_node *,prob.l);        //Array of Pointers to Nodes


		for(int r=0; r<prob.l; r++)
		{	
			int outlierFlag = 1;
			struct svm_node *x = Malloc(struct svm_node,mycols);            //Individual node

			for(int c=0; c < mycols; ++c)
			{
				x[c].index = c+1;
				x[c].value = objfeatures.at(r).at(c);
			}
			x[mycols].index = -1;

			double currVol = x[0].value; // current volume
			double v = svm_predict(m_svm_model,x);

			if(v ==-1 && currVol>0)
				id_outlier(r,classnumb) = 1;	
		}

	}

	
	for(unsigned int r=0; r< id_outlier.rows(); ++r)
	{	
		int outlierflag = 1;
		for(int c=0; c < id_outlier.cols() ; ++c) // Edit
		{	
			int x123 = id_outlier(r,c);
			if(x123!=1)
			{
				outlierflag = 0;
				break;
			}
		}
		if(outlierflag==1)
		{
			std::cout<<labels.at(r+1)<< " is an outlier"<<std::endl;
			outliers.push_back(labels.at(r+1));
		}

	}
	//svm_destroy_model(m_svm_model);
	return outliers;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// SPLIT IMAGE BASED ON THE OUTLIERS 
///////////////////////////////////////////////////////////////////////////////////////////////////////
model_nucleus_seg::OutputImageType::Pointer model_nucleus_seg::SplitImage(std::vector<unsigned short> outliers,InputImageType::Pointer rimage, OutputImageType::Pointer limage)
{
	//Filter to get the region associated with the id
	LabelStatsFilterType::Pointer lsfilter = LabelStatsFilterType::New();
	lsfilter->SetInput(rimage);
	lsfilter->SetLabelInput(limage); 
	lsfilter->Update();	
	
	unsigned long maxLabel=0;


	std::vector<unsigned short>::iterator otIterator;
	for(otIterator = outliers.begin(); otIterator != outliers.end(); otIterator++)
	{
		int id  = *otIterator;

		// region will contain the region corresponding to the 
		// label which is an outlier. Since region is rectangular,
		// it may contain a part of the region of other labels as well
		RegionType region1;	
		region1 = lsfilter->GetRegion(id);	

		InputImageType::SizeType size1 = region1.GetSize();
		InputImageType::IndexType start1 = region1.GetIndex();

		//roiImage is the grayscale image which has been under segmented
		InputImageType::Pointer roiImage  = InputImageType::New();
		roiImage  = returnROI<InputImageType>(rimage,region1);

		OutputImageType::Pointer lroiImage  = OutputImageType::New();
		lroiImage  = returnROI<OutputImageType>(limage,region1);


		//Reset the background pixels to zero
		IteratorType reseterator(lroiImage,lroiImage->GetRequestedRegion());
		IIteratorType fiterator(roiImage,roiImage->GetRequestedRegion());


		for (reseterator.GoToBegin(),fiterator.GoToBegin(); !reseterator.IsAtEnd();++reseterator,++fiterator )
		{
			if(reseterator.Get()!=id)
				fiterator.Set(0);
			else
				fiterator.Set(255);
		} 

		OutputImageType::Pointer pasteImage  = OutputImageType::New();
		
		// This function shatters the neculei into small fragments
		// which can later be merged.
		pasteImage = Shatter_Nuclei(roiImage,&maxLabel);
		IteratorType piterator(pasteImage,pasteImage->GetRequestedRegion());

		//Reset the background pixels in pasteImage to "background" pixels in the original segmented image
		//This is required because the region is rectangular and the regions are arbitrarily shaped
		for (reseterator.GoToBegin(),piterator.GoToBegin(); !reseterator.IsAtEnd();++reseterator,++piterator )
		{
			if(piterator.Get()==0 && reseterator.Get()!=id)
			{
				piterator.Set(reseterator.Get());
			}
		}

		PasteFilterType::Pointer pfilter = PasteFilterType::New();
		pfilter->SetDestinationImage(limage);
		pfilter->SetSourceImage( pasteImage );

		PasteFilterType::InputImageIndexType destIndex;
		destIndex = region1.GetIndex();
		pfilter->SetDestinationIndex( destIndex );
		pfilter->SetSourceRegion(pasteImage->GetLargestPossibleRegion());
		pfilter->Update();
		limage = pfilter->GetOutput();

	}		

	return limage;
}




///////////////////////////////////////////////////////////////////////////////////////////////////////
// SHATTERS NUCLEI 
///////////////////////////////////////////////////////////////////////////////////////////////////////
model_nucleus_seg::OutputImageType::Pointer model_nucleus_seg::Shatter_Nuclei(InputImageType::Pointer roiImage,unsigned long * maxLabel)
{
	RegionType region1;	
	InputImageType::SizeType size = roiImage->GetLargestPossibleRegion().GetSize();
	//Allocate Memory for Images
	OutputImageType::IndexType start;
	start.Fill(0);
	region1.SetSize(size);
	region1.SetIndex( start );
	bImagetemp = OutputImageType::New();
	bImagetemp->SetRegions( region1 );
	bImagetemp->Allocate();


	unsigned char *in_Image;
	in_Image = (unsigned char *) malloc (size[0]*size[1]*size[2]);
	//in_Image = (unsigned char *) malloc (size[0]*size[1]*size[2]);

	if( in_Image == NULL )
	{
		std::cout<<"Memory allocation for image failed\n";
		return 0;
	}
	memset(in_Image/*destination*/,0/*value*/,size[0]*size[1]*size[2]*sizeof(unsigned char)/*num bytes to move*/);

	IIteratorType pix_buf(roiImage,roiImage->GetRequestedRegion());
	int ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
		in_Image[ind]=(pix_buf.Get());
	
	
	if(*maxLabel==0)
	{
	//Define a filter to access all the ids
	LabelStatsFilterType::Pointer lsfilter = LabelStatsFilterType::New();
	lsfilter->SetInput(inputImage);
	lsfilter->SetLabelInput(bImage); 
	lsfilter->Update();	
	*maxLabel = lsfilter->GetNumberOfLabels()-1;
	}

	yousef_nucleus_seg *NucleusSeg = new yousef_nucleus_seg();
	NucleusSeg->readParametersFromFile("");
	NucleusSeg->setDataImage(in_Image,size[0],size[1],size[2],"");

	// We need the new fragments to have unique Ids which do not
	// correspond to any existing ids. So, get the max label id in 
	// the segmented image and add this label to all the new fragments

	unsigned short *output_img;
	//segmentation steps
	//1-Binarization
	NucleusSeg->runBinarization();
	NucleusSeg->runSeedDetection(this->minScale,this->maxScale);
	NucleusSeg->runClustering();
	//NucleusSeg->runAlphaExpansion();
	output_img=NucleusSeg->getClustImage();

	IteratorType biterator(bImagetemp,bImagetemp->GetRequestedRegion());
	
	unsigned long maxlocalID = 0;

	for(unsigned int i=0; i<size[0]*size[1]*size[2]; i++)
	{	
		if(output_img[i]==0)
			biterator.Set(output_img[i]);
		else
		{
			biterator.Set(output_img[i]+*maxLabel);
			if(output_img[i]+*maxLabel > maxlocalID)
				maxlocalID = output_img[i]+*maxLabel;
		}
		++biterator;	
	}
	

	*maxLabel = maxlocalID;
	return bImagetemp;  		
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// CALCULATE THE GENERATIVE MODEL SCORE USING LIBAGF
///////////////////////////////////////////////////////////////////////////////////////////////////////
double model_nucleus_seg::CalcGenModelScore(double** train, long ntrain,long nvar,long ntest,double **test) 
{
	FILE *fs;
	FILE *diagfs;			//print diagnostics to this file stream


	double *result;	//results of pdf estimation

	agf_command_opts opt_args;

	double *std, *ave;

	//diagnostics:
	agf_diag_param diag_param;

	long min_nd, max_nd, total_nd;
	double min_f, max_f, total_f;
	double min_W, max_W, total_W;

	//set defaults and parse command line options:
	opt_args.k=-1;
	opt_args.Wc=WC_DEFAULT;
	//if the initial filter variance is not set, set it to the total
	//variance of the data:	

	opt_args.normflag = 1;
	opt_args.var_0 = 1;

	std=NULL;
	ave=NULL;

	//fprintf(diagfs, "Objective total weight: Wc=%f\n", opt_args.Wc);
	//begin the classification scheme:
	result=new double[ntest];
	//initialize diagnostic values:
	min_nd=100;
	max_nd=0;
	total_nd=0;

	min_f=1;
	max_f=0;
	total_f=0;

	//	min_W=100*opt_args.Wc;
	max_W=0;
	total_W=0;

	for (long i=0; i<ntest; i++) 
	{
		result[i]=agf_calc_pdf(train, nvar, ntrain, test[i],opt_args.var_0, opt_args.Wc, &diag_param);
		if(result[i]!=result[i]) //if(_isnan(result[i]))
		{
			result[i] = 0;
		}
	}	

	return result[0];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//PREPARE THE DATA FOR SCORES 
///////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double **> model_nucleus_seg::PrepareDataforScores()
{
	std::vector<double **> train4class2;
	train4class2.resize(NUM_CLASS);

	classIndex.clear();
	classIndex.resize(NUM_CLASS);

	//Read the size of training file:
	int myrows = Num_Lines2(myfilenamePCA, 0) ;   // number of rows 			
	double **train;
	std::vector<int> classes2;	

	for(int i=0; i<myrows;++i)
	{
		classes2 = classIndex[myKnownClass[i]-1];
		classes2.push_back(i);
		classIndex[myKnownClass[i]-1] = classes2;
	}

	for(int i=0; i<NUM_CLASS;++i)
	{
		std::vector<int> classes2 = classIndex.at(i);		
		train = (double **) malloc(classes2.size()*sizeof(double *));
		for(int j = 0 ; j < classes2.size() ; ++j)
		{
			train[j] = (double *)malloc(NUM_COMP * sizeof(double));  

			for(int k = 0 ; k < NUM_COMP ; ++k )
			{
				train[j][k] = PCAdata[classes2.at(j)][k];	
			}
		}
		train4class2[i] = train;
	}
	return train4class2;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//SELECT THE BEST HYPOTHESIS AND CALL PERFORM MERGES
////////////////////////////////////////////////////////////////////////////////////////////////

void model_nucleus_seg::SelectHypothesis()
{
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_prob_name(lp, "Sel_Hypo");
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, allFeat.size());	

	long p1 = hypoMatrix.size()*allFeat.size()+1;
	
	int *ia = (int *)malloc(p1*sizeof(int));
	int *ja = (int *)malloc(p1*sizeof(int));
	double *ar = (double *)malloc(p1*sizeof(double));

	

	//Setting the upper and lower bounds of the variables

	for(unsigned int i=0 ; i < allFeat.size(); ++i)
	{
		glp_set_row_bnds(lp, i+1, GLP_UP, 0.0, 1.0);
	}

	glp_add_cols(lp, hypoMatrix.size());
	
	//Setting the co-eff of the objective function
	for(unsigned int i=0 ; i < hypoMatrix.size(); ++i)
	{
		glp_set_col_kind(lp,i+1,GLP_BV); 
		glp_set_col_bnds(lp, i+1, GLP_DB, 0.0, 1.0);

		std::vector <double> sumrow = hypoMatrix[i];
		double sum = 0.0;


		for(unsigned int j=0; j < sumrow.size(); ++j)
		{
			sum = sum + sumrow[j];
		}

		double range = i +1 ;
		//std::cout<<sum<<std::endl;
		glp_set_obj_coef(lp,range,sum);

	}
																				
	
	unsigned int ctr = 1;

	//Coefficients of the constrains ... B' matrix
	for(unsigned int i = 0 ; i< hypoMatrix.size() ; ++i)
	{
		for(unsigned int j = 0 ; j < allFeat.size() ; ++j)	
		{
			ia[ctr] = j+1;
			ja[ctr] = i+1;

			std::vector <double> transposedrow = hypoBMatrix[i];

			ar[ctr] =  transposedrow[j]  ; //B'
			++ctr;
		}
	}

	glp_load_matrix(lp, ctr-1, ia, ja, ar);
	glp_simplex(lp,NULL);
	glp_intopt(lp,NULL); 

	for(unsigned int i =0; i<hypoMatrix.size(); ++i)
	{
		if(glp_mip_col_val(lp,i+1) == 1)
		{
			mergeindices.push_back(i);
			//std::vector<double> row = hypoBMatrix[i];

			//for (unsigned int j =0; j< row.size() ; ++j)
			//{
			//	std::cout<< row[j] << "  " ;
			//}
			//std::cout<<i<<std::endl;
		}
	}

}

////Loads the XML file 
///////////////////////////////////////////////////////////////////////////////////////////////////////
//LOAD THE PROJECT DEFINITION FILE
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool model_nucleus_seg::LoadAssoc(std::string filename)
{
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return false;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "ProjectDefinition" ) != 0 )
		return false;

	std::string name = rootElement->Attribute("name");

	TiXmlElement * parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "Inputs" ) == 0 )
		{
			//inputs = ftk::ReadChannels(parentElement);
		}
		else if( strcmp( parent, "AssociationRules" ) == 0 )
		{
			associationRules = ReadAssociationRules(parentElement);
		}

		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();
	return true;
}



std::vector<ftk::AssociationRule> model_nucleus_seg::ReadAssociationRules(TiXmlElement * inputElement)
{
	std::vector<ftk::AssociationRule> returnVector;

	TiXmlElement * parameterElement = inputElement->FirstChildElement();
	while (parameterElement)
	{
		const char * parameter = parameterElement ->Value();
		ftk::AssociationRule assocRule("");
		if ( strcmp(parameter,"AssociationRule") == 0 )
		{
			assocRule.SetRuleName(parameterElement->Attribute("Name"));
			assocRule.SetSegmentationFileNmae(parameterElement->Attribute("SegmentationSource"));
			assocRule.SetTargetFileNmae(parameterElement->Attribute("Target_Image"));
			assocRule.SetOutDistance(atoi(parameterElement->Attribute("Outside_Distance")));
			assocRule.SetInDistance(atoi(parameterElement->Attribute("Inside_Distance")));

			if(strcmp(parameterElement->Attribute("Use_Whole_Object"),"True")==0)
				assocRule.SetUseWholeObject(true);
			else
				assocRule.SetUseWholeObject(false);

			if(strcmp(parameterElement->Attribute("Use_Background_Subtraction"),"True")==0)
				assocRule.SetUseBackgroundSubtraction(true);
			else
				assocRule.SetUseBackgroundSubtraction(false);

			if(strcmp(parameterElement->Attribute("Use_MultiLevel_Thresholding"),"True")==0){
				assocRule.SetUseMultiLevelThresholding(true);
				assocRule.SetNumberOfThresholds(atoi(parameterElement->Attribute("Number_Of_Thresholds")));
				assocRule.SetNumberIncludedInForeground(atoi(parameterElement->Attribute("Number_Included_In_Foreground")));
			}
			else{
				assocRule.SetUseMultiLevelThresholding(false);
				assocRule.SetNumberOfThresholds(1);
				assocRule.SetNumberIncludedInForeground(1);
			}

			if(strcmp(parameterElement->Attribute("Association_Type"),"MIN")==0)
				assocRule.SetAssocType(ftk::ASSOC_MIN);
			else if(strcmp(parameterElement->Attribute("Association_Type"),"MAX")==0)
				assocRule.SetAssocType(ftk::ASSOC_MAX);
			else if(strcmp(parameterElement->Attribute("Association_Type"),"TOTAL")==0)
				assocRule.SetAssocType(ftk::ASSOC_TOTAL);
			else if(strcmp(parameterElement->Attribute("Association_Type"),"SURROUNDEDNESS")==0)
				assocRule.SetAssocType(ftk::ASSOC_SURROUNDEDNESS);
			else
				assocRule.SetAssocType(ftk::ASSOC_AVERAGE);
		}
		returnVector.push_back(assocRule);
		parameterElement = parameterElement->NextSiblingElement();
	}
	return returnVector;
}



////Loads the Parameters file 
///////////////////////////////////////////////////////////////////////////////////////////////////////
//LOAD THE SEGMENTATION PARAMETERS FILE
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool model_nucleus_seg::LoadSegParams(std::string filename)
{
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return false;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "SegmentationDefinition" ) != 0 )
		return false;

	std::string name = rootElement->Attribute("name");//default

	TiXmlElement * parentElement = rootElement->FirstChildElement();
	const char * parent = parentElement->Value();
		if ( strcmp(parent,"SegParameters") == 0 )
		{
			TiXmlElement * parameterElement = parentElement->FirstChildElement();
			const char * parameter = parameterElement ->Value();
			this->NUM_CLASS = atoi(parameterElement->Attribute("NUM_CLASS"));
			this->NUM_COMP = atoi(parameterElement->Attribute("NUM_COMP"));
			this->NUM_FEAT = atoi(parameterElement->Attribute("NUM_FEAT"));			
			this->NUM_ASSOC_FEAT = atoi(parameterElement->Attribute("NUM_ASSOC_FEAT"));
			this->SPLIT = atoi(parameterElement->Attribute("SPLIT"));
			this->SPLIT_SENSITIVITY = atoi(parameterElement->Attribute("SPLIT_SENSTIVITY"));

			if(parameterElement->Attribute("MAX_VOL"))
				this->MAX_VOL = static_cast<double>(atoi(parameterElement->Attribute("MAX_VOL")));
			if(parameterElement->Attribute("MIN_VOL"))
				this->MIN_VOL = static_cast<double>(atoi(parameterElement->Attribute("MIN_VOL")));
			if(parameterElement->Attribute("MAX_DEPTH"))
				this->MAX_DEPTH =  atoi(parameterElement->Attribute("MAX_DEPTH"));
			if(parameterElement->Attribute("WC_DEFAULT"))
				this->WC_DEFAULT = static_cast<double>(atoi(parameterElement->Attribute("WC_DEFAULT")));
			

			this->prior.resize(this->NUM_CLASS);
			//this->volLimit.resize(this->NUM_CLASS);		
			//Assign Prior Probabilities
			for(unsigned long i=1; i<=this->NUM_CLASS;++i)
			{
				string Prior = "PRIOR"+convert2string(i);
				this->prior[i-1] = strtod((parameterElement->Attribute(Prior))->c_str(),NULL);			
			}
		}
		else 
		{
			printf("Wrong Format for Parameter File\n");
		}

	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//COMPUTE THE ASSOCIATIONS (CHANGES MADE TO THE ORIGINAL )
///////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector < std::vector<double> > model_nucleus_seg::ComputeAssociations(void)
{
	std::vector < double > vals;
	std::vector<std::vector < double > > afeats;
	afeats.resize(allFeat.size());

	int count = 0;
	for(std::vector<ftk::AssociationRule>::iterator ascit=associationRules.begin(); ascit!=associationRules.end(); ++ascit )
	{
		int seg_channel_number=-1;
		unsigned short inp_channel_number=-1;
		if( strcmp( ascit->GetSegmentationFileName().c_str(), "NUCLEAR" ) == 0 )
		{
			seg_channel_number = 0;
		} 
		else
		{
			std::cout<<"Please check region type for associative feature computation\n";
		}

		for( int j=0; j<(int)myFtkImage->GetImageInfo()->channelNames.size(); ++j )
			if( strcmp( myFtkImage->GetImageInfo()->channelNames.at(j).c_str(), ascit->GetTargetFileNmae().c_str() ) == 0 ){
				inp_channel_number=j;
				break;
			}
			if( inp_channel_number == -1 ){
				std::cout<<"Unable to access grayscale image while computing associative feature: "<<ascit->GetRuleName()<<std::endl;
			}
			OutputImageType::Pointer inp_im  = OutputImageType::New();	
			inp_im = myFtkImage->GetItkPtr<unsigned short>(0,inp_channel_number,ftk::Image::DEEP_COPY);
			
			assoc = new ftk::NuclearAssociationRules("",0,bImage, inp_im);
			assoc->AddAssociation( (*ascit).GetRuleName(), "", (*ascit).GetOutDistance(), (*ascit).GetInDistance(),	(*ascit).IsUseWholeObject(), (*ascit).IsUseBackgroundSubtraction(), (*ascit).IsUseMultiLevelThresholding(),(*ascit).GetNumberOfThresholds(), (*ascit).GetNumberIncludedInForeground(), (*ascit).GetAssocType(), (*ascit).get_path() );		

			//I have modified ComputeoneAssoc method to include multiple labels in getmergedfeaturestest.
			// For this i am changing int lbl to std::vector<int> idvec;

			std::vector<int> idvec;
			idvec.resize(1);

		
	//Need to calculate the threshold only once for an image

	if(imgindex.end() == find(imgindex.begin(), imgindex.end(), inp_channel_number))
	{
		if( associationRules[count].IsUseBackgroundSubtraction() )
			{
				if( associationRules[count].IsUseMultiLevelThresholding() )
					if( associationRules[count].GetNumberOfThresholds()>=associationRules[count].GetNumberIncludedInForeground() )
					{
						unsigned short threshval = returnthresh( inp_im, associationRules[count].GetNumberOfThresholds(), associationRules[count].GetNumberIncludedInForeground() );
						mymap[inp_channel_number] = threshval;
					}
					else
					{
						unsigned short threshval = returnthresh( inp_im, associationRules[count].GetNumberOfThresholds(), associationRules[count].GetNumberOfThresholds() );
						mymap[inp_channel_number] = threshval;
					}
				else
				{
					unsigned short threshval =returnthresh( inp_im, 1, 1 );
					mymap[inp_channel_number] = threshval;					
				}
			imgindex.push_back(inp_channel_number);		
			}
	}	
			
		// threshold value for the current channel
		vector<unsigned short>::iterator posn1 = find(imgindex.begin(), imgindex.end(), inp_channel_number);
		thresh	 = mymap[inp_channel_number];
				
			for(int j=0; j<allFeat.size(); ++j)
			{
				int lbl1 = 	allFeat[j].num;
				if(lbl1 == 0) continue;			
				idvec[0] = lbl1;
				vector<unsigned short>::iterator posn1 = find(labelIndex.begin(), labelIndex.end(), lbl1);

				float* bbox = allFeat[posn1 - labelIndex.begin()].BoundingBox;
				x_Size=inp_im->GetLargestPossibleRegion().GetSize()[0];
				y_Size=inp_im->GetLargestPossibleRegion().GetSize()[1];
				z_Size=inp_im->GetLargestPossibleRegion().GetSize()[2];
				
				vals = ComputeOneAssocMeasurement(inp_im,0,idvec,bbox);
				for(int l1 =0; l1<vals.size();l1++)
					afeats[j].push_back(vals[l1]);
			}
			
			count++;
	}

	std::cout << "Done Associations\n";
	return afeats;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
///COMPUTE ASSOCIATIVE MEASUREMENTS FOR ONE CELL
///////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> model_nucleus_seg::ComputeOneAssocMeasurement(itk::SmartPointer<OutputImageType> trgIm, int ruleID, std::vector<int>objID,float* b1)
{	
	//
	//fisrt, get the bounding box around the object
	//The bounding box is defined by the area around the object such that the object+outerdistance are included	
	int imBounds[6] = {0,x_Size,0,y_Size,0,z_Size};	

	//Get the Bounding Box here !!! 

	//< INSERT CODE > 
	int bbox[6] = {b1[0],b1[1],b1[2],b1[3],b1[4],b1[5]};

	std::vector< int > retBbox(0);
	int dist, val;
	dist = assoc->GetAssociationRules().at(ruleID).GetOutDistance();
	for(int dim=0; dim < myDimension*2; ++dim)
	{  
		if(dim%2 == 0) //even
		{			
			val = int(bbox[dim])-dist-1;
			if(val<imBounds[dim])
				val=imBounds[dim];
		}
		else //odd
		{
			val = int(bbox[dim])+dist+1;
			if(val>=imBounds[dim])
				val=imBounds[dim]-1;
		}
		retBbox.push_back( val );
	}

	//the bounding box defines the region of interest we need (from both segmentation and target images)
	//so, use it to get sub images
	DistImageType::Pointer subSegImg = DistImageType::New();

	OutputImageType::IndexType start;
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    

	OutputImageType::SizeType  size;
	size[0]  = retBbox[1]-retBbox[0]+1;  // size along X
	size[1]  = retBbox[3]-retBbox[2]+1;  // size along Y
	if(myDimension == 3)
		size[2]  = retBbox[5]-retBbox[4]+1;  // size along Z
	else
		size[2] = 1;

	DistImageType::RegionType region1;
	region1.SetSize( size );
	region1.SetIndex( start );
	subSegImg->SetRegions( region1 );	
	subSegImg->Allocate();
	subSegImg->FillBuffer(0);
	subSegImg->Update();	

	DistImageType::IndexType start2;
	start2[0] =   retBbox[0];  // first index on X
	start2[1] =   retBbox[2];  // first index on Y    
	if(myDimension == 3)
		start2[2] =   retBbox[4];  // first index on Z   
	else
		start2[2] = 0;
	OutputImageType::RegionType region2;
	region2.SetSize( size );
	region2.SetIndex( start2 );
	bImage->SetRequestedRegion(region2);
	trgIm->SetRequestedRegion(region2);
	IteratorType iterator1(bImage, bImage->GetRequestedRegion());
	DIteratorType iterator2(subSegImg, subSegImg->GetRequestedRegion());

	//in the sub-segmentation image, we need to mask out any pixel from another object	
	int counter = 0;
	while ( ! iterator1.IsAtEnd())
	{		
		int V = iterator1.Get();
		if(objID.end() != find(objID.begin(), objID.end(), V))
		{
			iterator2.Set(255.0);
		}
		else
			iterator2.Set(0.0);
		++iterator1;
		++iterator2;		
	}	

	//Compute the distance transform in the sub-segmentation image region	
	DTFilter::Pointer dt_obj= DTFilter::New() ;
	dt_obj->SetInput(subSegImg) ;	
	try{
		dt_obj->Update() ;
	}
	catch( itk::ExceptionObject & err ){
		std::cout << "Error in Distance Transform: " << err << std::endl; 
		//return 0;	
	}


	//now, mask out all the pixels (in the sub-seg image) that are not in the region of interest as defined by the association rule and get the intensities of the needed pixels from the target image. The intensities are saved into an std vector
	DIteratorType iterator3(dt_obj->GetOutput(), dt_obj->GetOutput()->GetRequestedRegion());
	IteratorType iterator4(trgIm, trgIm->GetRequestedRegion());
	std::vector<int> trgInt;
	int counter_in = 0;
	int counter_at = 0;
	int counter_ot = 0;
	while ( ! iterator3.IsAtEnd())
	{
		int V = (int)iterator3.Get();
		if(V<0)
			counter_in++;
		if(V==0)
			counter_at++;
		if(V>0)
			counter_ot++;
		//if it is outside with distance less than outDistance away
		if(V>0 && V<=assoc->GetAssociationRules().at(ruleID).GetOutDistance())
			trgInt.push_back(iterator4.Get());
		//if it is inside and the whole cell is used
		else if(V<=0 && assoc->GetAssociationRules().at(ruleID).IsUseWholeObject())
			trgInt.push_back(iterator4.Get());
		//if it is inside with distance less than in Distance
		else if(V<=0 && abs(V)<=assoc->GetAssociationRules().at(ruleID).GetInDistance())
			trgInt.push_back(iterator4.Get());

		++iterator3;
		++iterator4;
	}
//	if(!trgInt.size())
	std::vector<double> assocs;
   
	assocs.push_back(ComputeTotal(trgInt));
    assocs.push_back(ComputeAverage(trgInt));
    //assocs.push_back(ComputeSurroundedness( trgIm, this->bImage, 8, thresh, dist,objID));

	//we will never go here, just silencing a compiler warning
	return assocs;
}


double model_nucleus_seg::ComputeTotal(std::vector<int> LST)
{
	float tl = 0;
	for(unsigned int i=0; i<LST.size(); i++)
	{
		if( LST[i] >= thresh )
			tl += LST[i];
	}
	return tl;
}


//compute average
double model_nucleus_seg::ComputeAverage(std::vector<int> LST)
{
	float av = 0;
	float lst_sz = (float)LST.size();
	for(unsigned int i=0; i<LST.size(); i++)
	{
		if( LST[i] >= thresh )
			av += LST[i];
		else
			--lst_sz;
	}
	if( lst_sz )
		av/=lst_sz;
	return av;
}



bool model_nucleus_seg::sorthelp (double i,double j) 
{ 
	return (i<j); 
}


bool model_nucleus_seg::uniquehelp (unsigned short i, unsigned short j) {
	return (i==j);
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
//NORMALIZE A SET OF VECTORS BY THE STD DEVIATIONS OF EACH OF THE VARIABLES
///////////////////////////////////////////////////////////////////////////////////////////////////////
void model_nucleus_seg::norm_vec_std(double **mat1, long m, long n1, double **mat2, long n2) 
{
	double *ave;
	double *std;
	double anom;		//intermediate output

	//allocate space:
	ave=new double[m];
	std=new double[m];

	for (long i=0; i<m; i++) {
		ave[i]=0;
		for (long j=0; j<n1; j++) ave[i]+=mat1[j][i];
		ave[i]/=n1;
		std[i]=0;
		for (long j=0; j<n1; j++) {
			anom=mat1[j][i]-ave[i];
			std[i]+=anom*anom;
			mat1[j][i]=anom;
		}
		std[i]=sqrt(std[i]/(n1-1));

		//normalize the first matrix:
		for (long j=0; j<n1; j++) mat1[j][i]=mat1[j][i]/std[i];

		//normalize the second matrix:
		for (long j=0; j<n2; j++) mat2[j][i]=(mat2[j][i]-ave[i])/std[i];
	}
	//clean up:
	delete [] ave;
	delete [] std;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// GET THE NUMBER OF ROWS AND COLUMNS IN THE TRAINING FILE
///////////////////////////////////////////////////////////////////////////////////////////////////////

int model_nucleus_seg::Num_Lines2(char *root, int mode) {
	int   row_no, column_no, zctr1 = 0, zctr2 = 0, zctr3 = 0;
	char  file_name[80], temp[80], ch;
	FILE  *fpinn;

	strcpy(file_name, root);
	if ((fpinn = fopen(file_name, "r")) == NULL) {
		printf("\n  The file %s does not exist !\n", file_name);
		printf("      Press <ENTER> to exit     \n\n");
		getchar();
		exit(1);
	}
	while (!feof(fpinn)) {
		fscanf(fpinn, "%c",&ch);
		zctr2++;
		//--- we don't want to count empty rows
		if (isspace(ch)) zctr2--;
		//
		if (ch=='\n' && zctr2 > 0){
			zctr3++;
			zctr2=0;
		}
		if (feof(fpinn) && zctr2 > 0){//--- this part required for linux
			zctr3++;
			zctr1++;
		}
	}
	rewind(fpinn);
	while (!feof(fpinn)){ //--- counts all data points
		fscanf(fpinn, "%s",&temp);
		zctr1++;
	}
	fclose(fpinn);
	zctr1     = zctr1-1;
	row_no    = zctr3;
	column_no = zctr1/row_no;
	//printf("\n\ncolumn no = %d  row_no= %d\n\n", column_no,   row_no);
	if (mode == 0) return row_no;
	if (mode == 1) return column_no;
	return -999;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// FILE READING/WRITING 
///////////////////////////////////////////////////////////////////////////////////////////////////////

FILE* model_nucleus_seg::FDeclare2(char *root, char *extension, char key) {
	char string1[80];
	FILE *fpot;

	strcpy(string1, root);
	if (strcmp(extension,"")!=0) strcat(string1, ".");
	strcat(string1, extension);
	if ((key != 'w')&&(key!='r')) {
		printf("either read or write, wrong key %c \n", key);getchar();
		exit(1);
	}
	if (key == 'r') {
		//if (( (FILE *) fpot = fopen(string1, "r")) == NULL) {
		//changed for DOS?UNIX compatibility
		if ((fpot = fopen(string1, "r")) == NULL) {
			printf("\n");
			printf("input file %s does not exist\n", string1);
			exit(1); getchar();
		}
	}
	if (key == 'w') {
		if ((fpot = fopen(string1, "w")) == NULL) {
			printf("\n");
			printf("output file %s does not exist\n", string1);
			exit(1);getchar();
		}
	}
	return fpot;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//CALCULATE THE AVERAGES AND STANDARD DEVIATIONS OF A SET OF VECTORS
///////////////////////////////////////////////////////////////////////////////////////////////////////

void model_nucleus_seg::calc_norm(double **mat, long D,long n, double *ave,double *stdev)
{
	double anom;		//intermediate output

	for (long i=0; i<D; i++) {
		stdev[i]=0;
		ave[i]=0;
		for (long j=0; j<n; j++) 
		{
			ave[i]+=mat[j][i];
			//std::cout<<mat[j][i] <<std::endl;
		}
		ave[i]/=n;		
		for (long j=0; j<n; j++) {
			anom=mat[j][i]-ave[i];	
			stdev[i]+=anom*anom;
		}
		stdev[i]=sqrt(stdev[i]/(n-1));

	}

}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// CONVERT THE INPUT DATA TYPE TO STRING 
// NEED TO TEMPLATE
///////////////////////////////////////////////////////////////////////////////////////////////////////
std::string model_nucleus_seg::convert2string(unsigned long id)
{
	std::stringstream out;
	out << id;
	std::string s = out.str();
	return s;
}



unsigned short model_nucleus_seg::returnthresh( OutputImageType::Pointer input_image, int num_bin_levs, int num_in_fg ){
	//Instantiate the different image and filter types that will be used
	typedef itk::ImageRegionConstIterator< InputImageType > ConstIteratorType;
	//typedef itk::Statistics::ScalarImageToHistogramGenerator< USImageType > ScalarImageToHistogramGeneratorType;
	//typedef ScalarImageToHistogramGeneratorType::HistogramType HistogramType;
	typedef itk::Statistics::Histogram< float, 1 > HistogramType;
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;

	std::cout<<"Starting threshold computation\n";

	//Create a temporary histogram container:
	const int numBins = 256;
	double tempHist[numBins];
	for(int i=0; i<numBins; ++i)
	{
		tempHist[i] = 0;
	}

	//Populate the histogram (assume pixel type is actually uchar):
	IteratorType it( input_image, input_image->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
	{
		unsigned char pix = it.Get();
		if(pix <= 255)
		{
			++tempHist[pix];
		}
	}
	
	//Find max value in the histogram
	double floatIntegerMax = itk::NumericTraits<unsigned short>::max();
	double max = 0.0;
	for(int i=0; i<numBins; ++i)
	{
		if( tempHist[i] > max )
			max = tempHist[i];
	}

	double scaleFactor = 1;
	if(max >= floatIntegerMax)
	{
		scaleFactor = floatIntegerMax / max;
	}

	HistogramType::Pointer histogram = HistogramType::New() ;
	// initialize histogram
	HistogramType::SizeType size;
	HistogramType::MeasurementVectorType lowerBound ;
	HistogramType::MeasurementVectorType upperBound ;

	lowerBound[0] = 0.0;
	upperBound[0] = (float)255.0;
	size.Fill(numBins);

	histogram->Initialize(size, lowerBound, upperBound ) ;

	int i=0;
	for (HistogramType::Iterator iter = histogram->Begin(); iter != histogram->End(); ++iter )
	{
		float norm_freq = (float)(tempHist[i] * scaleFactor);
		iter.SetFrequency(norm_freq);
		++i;
	}

	//std::cout<<"Histogram computed\n";

	//ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = ScalarImageToHistogramGeneratorType::New();
	//scalarImageToHistogramGenerator->SetNumberOfBins( 256 );
	//scalarImageToHistogramGenerator->SetInput( input_image);
	//scalarImageToHistogramGenerator->Compute();
	CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetNumberOfThresholds( num_bin_levs );
	//calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
	calculator->SetInputHistogram( histogram );
	calculator->Update();
	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput(); 
	CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();

	float thresh;

	for(int i=0; i < num_in_fg; ++itNum, ++i)
		thresh = (static_cast<float>(*itNum));

	std::cout<<"Threshold computed: "<<thresh<<std::endl;

	return (unsigned short)(thresh+0.5);

}


//Computing the surroundedness








//}
// printing the hypothesis matrix


//			 				FILE *file; 
//							file = fopen("/Users/ragzpad/hippo.txt","a+");
//									
//							for (unsigned long r = 1; r<mhmRows.size();++r)
//							{
//								fprintf(file,"%0.3lf    ", mhmRows[r]);
//								//fprintf(file,"%0.2lf    ", mhmBRows[r]);
//									
//							}
//									fprintf(file,"%d", hypothesis.size());
//								fprintf(file,"\n");
//							
//							fclose(file);
