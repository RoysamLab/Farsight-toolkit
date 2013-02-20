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

#include "yousef_seg.h"
#include <fstream>

using namespace std;

//Constructor
yousef_nucleus_seg::yousef_nucleus_seg()
{
	dataImagePtr = NULL;
	binImagePtr = NULL;

	// Added for 16 Bit And GMM estimation:
	dataImagePtr16 = NULL;
	initial_binaryImage = NULL;


	seedImagePtr = NULL;
	logImagePtr = NULL;
	clustImagePtr = NULL;
	segImagePtr = NULL;
	mySeeds.clear();

	myConnComp = NULL;
	m_pData = NULL;
	
	//int numStacks = 0;
	//int numRows = 0;
	//int numColumns = 0;		
	
	this->useDistMap = 0;	

	autoParamEstimation = false;
}

//Destructor
yousef_nucleus_seg::~yousef_nucleus_seg()
{
	clearBinImagePtr();
	clearSeedImagePtr();
	clearLogImagePtr();
	clearSegImagePtr();
	clearClustImagePtr();
	clearMyConnComp();
	mySeeds.clear();
}

//*******************************************************************************
// Set Functions
//*******************************************************************************
void yousef_nucleus_seg::setParams(int *params)
{
	shift = *params; params++;
	adaptive_bin = *params; params++;
	sigma = *params; params++;
	scaleMin = *params; params++;
	scaleMax = *params; params++;
	regionXY = *params; params++;
	regionZ = *params; params++;
	finalizeSegmentation = *params; params++;
	sampling_ratio_XY_to_Z = *params; params++;
	useDistMap = *params; params++;
	refineRange = *params; params++;
	minObjSize	= *params;
}

void yousef_nucleus_seg::setParamsForSeedDetection(int highsensitivity, double sMin, double sMax, double rXY,  double rZ, int usedistMap, int samplingRatio, int minSize)
{
	shift = highsensitivity;
	scaleMin = sMin;
	scaleMax = sMax;
	regionXY = rXY;
	regionZ = rZ;
	useDistMap = usedistMap;
	sampling_ratio_XY_to_Z = samplingRatio;
	minObjSize = minSize;
}


void yousef_nucleus_seg::setDataImage( unsigned char *imgPtr,  int x, int y, int z, const char *filename )
{
	numStacks = z;
	numRows = y;//y();			//y-direction
	numColumns = x; 		//x-direction

	dataFilename = filename;

	//delete[] dataImagePtr; I will never delete the dataImagePtr because it is created outside this class
	dataImagePtr = imgPtr;
}
void yousef_nucleus_seg::setDataImage( unsigned short *imgPtr,  int x, int y, int z, const char *filename )
{
	numStacks = z;
	numRows = y;//y();			//y-direction
	numColumns = x; 		//x-direction

	dataFilename = filename;
	dataImagePtr16 = imgPtr;
}
void yousef_nucleus_seg::setInitialBinaryImage(unsigned short* imgPtr)
{
	initial_binaryImage = imgPtr;
}

//********************************************************************************************
// get Functions
//********************************************************************************************
std::vector<int> yousef_nucleus_seg::getImageSize()
{
	std::vector<int> retVal(3);
	retVal[0] = numStacks;
	retVal[1] = numRows;
	retVal[2] = numColumns;

	return retVal;
}


//*********************************************************************************************
// internal module functions
//*********************************************************************************************
void yousef_nucleus_seg::runGradAnisDiffSmoothing()
{
	//First check to be sure that we have a dataImage to use
	if (!dataImagePtr)
		return;
	
	//Now clear all subsequent variables 
	numConnComp = 0;
	clearBinImagePtr();
	clearSeedImagePtr();
	clearLogImagePtr();
	clearSegImagePtr();
	clearClustImagePtr();
	clearMyConnComp();
	mySeeds.clear();

	std::cout<<"Starting anisotropic diffusion...";
	int ok = runGrAnisDiff(dataImagePtr, numRows, numColumns, numStacks, 5, .2, 2);
	
	if(ok)
	{
		std::cout<<"done\n";
	}
	else
	{
		std::cout<<"failed!\n segmentation will be applied on the raw image\n";
	}
}
bool yousef_nucleus_seg::RunKmeansClustering()
{
	if ( !dataImagePtr16 )
		return false;

	gmm_params_.background_mean = 0.0;
	gmm_params_.background_stdev = 0.0;
	gmm_params_.foreground_mean = 0.0;
	gmm_params_.foreground_stdev = 0.0;
	gmm_params_.background_prior = 0.0;


	// Setup Label Array for Class Labels
	unsigned short *labelArray;
	labelArray = (unsigned short *) malloc(numRows*numColumns*numStacks*sizeof(unsigned short));
	std::cout<<"trying to allocate memory to Label Array..\n";
	if( labelArray == NULL )
	{
		std::cout<<"Memory allocation for image failed\n";
		return false;
	}
	memset(labelArray/*destination*/,0/*value*/,numRows*numColumns*numStacks*sizeof(unsigned short)/*num bytes to move*/);

	// Start Kmeans:
	float mean_bg = 10.0;
	float mean_fg = 2000.0;

	bool converged = false;
	unsigned short num_iter = 0;
	std::cout<<"About Start Kmeans Iteration\n";

	while(!converged)
	{
		num_iter += 1; 
		unsigned short numBackPix = 0;
		unsigned short numForgPix = 0;
		float sumBackPix = 0.0;
		float sumForgPix = 0.0;

		for(unsigned short index = 0; index < (numRows*numColumns*numStacks);++index)
		{
			float dist_to_bg =  fabs((float)dataImagePtr16[index] - mean_bg) ;
			float dist_to_fg =  fabs((float)dataImagePtr16[index] - mean_fg) ;
			if(dist_to_bg < dist_to_fg)
			{
				labelArray[index] = 1;			// pixel is background
				numBackPix += 1;
				sumBackPix += (float)dataImagePtr16[index];
			}
			else
			{
				labelArray[index] = 2;			// pixel is foreground
				numForgPix +=1;
				sumForgPix += (float)dataImagePtr16[index];
			}
		}
		
		// add checking for number of pixels later:

		float updated_mean_bg = sumBackPix/numBackPix;
		float updated_mean_fg = sumForgPix/numForgPix;

		float changeBg = fabs(updated_mean_bg -mean_bg);
		float changeFg = fabs(updated_mean_fg -mean_fg);

		if(changeBg<0.02 || changeFg<0.02)
		{
			converged = true;
		}
		else
		{
			mean_bg = updated_mean_bg;
			mean_fg = updated_mean_fg;
		}
	}
	std::cout<<"Converged After: "<<num_iter<<" iterations\n";

	// set mean parameters:
	gmm_params_.background_mean = mean_bg;
	gmm_params_.foreground_mean = mean_fg;

	// compute the other parameters:
	unsigned long long numForgPix = 0;
	unsigned long long numBackPix = 0;

	float intensityForgSum2 = 0.0;
	float intensityBackSum2 = 0.0;

	// Compute  Std of the two Gaussians
	for(unsigned short index = 0; index < (numRows*numColumns*numStacks);++index)
	{
		if(labelArray[index]==2)
		{
			numForgPix +=1;
			intensityForgSum2 +=  (float) dataImagePtr16[index]*dataImagePtr16[index];	// E{x^2}
		}
		else
		{
			intensityBackSum2 +=  (float) dataImagePtr16[index]*dataImagePtr16[index];	// E{x^2}
		}
	}

	numBackPix = (unsigned long long)(numRows*numColumns*numStacks) - numForgPix;
	gmm_params_.foreground_stdev  = sqrt(((double)intensityForgSum2/(double)numForgPix) - (gmm_params_.foreground_mean*gmm_params_.foreground_mean));
	gmm_params_.background_stdev  = sqrt(((double)intensityBackSum2/(double)numBackPix) - (gmm_params_.background_mean*gmm_params_.background_mean));

	gmm_params_.foreground_prior = (double)numForgPix/(double)(numForgPix+numBackPix);
	gmm_params_.background_prior = (double)numBackPix/(double)(numForgPix+numBackPix);

	std::cout<< "fg_mean:"<<gmm_params_.foreground_mean<<std::endl;
	std::cout<< "fg_stdev:"<<gmm_params_.foreground_stdev<<std::endl;
	std::cout<< "fg_prior:"<<gmm_params_.foreground_prior<<std::endl;
	std::cout<< "\n";
	std::cout<< "bg_mean:"<<gmm_params_.background_mean<<std::endl;
	std::cout<< "bg_stdev:"<<gmm_params_.background_stdev<<std::endl;	
	std::cout<< "bg_prior:"<<gmm_params_.background_prior<<std::endl;

	if(gmm_params_.foreground_prior < 0.0005)
	{
		std::cout<<"Check Initial Binarization Because the foreground prior is about:"<<gmm_params_.foreground_prior<<std::endl;
		return false;
	}


	delete labelArray;
	return true;

}
bool yousef_nucleus_seg::EstimateGMMParameters()
{

	//if ( !dataImagePtr16 || !initial_binaryImage )
	if ( !dataImagePtr16 )
		return false;

	gmm_params_.background_mean = 0.0;
	gmm_params_.background_stdev = 0.0;
	gmm_params_.foreground_mean = 0.0;
	gmm_params_.foreground_stdev = 0.0;
	gmm_params_.background_prior = 0.0;



	//unsigned long long num_foreground_pix = 0;
	//unsigned long long num_background_pix = 0;

	//unsigned long long sum_foreground_intensity = 0;
	//unsigned long long sum_background_intensity = 0;

	//unsigned long long sum2_foreground_intensity = 0;
	//unsigned long long sum2_background_intensity = 0;


	//// Compute Mean and Std of the two Gaussians
	//for(unsigned int index = 0; index < (numRows*numColumns*numStacks);++index)
	//{
	//	if(initial_binaryImage[index]>0)
	//	{
	//		num_foreground_pix +=1;
	//		sum_foreground_intensity += (unsigned long long) dataImagePtr16[index];							// E{x}
	//		sum2_foreground_intensity +=  (unsigned long long) dataImagePtr16[index]*dataImagePtr16[index];	// E{x^2}
	//	}
	//	else
	//	{
	//		sum_background_intensity += (unsigned long long) dataImagePtr16[index];							// E{x}
	//		sum2_background_intensity +=  (unsigned long long) dataImagePtr16[index]*dataImagePtr16[index];	// E{x^2}
	//	}
	//}

	//std::cout<< "sum_foreground_intensity:"<<sum_foreground_intensity<<std::endl;
	//std::cout<< "sum2_foreground_intensity:"<<sum2_foreground_intensity<<std::endl;
	//std::cout<< "sum_background_intensity:"<<sum_background_intensity<<std::endl;
	//std::cout<< "sum2_background_intensity:"<<sum2_background_intensity<<std::endl;

	//num_background_pix = (unsigned long long)(numRows*numColumns*numStacks) - num_foreground_pix;

	//gmm_params_.foreground_mean = (double)sum_foreground_intensity/(double)num_foreground_pix;
	//gmm_params_.foreground_stdev  = sqrt(((double)sum2_foreground_intensity/(double)num_foreground_pix) - (gmm_params_.foreground_mean*gmm_params_.foreground_mean));

	//gmm_params_.background_mean = (double)sum_background_intensity/(double)num_background_pix;
	//gmm_params_.background_stdev  = sqrt(((double)sum2_background_intensity/(double)num_background_pix) - (gmm_params_.background_mean*gmm_params_.background_mean));

	//gmm_params_.foreground_prior = (double)num_foreground_pix/(double)(num_foreground_pix+num_background_pix);
	//gmm_params_.background_prior = (double)num_background_pix/(double)(num_foreground_pix+num_background_pix);

	//std::cout<< "fg_mean:"<<gmm_params_.foreground_mean<<std::endl;
	//std::cout<< "fg_stdev:"<<gmm_params_.foreground_stdev<<std::endl;
	//std::cout<< "fg_prior:"<<gmm_params_.foreground_prior<<std::endl;
	//std::cout<< "\n";
	//std::cout<< "bg_mean:"<<gmm_params_.background_mean<<std::endl;
	//std::cout<< "bg_stdev:"<<gmm_params_.background_stdev<<std::endl;	
	//std::cout<< "bg_prior:"<<gmm_params_.background_prior<<std::endl;

	//if(gmm_params_.foreground_prior < 0.0005)
	//{
	//	std::cout<<"Check Initial Binarization Because the foreground prior is about:"<<gmm_params_.foreground_prior<<std::endl;
	//	return false;
	//}

	/////////////// estimate using ITK EM://///////////
	std::cout<<"Started GMM Parameter Estimation...\n";
	unsigned int numberOfClasses = 2;
	typedef itk::Vector< double, 1 > MeasurementVectorType;
	typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
	SampleType::Pointer sample = SampleType::New();

	for(unsigned int index = 0; index < (numRows*numColumns*numStacks);++index)
	{
		sample->PushBack( (double)dataImagePtr16[index]);
	}

	 
  typedef itk::Array< double > ParametersType;
  ParametersType params( 6 );
 
  // Create the first set of initial parameters
  std::vector< ParametersType > initialParameters( numberOfClasses );
  params[0] = 10.0; // mean of background
  params[1] = 1000.0; // variance of background
  initialParameters[0] = params;
 
  // Create the second set of initial parameters
  params[0] = 500.0; // mean of foreground
  params[1] = 1000.0; // mean of background
  initialParameters[1] = params;

  typedef itk::Statistics::GaussianMixtureModelComponent< SampleType >  ComponentType;
   // Create the components
  std::vector< ComponentType::Pointer > components;
  for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
    {
    components.push_back( ComponentType::New() );
    (components[i])->SetSample( sample );
    (components[i])->SetParameters( initialParameters[i] );
    }
 
  typedef itk::Statistics::ExpectationMaximizationMixtureModelEstimator< SampleType > EstimatorType;
  EstimatorType::Pointer estimator = EstimatorType::New();
 
  estimator->SetSample( sample );
  estimator->SetMaximumIteration( 200 );
 
  itk::Array< double > initialProportions(numberOfClasses);
  initialProportions[0] = 0.8;	// background prior
  initialProportions[1] = 0.2;	// foreground prior
 
  estimator->SetInitialProportions( initialProportions );
 
  for ( unsigned int i = 0 ; i < numberOfClasses ; i++)
    {
    estimator->AddComponent( (ComponentType::Superclass*)(components[i]).GetPointer() );
    }
 
  estimator->Update();
  // Output the results
  //for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
  //  {
		//std::cout << "Cluster[" << i << "]" << std::endl;
		//std::cout << "    Parameters:" << std::endl;
		//std::cout << "         " << (components[i])->GetFullParameters()[1]
		//		  << std::endl;
		//std::cout << "    Proportion: ";
		//std::cout << "         " << estimator->GetProportions()[i] << std::endl;
  //  }
   
	gmm_params_.foreground_mean = (components[1])->GetFullParameters()[0];
	gmm_params_.foreground_stdev  = sqrt((components[1])->GetFullParameters()[1]);

	gmm_params_.background_mean = (components[0])->GetFullParameters()[0];
	gmm_params_.background_stdev  = sqrt((components[0])->GetFullParameters()[1]);

	gmm_params_.background_prior = estimator->GetProportions()[0];
	gmm_params_.foreground_prior = estimator->GetProportions()[1];
	

	std::cout<< "fg_mean:"<<gmm_params_.foreground_mean<<std::endl;
	std::cout<< "fg_stdev:"<<gmm_params_.foreground_stdev<<std::endl;
	std::cout<< "fg_prior:"<<gmm_params_.foreground_prior<<std::endl;
	std::cout<< "\n";
	std::cout<< "bg_mean:"<<gmm_params_.background_mean<<std::endl;
	std::cout<< "bg_stdev:"<<gmm_params_.background_stdev<<std::endl;	
	std::cout<< "bg_prior:"<<gmm_params_.background_prior<<std::endl;

	if(gmm_params_.foreground_prior < 0.00005)
	{
		std::cout<<"Check Initial Binarization Because the foreground prior is about:"<<gmm_params_.foreground_prior<<std::endl;
		return false;
	}
	std::cout<<"Exiting GMM Parameter Estimation...\n";

	//scanf("%*d");
	return true;
}




void yousef_nucleus_seg::runBinarization(unsigned short number_of_bins)
{
	//std::cout<<std::endl<<"RECENT CHANGES: Min object size was hard coded to 50, now the size is correctly read from the project definition, so make sure to include this parameter. The minimum number of objects in an image is hard coded to 3, this to avoid spurious sedd detection results in noise tiles";
	
	//try this for now
	//runGradAnisDiffSmoothing();

	//First check to be sure that we have a dataImage to use
	if (!dataImagePtr)
		return;
		
	//Now clear all subsequent variables (dependent upon this binary image
	std::cout << "Clearing binarization stuff" << std::endl;
	numConnComp = 0;
	clearBinImagePtr();
	clearSeedImagePtr();
	clearLogImagePtr();
	clearSegImagePtr();
	clearClustImagePtr();
	clearMyConnComp();
	mySeeds.clear();
	
	
	//***************************************************
	// TEST CODE:
	
	/*typedef unsigned char PixelType;
	typedef itk::Image< PixelType,  3 >   InputImageType;
	
	InputImageType::Pointer im;
	im = InputImageType::New();
	InputImageType::PointType origin;
    origin[0] = 0; 
    origin[1] = 0;    
	origin[2] = 0;    
    im->SetOrigin( origin );
	
    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    InputImageType::SizeType  size;
    size[0]  = numColumns;  // size along X
    size[1]  = numRows;  // size along Y
	size[2]  = numStacks;  // size along Z
	
    InputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    im->SetRegions( region );
    im->Allocate();
    im->FillBuffer(0);
	im->Update();
	
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	for(int i=0; i<numRows*numColumns*numStacks; i++)
	{		
		iterator1.Set(dataImagePtr[i]);
		++iterator1;	
	}
	
	typedef itk::ImageFileWriter< InputImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(im);
	writer->SetFileName("input_test.tif");
	writer->Update();*/
	
	//******
	//*******
	//**************

	//By Yousef on 9-3-2009
	//subtract the gradient image from the input image
	//subtractGradientImage(dataImagePtr, numRows, numColumns, numStacks, sampling_ratio_XY_to_Z);
	//allocate space for the binary image
	std::cout << "Allocating " << numStacks*numRows*numColumns*sizeof(unsigned short) / (1024.0 * 1024) << " MB of memory for binImagePtr" << std::endl;
	binImagePtr = new unsigned short[numStacks*numRows*numColumns];	

	std::cout<<"Start Binarization ..."<<std::endl;
	int ok = 0;
	if (numStacks == 0)
	{}
	/*else if (numStacks == 1)
	{
		ok = Cell_Binarization_2D(dataImagePtr,binImagePtr, numRows, numColumns, shift); //Do Binarization		
	}
	else
	{*/
		ok = Cell_Binarization_3D(dataImagePtr,binImagePtr, numRows, numColumns, numStacks, shift, 1, number_of_bins);		//Do Binarization		
	//}

	if(ok)
	{
		std::cout << "Entering getConnCompImage: HERE: min obj size: " <<minObjSize<< std::endl;
		numConnComp = getConnCompImage(binImagePtr, 26, minObjSize, numRows, numColumns, numStacks,1);			//Find connected components
		std::cout << "Entering getConnCompInfo3D" << std::endl;
		getConnCompInfo3D();																			//Populate myConnComp
		std::cout << "Cell Binarized.. with " << numConnComp << " connected components" << endl;	
	}
	else
	{
		cerr << "Binarization Failed!!" << endl;
	}
	
	/*typedef unsigned short PixelType2;
	typedef itk::Image< PixelType2,  3 >   InputImageType2;
	
	InputImageType2::Pointer im2;
	im2 = InputImageType2::New();
	InputImageType2::PointType origin2;
    origin2[0] = 0; 
    origin2[1] = 0;    
	origin2[2] = 0;    
    im2->SetOrigin( origin2 );
	
    InputImageType::IndexType start2;
    start2[0] =   0;  // first index on X
    start2[1] =   0;  // first in dex on Y    
	start2[2] =   0;  // first index on Z    
    InputImageType2::SizeType  size2;
    size2[0]  = numColumns;  // size along X
    size2[1]  = numRows;  // size along Y
	size2[2]  = numStacks;  // size along Z
	
    InputImageType2::RegionType region2;
    region2.SetSize( size2 );
    region2.SetIndex( start2 );
    
    im2->SetRegions( region2 );
    im2->Allocate();
    im2->FillBuffer(0);
	im2->Update();
	
	typedef itk::ImageRegionIteratorWithIndex< InputImageType2 > IteratorType2;
	IteratorType2 iterator12(im2,im2->GetRequestedRegion());
	for(int i=0; i<numRows*numColumns*numStacks; i++)
	{		
		iterator12.Set(binImagePtr[i] * 5000);
		++iterator12;	
	}
	
	typedef itk::ImageFileWriter< InputImageType2 > WriterType2;
	WriterType2::Pointer writer2 = WriterType2::New();
	writer2->SetInput(im2);
	writer2->SetFileName("binImage.tif");
	writer2->Update();

	std::cout << "Bin Image written to binImage.tif" << endl;*/
}
bool yousef_nucleus_seg::runBinarization16()
{
	
	//First check to be sure that we have a dataImage to use
	if (!dataImagePtr16)
	{
		std::cout<<"There is no 16-bit input\n";
		return false;
	}	
	//Now clear all subsequent variables (dependent upon this binary image
	std::cout << "Clearing binarization stuff" << std::endl;
	numConnComp = 0;
	clearBinImagePtr();
	clearSeedImagePtr();
	clearLogImagePtr();
	clearSegImagePtr();
	clearClustImagePtr();
	clearMyConnComp();
	mySeeds.clear();
	

	std::cout << "Allocating " << numStacks*numRows*numColumns*sizeof(unsigned short) / (1024.0 * 1024) << " MB of memory for binImagePtr" << std::endl;
	binImagePtr = new unsigned short[numStacks*numRows*numColumns];	

	std::cout<<"Start Binarization ..."<<std::endl;
	int ok = 0;

	ok = Cell_Binarization_3D(dataImagePtr16,binImagePtr, numRows, numColumns, numStacks,\
							  gmm_params_.background_mean,gmm_params_.background_stdev,gmm_params_.background_prior,\
							  gmm_params_.foreground_mean,gmm_params_.foreground_stdev,gmm_params_.foreground_prior);



	if(ok)
	{
		std::cout << "Entering getConnCompImage: HERE: min obj size: " <<minObjSize<< std::endl;
		numConnComp = getConnCompImage(binImagePtr, 26, minObjSize, numRows, numColumns, numStacks,1);			//Find connected components
		std::cout << "Entering getConnCompInfo3D" << std::endl;
		getConnCompInfo3D();																			//Populate myConnComp
		std::cout << "Cell Binarized.. with " << numConnComp << " connected components" << endl;	
	}
	else
	{
		cerr << "Binarization Failed!!" << endl;
		return false;
	}
	std::cout<<"I'm done with cell binarization\n";
	
	//for(int i=0; i<512*512; i++)
	//{		
	//	if(binImagePtr[i]>200)
	//	{
	//		std::cout<<"yes it is\n";
	//	}
	//	else
	//	{
	//		std::cout<<"no it is not\n";
	//	}
	
	//}
	return true;
}

void yousef_nucleus_seg::setBinImage(unsigned short* ptr)
{
	binImagePtr = ptr;
	std::cout << "Entering getConnCompImage" << std::endl;
	numConnComp = getConnCompImage(binImagePtr, 26, minObjSize, numRows, numColumns, numStacks,1);			//Find connected components
	std::cout << "Entering getConnCompInfo3D" << std::endl;
	getConnCompInfo3D();																			//Populate myConnComp
	std::cout << "Cell Binarized.. with " << numConnComp << " connected components" << endl;	
}

void yousef_nucleus_seg::runSeedDetection()
{
	//Check for required images
	if ( !dataImagePtr || !binImagePtr )
		return;

	std::cout<< scaleMin<<"\t"<< scaleMax<<"\t"<< regionXY <<"\t"<<regionZ<<std::endl;
	std::cout<< numStacks<<"\t"<< numRows<<"\t"<< numColumns <<std::endl;
	//Now clear all subsequent variables
	clearSeedImagePtr();
	clearLogImagePtr();
	clearClustImagePtr();
	clearSegImagePtr();
	mySeeds.clear();

	//allocate space for the binary image of seed points
	//allocate inside the 3-D seeds detection function in order to save memory for the intermediate steps
	//seedImagePtr = 0;
	//seedImagePtr = new unsigned short[numStacks*numRows*numColumns];
	
	//copy the binary image into the seeds image for now
	//memcpy(seedImagePtr/*destination*/, binImagePtr/*source*/, numStacks*numRows*numColumns*sizeof(int)/*num bytes to move*/);

	//need to pass a float pointer with input image in it, so create it here
	float *imgPtr = new float[numStacks*numRows*numColumns];
	ucharToFloat(dataImagePtr /*from*/, imgPtr /*to*/, numRows, numColumns, numStacks, 1 /*invert*/);
	//ushortToFloat(binImagePtr /*from*/, imgPtr /*to*/, numRows, numColumns, numStacks, 1 /*invert*/);

	//allocate space for the laplacian of gaussian
	//allocate inside the 3-D seeds detection function in order to save memory for the intermediate steps
	//logImagePtr = new float[numStacks*numRows*numColumns];
	
	//Now do seed detection
	int ok = 0;
	if (numStacks == 1)
	{		
		seedImagePtr = new unsigned short[numStacks*numRows*numColumns];		
		logImagePtr = new float[numStacks*numRows*numColumns];
		ok = detectSeeds2D( imgPtr, logImagePtr, seedImagePtr, numRows, numColumns, &scaleMin, &scaleMax, &regionXY, binImagePtr, autoParamEstimation );		
	}
	else
	{	
		minLoGImg = 10000;
		ok = Seeds_Detection_3D( imgPtr, &logImagePtr, &seedImagePtr, numRows, numColumns, numStacks, &scaleMin, &scaleMax, &regionXY, &regionZ, getSamplingRatio(), binImagePtr, useDistMap, &minLoGImg, autoParamEstimation );						
	}		
	delete [] imgPtr;	//cleanup
	if(!ok)
		cerr << "Seed detection Failed!!" << endl;
	else
		//Make sure all seeds are in foreground and extract vector of seeds
    //cout << "zackdebug: extracing seeds" << endl;
		ExtractSeeds();
	//added by Yousef on 9/2/2009
	//In case we did parameter estimation, write the parameters into a file
	if(autoParamEstimation)
	{
		//Write the automatically estimated parameters into a file
    //cout << "zackdebug: writing parameters to file" << endl;
		//writeParametersToFile();
	}
}

void yousef_nucleus_seg::runSeedDetection16()
{
	//Check for required images
	if ( !dataImagePtr16 || !binImagePtr )
		return;

	std::cout<< scaleMin<<"\t"<< scaleMax<<"\t"<< regionXY <<"\t"<<regionZ<<std::endl;
	std::cout<< numStacks<<"\t"<< numRows<<"\t"<< numColumns <<std::endl;
	//Now clear all subsequent variables
	clearSeedImagePtr();
	clearLogImagePtr();
	clearClustImagePtr();
	clearSegImagePtr();
	mySeeds.clear();

	//allocate space for the binary image of seed points
	//allocate inside the 3-D seeds detection function in order to save memory for the intermediate steps
	//seedImagePtr = 0;
	//seedImagePtr = new unsigned short[numStacks*numRows*numColumns];
	
	//copy the binary image into the seeds image for now
	//memcpy(seedImagePtr/*destination*/, binImagePtr/*source*/, numStacks*numRows*numColumns*sizeof(int)/*num bytes to move*/);

	//need to pass a float pointer with input image in it, so create it here
	float *imgPtr = new float[numStacks*numRows*numColumns];
	//ucharToFloat(dataImagePtr /*from*/, imgPtr /*to*/, numRows, numColumns, numStacks, 1 /*invert*/);
	ushortToFloat(binImagePtr /*from*/, imgPtr /*to*/, numRows, numColumns, numStacks, 1 /*invert*/);

	//allocate space for the laplacian of gaussian
	//allocate inside the 3-D seeds detection function in order to save memory for the intermediate steps
	//logImagePtr = new float[numStacks*numRows*numColumns];
	
	//Now do seed detection
	int ok = 0;
	if (numStacks == 1)
	{		
		seedImagePtr = new unsigned short[numStacks*numRows*numColumns];		
		logImagePtr = new float[numStacks*numRows*numColumns];
		ok = detectSeeds2D( imgPtr, logImagePtr, seedImagePtr, numRows, numColumns, &scaleMin, &scaleMax, &regionXY, binImagePtr, autoParamEstimation );		
	}
	else
	{	
		std::cout<<"seed detection is not supposed to work for 3-D images\n";
		//minLoGImg = 10000;
		//ok = Seeds_Detection_3D( imgPtr, &logImagePtr, &seedImagePtr, numRows, numColumns, numStacks, &scaleMin, &scaleMax, &regionXY, &regionZ, getSamplingRatio(), binImagePtr, useDistMap, &minLoGImg, autoParamEstimation );						
	}		
	delete [] imgPtr;	//cleanup
	if(!ok)
		cerr << "Seed detection Failed!!" << endl;
	else
		//Make sure all seeds are in foreground and extract vector of seeds
    //cout << "zackdebug: extracing seeds" << endl;
		ExtractSeeds();
	//added by Yousef on 9/2/2009
	//In case we did parameter estimation, write the parameters into a file
	if(autoParamEstimation)
	{
		//Write the automatically estimated parameters into a file
    //cout << "zackdebug: writing parameters to file" << endl;
		//writeParametersToFile();
	}
}





void yousef_nucleus_seg::runSeedDetection(int minScale,int maxScale)
{
	//Check for required images
	if ( !dataImagePtr || !binImagePtr )
		return;

	//Now clear all subsequent variables
	clearSeedImagePtr();
	clearLogImagePtr();
	clearClustImagePtr();
	clearSegImagePtr();
	mySeeds.clear();
	
	scaleMin = static_cast<double>(minScale);
	scaleMax = static_cast<double>(maxScale);
	autoParamEstimation = false;
	//allocate space for the binary image of seed points
	//allocate inside the 3-D seeds detection function in order to save memory for the intermediate steps
	//seedImagePtr = 0;
	//seedImagePtr = new unsigned short[numStacks*numRows*numColumns];
	
	//copy the binary image into the seeds image for now
	//memcpy(seedImagePtr/*destination*/, binImagePtr/*source*/, numStacks*numRows*numColumns*sizeof(int)/*num bytes to move*/);

	//need to pass a float pointer with input image in it, so create it here
	float *imgPtr = new float[numStacks*numRows*numColumns];
	ucharToFloat(dataImagePtr /*from*/, imgPtr /*to*/, numRows, numColumns, numStacks, 1 /*invert*/);
	//ushortToFloat(binImagePtr /*from*/, imgPtr /*to*/, numRows, numColumns, numStacks, 1 /*invert*/);

	//allocate space for the laplacian of gaussian
	//allocate inside the 3-D seeds detection function in order to save memory for the intermediate steps
	//logImagePtr = new float[numStacks*numRows*numColumns];
	
	//Now do seed detection
	int ok = 0;
	if (numStacks == 1)
	{		
		seedImagePtr = new unsigned short[numStacks*numRows*numColumns];		
		logImagePtr = new float[numStacks*numRows*numColumns];
		ok = detectSeeds2D( imgPtr, logImagePtr, seedImagePtr, numRows, numColumns, &scaleMin, &scaleMax, &regionXY, binImagePtr, autoParamEstimation );		
	}
	else
	{	
		minLoGImg = 10000;
		ok = Seeds_Detection_3D( imgPtr, &logImagePtr, &seedImagePtr, numRows, numColumns, numStacks, &scaleMin, &scaleMax, &regionXY, &regionZ, getSamplingRatio(), binImagePtr, useDistMap, &minLoGImg, autoParamEstimation );						
	}		
	delete [] imgPtr;	//cleanup
	if(!ok)
		cerr << "Seed detection Failed!!" << endl;
	else
		//Make sure all seeds are in foreground and extract vector of seeds
    //cout << "zackdebug: extracing seeds" << endl;
		ExtractSeeds();
	//added by Yousef on 9/2/2009
	//In case we did parameter estimation, write the parameters into a file
	if(autoParamEstimation)
	{
		//Write the automatically estimated parameters into a file
    //cout << "zackdebug: writing parameters to file" << endl;
		//writeParametersToFile();
	}
}











void yousef_nucleus_seg::runClustering()
{
	//fitMixGaussians();
	//return;
	
	//runGradWeightedDistance();
	//return;

	//Check for required images
	if( !dataImagePtr || !logImagePtr || !seedImagePtr || !binImagePtr )
		return;

	//Now clear all subsequent variables				
	clearSegImagePtr();
	clearClustImagePtr();

	//Allocate space
	clustImagePtr = new unsigned short[numStacks*numRows*numColumns];

	/*if (numStacks == 1)
	{
		std::cout << "Starting Initial Clustering" << std::endl;
		//ExtractSeeds();
		int *seed_xmclust, *seed_ymclust;
		int numseedsmclust = (int)mySeeds.size();
		seed_xmclust = (int *) malloc(mySeeds.size()*sizeof(int));
		seed_ymclust = (int *) malloc(mySeeds.size()*sizeof(int));
		for (int i=0; i<((int)mySeeds.size()); ++i)
		{
			seed_ymclust[i] = mySeeds[i].y();
			seed_xmclust[i] = mySeeds[i].x();
		}
		local_max_clust_2D(logImagePtr, numRows, numColumns, regionXY, clustImagePtr, seed_xmclust, seed_ymclust, numseedsmclust, binImagePtr);
		free( seed_xmclust );
		free( seed_ymclust );
	}
	else
	{*/
		std::cout << "Starting Initial Clustering" << std::endl;
		std::cout << "scale_xy = " << regionXY << std::endl;
		std::cout << "scale_z = " << regionZ << std::endl;
		local_max_clust_3D(logImagePtr/*LoG*/, seedImagePtr/*local max vals*/, binImagePtr/*binary mask*/,clustImagePtr/*output*/,numRows, numColumns, numStacks, regionXY, regionZ);		
	//}	
}

void yousef_nucleus_seg::runClustering16()
{
	//fitMixGaussians();
	//return;
	
	//runGradWeightedDistance();
	//return;

	//Check for required images
	if( !dataImagePtr16 || !logImagePtr || !seedImagePtr || !binImagePtr )
		return;

	//Now clear all subsequent variables				
	clearSegImagePtr();
	clearClustImagePtr();

	//Allocate space
	clustImagePtr = new unsigned short[numStacks*numRows*numColumns];

	/*if (numStacks == 1)
	{
		std::cout << "Starting Initial Clustering" << std::endl;
		//ExtractSeeds();
		int *seed_xmclust, *seed_ymclust;
		int numseedsmclust = (int)mySeeds.size();
		seed_xmclust = (int *) malloc(mySeeds.size()*sizeof(int));
		seed_ymclust = (int *) malloc(mySeeds.size()*sizeof(int));
		for (int i=0; i<((int)mySeeds.size()); ++i)
		{
			seed_ymclust[i] = mySeeds[i].y();
			seed_xmclust[i] = mySeeds[i].x();
		}
		local_max_clust_2D(logImagePtr, numRows, numColumns, regionXY, clustImagePtr, seed_xmclust, seed_ymclust, numseedsmclust, binImagePtr);
		free( seed_xmclust );
		free( seed_ymclust );
	}
	else
	{*/
		std::cout << "Starting Initial Clustering" << std::endl;
		std::cout << "scale_xy = " << regionXY << std::endl;
		std::cout << "scale_z = " << regionZ << std::endl;
		local_max_clust_3D(logImagePtr/*LoG*/, seedImagePtr/*local max vals*/, binImagePtr/*binary mask*/,clustImagePtr/*output*/,numRows, numColumns, numStacks, regionXY, regionZ);		
	//}	
}

//added by Yousef on 09/06/2009
void yousef_nucleus_seg::fitMixGaussians()
{
	//First check for necessary prerequisites
	if( !dataImagePtr || !seedImagePtr || !binImagePtr)
	{
		return;
	}

	//Now clear all subsequent variables				
	clearClustImagePtr();
	clustImagePtr = new unsigned short[numStacks*numRows*numColumns];
	memset(clustImagePtr/*destination*/,0/*value*/,numStacks*numRows*numColumns*sizeof(unsigned short)/*num bytes to move*/);
	
	//std::cerr<<"Finalizing Segmentation"<<std::endl;

	//First, add minimum plus 1 to the LoG image to insure that the minimum is 1
	//but we need to check if the minimum is negative first
	if(minLoGImg<=0)
	{
		minLoGImg = -minLoGImg;
		for(int i=0; i<numStacks*numRows*numColumns; i++)
			logImagePtr[i]+= (minLoGImg+1);
	}

	//Now, we apply the next steps into the connected components one by one
	int ind, x_len, y_len, z_len, val;
	//int min_lbl, max_lbl;	

	int objects_count = 0;
	for(int n=0; n<numConnComp; n++)
	{
		std::cout<<"Processing Connected Component #"<<n+1<<"...";
		//Now, get the subimages (the bounding box) for the current connected component
		ind = 0;
		x_len = myConnComp[n].x2 - myConnComp[n].x1 + 1;
		y_len = myConnComp[n].y2 - myConnComp[n].y1 + 1;
		z_len = myConnComp[n].z2 - myConnComp[n].z1 + 1;

		//first, get the number of seeds in this connected component
		int num_seeds_cc = 0;
		int num_points_cc = 0;
		for(int k=myConnComp[n].z1; k<=myConnComp[n].z2; k++)		 
		{			
			for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++)
			{				
				for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++)
				{					
					val = binImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];
					if(val != (n+1))
					{
						continue;
					}
					num_points_cc++;
					val = seedImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];
					if(val > 0)
						num_seeds_cc++;
				}
			}
		}
		//if you have only one or no seeds, then just set the object to the connected component
		if(num_seeds_cc<2)
		{	
			objects_count++;
			for(int k=myConnComp[n].z1; k<=myConnComp[n].z2; k++)		 
			{			
				for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++)
				{				
					for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++)
					{					
						val = binImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];
						if(val != (n+1))
						{
							continue;
						}
						clustImagePtr[(k*numRows*numColumns)+(j*numColumns)+i] = objects_count;
					}
				}
			}			
			std::cout<<"done"<<std::endl;
			continue;
		}



		//float* sublogImg = new float[x_len*y_len*z_len];
		//unsigned short* subclustImg = new unsigned short[x_len*y_len*z_len];	
		//std::vector<int> labelsList;
		std::vector<std::vector<double> > X;
		int** seeds_cc = new int*[num_seeds_cc];
		int seeds_count_cc=0;
		int new_obj_count = 0;
		for(int k=myConnComp[n].z1; k<=myConnComp[n].z2; k++)		 //come back
		{
			for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++)
			{				
				for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++)
				{					
					val = binImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];
					//The bounding box could contain points from other neighbor connected components which need to be removed
					if(val != (n+1))
					{
						continue;
					}
					std::vector<double> xx;
					xx.push_back(i);
					xx.push_back(j);
					xx.push_back(k);
					xx.push_back(dataImagePtr[(k*numRows*numColumns)+(j*numColumns)+i]);
					X.push_back(xx);
					val = seedImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];
					if(val > 0)
					{
						seeds_cc[seeds_count_cc] = new int[3];
						seeds_cc[seeds_count_cc][0] = i;
						seeds_cc[seeds_count_cc][1] = j;
						seeds_cc[seeds_count_cc][2] = k;
						seeds_count_cc++;
					}
				}
			}
		}
		//fit a mixture of gaussians
		EM_Gmm(&X, seeds_cc,num_points_cc, num_seeds_cc);
		//read point to cell assigenement values
		for(int i=0; i<num_points_cc; i++)
		{
			int x_cc = (int) X[i][0];
			int y_cc = (int) X[i][1];
			int z_cc = (int) X[i][2];
			int ass_cc = (int) X[i][3];
			ass_cc += objects_count;
			clustImagePtr[(z_cc*numRows*numColumns)+(y_cc*numColumns)+x_cc] = ass_cc;
			if(ass_cc > new_obj_count)
				new_obj_count = ass_cc;
		}
		X.empty();
		objects_count = new_obj_count;
		std::cout<<"done"<<std::endl;
	}
	std::cout<<"Initial segmentation done with "<<objects_count<<" objects"<<std::endl;
}
/////////
void yousef_nucleus_seg::ExtractSeeds()
{
	mySeeds.clear();

	unsigned short seedVal;
	unsigned short binVal;
	int curNode;
	int id = 1;

	for (int k=0; k<numStacks; ++k)
	{
		for (int j=0; j<numRows; ++j)
		{
			for (int i=0; i<numColumns; ++i)
			{
				curNode = (k*numRows*numColumns)+(j*numColumns)+i;
				seedVal = seedImagePtr[curNode];
				binVal = binImagePtr[curNode];

				if (seedVal > 0)		//meaning a local maximum here (seed exists)
				{
					if(binVal == 0)		//I'm in the background
					{
						//Set this pixel to a -1
						//seedImagePtr[curNode] = -1;
						seedImagePtr[curNode] = 65535;
					}
					else
					{
						//Keep the pixel, put it in mySeeds, and change value to id
						mySeeds.push_back(Seed(i,j,k,id,binVal));
						seedImagePtr[curNode] = id;
						//std::cerr << "seed " << id << " Added at x=" << i << " y=" << j << " z=" << k << " cc = " << binVal << std::endl;
						++id;
					}
				}
			}
		}
	}
	std::cout << id-1 << " HERE: seeds were detected"<<std::endl;

}

void yousef_nucleus_seg::outputSeeds(void)
{
	if(mySeeds.size() <= 0)
		return;

	int len = (int)dataFilename.length();
    len = len-4;
	std::string outFName = dataFilename.substr(0,len);
	outFName = outFName + "_seedPoints.txt";
	FILE* fid = fopen(outFName.c_str(),"w");

	for (unsigned int i=0; i<mySeeds.size(); ++i)
	{
		fprintf( fid,"%d %d %d\n",mySeeds[i].x(),mySeeds[i].y(),mySeeds[i].z() );
	}
	fclose(fid);
}


int yousef_nucleus_seg::getConnCompImage(unsigned short *IM, int connectivity, int minSize, size_t r, size_t c, size_t z, int runConnComp)
{
	typedef unsigned short InputPixelType;
	typedef unsigned short OutputPixelType;
	typedef itk::Image< InputPixelType,  3 >   InputImageType;
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType;
	
	InputImageType::Pointer im;
	im = InputImageType::New();
	InputImageType::PointType origin;
    origin[0] = 0; 
    origin[1] = 0;    
	origin[2] = 0;    
    im->SetOrigin( origin );

    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    InputImageType::SizeType  size;
    size[0]  = c;  // size along X
    size[1]  = r;  // size along Y
	size[2]  = z;  // size along Z
  
    InputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    im->SetRegions( region );
    im->Allocate();
    im->FillBuffer(0);
	im->Update();
	
	//copy the input image into the ITK image
	std::cout << "Copying input image into ITKImage" << std::endl;
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	for(size_t i=0; i<r*c*z; i++)
	{		
		iterator1.Set(IM[i]);
		++iterator1;	
	}
	
	
	//typedef itk::ImageFileWriter< InputImageType > WriterType;
	//WriterType::Pointer writer = WriterType::New();
	//writer->SetInput(im);
	//writer->SetFileName("im_image.tif");
	//try
	//{
	//	writer->Update();
	//}
	//catch (itk::ExceptionObject &err)
	//{
	//	std::cerr << "Error in ImageFileWriter: " << err << std::endl;
	//	return -1;
	//}

	typedef itk::ConnectedComponentImageFilter< InputImageType, OutputImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelType;
	RelabelType::Pointer relabel = RelabelType::New();

	if(runConnComp == 1)
	{
		std::cout << "Computing the labeled connected component image" << std::endl;
		//Compute the labeled connected component image		
		filter->SetInput (im);
		filter->SetFullyConnected( connectivity );		
		//use the connected component image as the input to the relabel component filter		
		relabel->SetInput( filter->GetOutput() );
	}
	else
	{
		//use the input image as the input to the relabel component filter 		
		//relabel->SetInput( im );
	}

	std::cout << "Setting minimum object size" << std::endl;
	//set the minimum object size
	relabel->SetMinimumObjectSize( minSize );

	try
    {
		relabel->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
		std::cerr << "Relabel: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return -1;
    }

	
    //write the output of the labeling CC filter into our input image
	std::cout << "Writing output of relabel coloring filter to input image" << std::endl;
	// At least 3 connected components in the image, otherwise for the noise tiles, the seed detection will end up detecting many seeds. This parame
	int numObj = relabel->GetNumberOfObjects();
	if( numObj > 2 )
	{
		int mx;
		IteratorType iterator2(relabel->GetOutput(),relabel->GetOutput()->GetRequestedRegion());
		for(size_t i=0; i<r*c*z; i++)
		{			
			
			mx = iterator2.Get();	
			if(mx == -1)
				IM[i] = 0;
			else
				IM[i] = mx;
			++iterator2;	
		}
	}
	else
	{
		for( size_t i=0; i<r*c*z; ++i )
		{
			IM[i] = 0;
		}
	}
			

	//return the number of CCs
	std::cout << "Returning number of connected components" << std::endl;
	return numObj;
}

//By Yousef: 5-21-2008
//Get the connected component information of the binary (conncomp) image
void yousef_nucleus_seg::getConnCompInfo3D()
{
	int val; 
	size_t curNode;
	myConnComp = new ConnComp[numConnComp];
	for( int i=0; i<numConnComp; i++){
		myConnComp[i].x1 = 	numColumns+1;
		myConnComp[i].y1 =  numRows+1;
		myConnComp[i].x2 =  myConnComp[i].y2 = -1;
		if (numStacks == 1){
			myConnComp[i].z1 = 1;
			myConnComp[i].z2 = 1;
		}
		else{
			myConnComp[i].z1 = numStacks+1;
			myConnComp[i].z2 = -1;
		}
	}
	if (numStacks == 1){
		for ( int j=0; j<numRows; ++j ){
			for ( int i=0; i<numColumns; ++i ){
				curNode = (j*numColumns)+i;
				val = binImagePtr[curNode];
				if(val>0){
					val = val-1;
					myConnComp[val].y1 = (int) std::min((double)j,(double)myConnComp[val].y1);
					myConnComp[val].y2 = (int) std::max((double)j,(double)myConnComp[val].y2);
					myConnComp[val].x1 = (int) std::min((double)i,(double)myConnComp[val].x1);
					myConnComp[val].x2 = (int) std::max((double)i,(double)myConnComp[val].x2);
				}
			}
		}
	}
	else
	{
		for (int k=0; k<numStacks; ++k)
		{
			for (int j=0; j<numRows; ++j)
			{
				for (int i=0; i<numColumns; ++i)
				{
					curNode = (k*numRows*numColumns)+(j*numColumns)+i;
					val = binImagePtr[curNode];
					//if(val == 14)
					//	int uudf=1;
					if(val>0)
					{
						val = val-1;
						myConnComp[val].y1 = (int) std::min((double)j,(double)myConnComp[val].y1);
						myConnComp[val].y2 = (int) std::max((double)j,(double)myConnComp[val].y2);
						myConnComp[val].x1 = (int) std::min((double)i,(double)myConnComp[val].x1);
						myConnComp[val].x2 = (int) std::max((double)i,(double)myConnComp[val].x2);
						myConnComp[val].z1 = (int) std::min((double)k,(double)myConnComp[val].z1);
						myConnComp[val].z2 = (int) std::max((double)k,(double)myConnComp[val].z2);
					}
				}
			}
		}
	}
}

//THIS IS THE FINAL STAGE TO SEGMENTATION
// MAYBE SHOULD MOVE THIS FUNCTION TO THE ALPHA EXPANSION FOLDER
void yousef_nucleus_seg::runAlphaExpansion(){
	if (numStacks == 1){
		runAlphaExpansion2D();
	}
	else{
		runAlphaExpansion3D();
	}
}

void yousef_nucleus_seg::runAlphaExpansion16(){
	if (numStacks == 1){
		runAlphaExpansion2D16();
	}
	else{
		std::cout<<"does not work for 3-D data\n";
	}
}
void yousef_nucleus_seg::runAlphaExpansion2D16(){
	//First check for necessary prerequisites
	if( !dataImagePtr16 || !logImagePtr || !seedImagePtr || !binImagePtr || !myConnComp ){
		return;
	}

	//Now clear all subsequent variables
	clearSegImagePtr();

	std::cout<<"Finalizing Segmentation"<<std::endl;

	//Now, we apply the next steps into the connected components one by one
	int ind, x_len, y_len, val;

	segImagePtr = new unsigned short[numRows*numColumns];
	memset(segImagePtr/*destination*/,0/*value*/,numStacks*numRows*numColumns*sizeof(unsigned short)/*num bytes to move*/);

	for( int n=0; n<numConnComp; n++ )
	{
		std::cout<<"Processing Connected Component #"<<n+1<<"...";
		//Now, get the subimages (the bounding box) for the current connected component
		ind = 0;
		x_len = myConnComp[n].x2 - myConnComp[n].x1 + 1;
		y_len = myConnComp[n].y2 - myConnComp[n].y1 + 1;
		float* sublogImg = new float[x_len*y_len];
		unsigned short* subclustImg = new unsigned short[x_len*y_len];	
		std::vector<int> labelsList;
		
		for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++){
			for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++){	
				val = binImagePtr[(j*numColumns)+i];
				//The bounding box could contain points from other neighbor connected components which need to be removed
				if(val != (n+1)){
					sublogImg[ind] = 0;
					subclustImg[ind] = 0;
				}
				else{
					//for now, write the value in the clustering image into the segmentation image
					segImagePtr[(j*numColumns)+i] = clustImagePtr[(j*numColumns)+i];
					subclustImg[ind] = clustImagePtr[(j*numColumns)+i];
					//Do the same for the LoG image
					sublogImg[ind] = logImagePtr[(j*numColumns)+i];
					int found = 0;
					for(unsigned int l=0; l<labelsList.size(); l++){
						if(labelsList[l] == subclustImg[ind])
							found = 1;
					}
					if(found == 0) //add the label of the cluster
						labelsList.push_back(subclustImg[ind]);
				}
				++ind;
			}
		}
		//Now, if this connected component has one cell (one label) only, 
		//then take the clustering results of that connected as the final segmentation
		if(labelsList.size() == 1){
			std::cout<<"Done with only one object"<<std::endl;
			delete[] sublogImg;
			delete[] subclustImg;
			continue;
		}
		std::cout<<std::endl<<"    "<<labelsList.size()<<" objects found"<<std::endl;
		//If you reach here, it means that the current connected component contains two or more cells
		//First, sort the labels list
		std::cout<<"    "<<"sorting labels"<<std::endl;
		for(unsigned int l1=0; l1<labelsList.size(); l1++){
			for(unsigned int l2=l1+1; l2<labelsList.size(); l2++){
				if(labelsList[l2]<labelsList[l1]){
					int tmp = labelsList[l1];
					labelsList[l1] = labelsList[l2];
					labelsList[l2] = tmp;
				}
			}
		}
		

		//Relabel the clustering sub-image starting from 1
		//also get the original sub-image (bounding box) that will be used as the contrast term		
		ind = -1;		
		float* subDataImg = new float[x_len*y_len];

		for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++){
			for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++){
				ind++;
				subDataImg[ind] = (float)dataImagePtr16[(j*numColumns)+i];
				val = binImagePtr[(j*numColumns)+i];
				if(val != (n+1))
					continue;
				else{
					for(unsigned int l=0; l<labelsList.size(); l++){
						if(labelsList[l] == subclustImg[ind]){
							subclustImg[ind] = l+1;
							break;
						}
					}
				}
			}
		}

		
		alpha_expansion_2d( subDataImg, sublogImg, subclustImg, x_len, y_len );

		//relable and copy the output of the alpha expansion which is stored in the subclustImg to the final segmented image
		ind = 0;		
		for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++)
		{
			for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++)
			{		
				val = subclustImg[ind];
				if(val>0)
					segImagePtr[(j*numColumns)+i] = val;
				
				++ind;
			}
		}
		//std::cerr<<"Done with "<<labelsList.size()<<" objects"<<std::endl;
		std::cout<<"Done"<<std::endl;
		delete [] sublogImg;
		delete [] subclustImg;
		delete [] subDataImg;
	}
	//relabel the cells
	int numOfObjs = getRelabeledImage(segImagePtr, 8, minObjSize, numRows, numColumns,numStacks, 1);		
    numOfObjs--;
	std::cout << "done with " << numOfObjs<<" found"<<std::endl;
	std::cout << "Creating Final Label Image" << std::endl;	
}

void yousef_nucleus_seg::runAlphaExpansion2D(){
	//First check for necessary prerequisites
	if( !dataImagePtr || !logImagePtr || !seedImagePtr || !binImagePtr || !myConnComp ){
		return;
	}

	//Now clear all subsequent variables
	clearSegImagePtr();

	std::cout<<"Finalizing Segmentation"<<std::endl;

	//Now, we apply the next steps into the connected components one by one
	int ind, x_len, y_len, val;

	segImagePtr = new unsigned short[numRows*numColumns];
	memset(segImagePtr/*destination*/,0/*value*/,numStacks*numRows*numColumns*sizeof(unsigned short)/*num bytes to move*/);

	for( int n=0; n<numConnComp; n++ )
	{
		std::cout<<"Processing Connected Component #"<<n+1<<"...";
		//Now, get the subimages (the bounding box) for the current connected component
		ind = 0;
		x_len = myConnComp[n].x2 - myConnComp[n].x1 + 1;
		y_len = myConnComp[n].y2 - myConnComp[n].y1 + 1;
		float* sublogImg = new float[x_len*y_len];
		unsigned short* subclustImg = new unsigned short[x_len*y_len];	
		std::vector<int> labelsList;
		
		for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++){
			for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++){	
				val = binImagePtr[(j*numColumns)+i];
				//The bounding box could contain points from other neighbor connected components which need to be removed
				if(val != (n+1)){
					sublogImg[ind] = 0;
					subclustImg[ind] = 0;
				}
				else{
					//for now, write the value in the clustering image into the segmentation image
					segImagePtr[(j*numColumns)+i] = clustImagePtr[(j*numColumns)+i];
					subclustImg[ind] = clustImagePtr[(j*numColumns)+i];
					//Do the same for the LoG image
					sublogImg[ind] = logImagePtr[(j*numColumns)+i];
					int found = 0;
					for(unsigned int l=0; l<labelsList.size(); l++){
						if(labelsList[l] == subclustImg[ind])
							found = 1;
					}
					if(found == 0) //add the label of the cluster
						labelsList.push_back(subclustImg[ind]);
				}
				++ind;
			}
		}
		//Now, if this connected component has one cell (one label) only, 
		//then take the clustering results of that connected as the final segmentation
		if(labelsList.size() == 1){
			std::cout<<"Done with only one object"<<std::endl;
			delete[] sublogImg;
			delete[] subclustImg;
			continue;
		}
		std::cout<<std::endl<<"    "<<labelsList.size()<<" objects found"<<std::endl;
		//If you reach here, it means that the current connected component contains two or more cells
		//First, sort the labels list
		std::cout<<"    "<<"sorting labels"<<std::endl;
		for(unsigned int l1=0; l1<labelsList.size(); l1++){
			for(unsigned int l2=l1+1; l2<labelsList.size(); l2++){
				if(labelsList[l2]<labelsList[l1]){
					int tmp = labelsList[l1];
					labelsList[l1] = labelsList[l2];
					labelsList[l2] = tmp;
				}
			}
		}
		

		//Relabel the clustering sub-image starting from 1
		//also get the original sub-image (bounding box) that will be used as the contrast term		
		ind = -1;		
		float* subDataImg = new float[x_len*y_len];

		for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++){
			for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++){
				ind++;
				subDataImg[ind] = (float)dataImagePtr[(j*numColumns)+i];
				val = binImagePtr[(j*numColumns)+i];
				if(val != (n+1))
					continue;
				else{
					for(unsigned int l=0; l<labelsList.size(); l++){
						if(labelsList[l] == subclustImg[ind]){
							subclustImg[ind] = l+1;
							break;
						}
					}
				}
			}
		}

		
		alpha_expansion_2d( subDataImg, sublogImg, subclustImg, x_len, y_len );

		//relable and copy the output of the alpha expansion which is stored in the subclustImg to the final segmented image
		ind = 0;		
		for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++)
		{
			for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++)
			{		
				val = subclustImg[ind];
				if(val>0)
					segImagePtr[(j*numColumns)+i] = val;
				
				++ind;
			}
		}
		//std::cerr<<"Done with "<<labelsList.size()<<" objects"<<std::endl;
		std::cout<<"Done"<<std::endl;
		delete [] sublogImg;
		delete [] subclustImg;
		delete [] subDataImg;
	}
	//relabel the cells
	int numOfObjs = getRelabeledImage(segImagePtr, 8, minObjSize, numRows, numColumns,numStacks, 1);		
    numOfObjs--;
	std::cout << "done with " << numOfObjs<<" found"<<std::endl;
	std::cout << "Creating Final Label Image" << std::endl;	
}
void yousef_nucleus_seg::runAlphaExpansion3D()
{
	//First check for necessary prerequisites
	if( !dataImagePtr || !logImagePtr || !seedImagePtr || !binImagePtr || !myConnComp )
	{
		return;
	}

	//Now clear all subsequent variables				
	clearSegImagePtr();
	//added by Yousef on 8/28/2009
	//We no longer need to use the seeds image
	clearSeedImagePtr(); //This assumes that you can't go back to run clustering again in the segmentation wizard!

	std::cout<<"Finalizing Segmentation"<<std::endl;

	//by yousef on 9/2/2009
	int maxNumColors = 0;

	//First, add minimum plus 1 to the LoG image to insure that the minimum is 1
	//but we need to check if the minimum is negative first
	if(minLoGImg<=0)
	{
		minLoGImg = -minLoGImg;
		for(int i=0; i<numStacks*numRows*numColumns; i++)
			logImagePtr[i]+= (minLoGImg+1);
	}

	//Now, we apply the next steps into the connected components one by one
	int ind, x_len, y_len, z_len, val;
	//int min_lbl, max_lbl;
	segImagePtr = new unsigned short[numStacks*numRows*numColumns];
	//memcpy(segImagePtr/*destination*/, clustImagePtr/*source*/, numStacks*numRows*numColumns*sizeof(int)/*num bytes to move*/);
	memset(segImagePtr/*destination*/,0/*value*/,numStacks*numRows*numColumns*sizeof(unsigned short)/*num bytes to move*/);

	for(int n=0; n<numConnComp; n++)
	{
		std::cout<<"Processing Connected Component #"<<n+1<<"...";
		//Now, get the subimages (the bounding box) for the current connected component
		ind = 0;
		x_len = myConnComp[n].x2 - myConnComp[n].x1 + 1;
		y_len = myConnComp[n].y2 - myConnComp[n].y1 + 1;
		z_len = myConnComp[n].z2 - myConnComp[n].z1 + 1;
		float* sublogImg = new float[x_len*y_len*z_len];
		unsigned short* subclustImg = new unsigned short[x_len*y_len*z_len];	
		std::vector<int> labelsList;
		
		for(int k=myConnComp[n].z1; k<=myConnComp[n].z2; k++)		
		{
			for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++)
			{				
				for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++)
				{					
					val = binImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];
					//The bounding box could contain points from other neighbor connected components which need to be removed
					if(val != (n+1))
					{
						sublogImg[ind] = 0;
						subclustImg[ind] = 0;
					}
					else
					{
						//for now, write the value in the clustering image into the segmentation image
						segImagePtr[(k*numRows*numColumns)+(j*numColumns)+i] = clustImagePtr[(k*numRows*numColumns)+(j*numColumns)+i]; 
						subclustImg[ind] = clustImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];
						//Do the same for the LoG image
						sublogImg[ind] = logImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];
						int found = 0;
						for(unsigned int l=0; l<labelsList.size(); l++)
						{
							if(labelsList[l] == subclustImg[ind])
								found = 1;
						}
						if(found == 0) //add the label of the cluster						
							labelsList.push_back(subclustImg[ind]);												
					}
					ind++;
				}
			}			
		}
				

		//Now, if this connected component has one cell (one label) only, 
		//then take the clustering results of that connected as the final segmentation
		if(labelsList.size() == 1)
		{
			std::cout<<"Done with only one object"<<std::endl;
			delete[] sublogImg;
			delete[] subclustImg;
			continue;
		}
		
		std::cout<<std::endl<<"    "<<labelsList.size()<<" objects found"<<std::endl;
		//If you reach here, it means that the current connected component contains two or more cells
		//First, sort the labels list
		std::cout<<"    "<<"sorting labels"<<std::endl;
		for(unsigned int l1=0; l1<labelsList.size(); l1++)
		{
			for(unsigned int l2=l1+1; l2<labelsList.size(); l2++)
			{
				if(labelsList[l2]<labelsList[l1])
				{
					int tmp = labelsList[l1];
					labelsList[l1] = labelsList[l2];
					labelsList[l2] = tmp;
				}
			}
		}
		//Now, create another list that holds the new IDs of the cells from 1 to the number of cells
		//std::vector<int> cellIDsList;
		//for(int l1=0; l1<labelsList.size(); l1++)
			//cellIDsList.push_back(l1+1);
		//Relabel the clustering sub-image starting from 1
		//also get the original sub-image (bounding box) that will be used as the contrast term		
		ind = -1;		
		float* subDataImg = new float[x_len*y_len*z_len];
		for(int k=myConnComp[n].z1; k<=myConnComp[n].z2; k++)		
		{
			for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++)
			{
				for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++)   							
				{
					ind++;
					subDataImg[ind] = (float)dataImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];/*(float)dataImagePtr[(k*numRows*numColumns)+(i*numRows)+j];*/
					val = binImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];
					if(val != (n+1))
					{							
						continue;
					}
					else
					{												
						for(unsigned int l=0; l<labelsList.size(); l++)
						{
							if(labelsList[l] == subclustImg[ind])
							{
								subclustImg[ind] = l+1;
								break;
							}		
						}											
					}
				}
			}
		}

		//Call the module that does graph coloring, ML estimation, and graph learning		
		int NC = 1000;
		unsigned short* subsegImg = new unsigned short[x_len*y_len*z_len];			
		float* Dterms  = multiColGraphLearning(sublogImg, subclustImg, subsegImg, y_len, x_len, z_len, &NC,refineRange);		

		if(NC>maxNumColors)
			maxNumColors = NC;

		std::cout<<"    Starting alpha-expansion..";		
		start_alpha_expansion(subDataImg, subsegImg, Dterms, y_len, x_len, z_len, NC+1);										

		//relable and copy the output of the alpha expansion to the segmentaion image
		ind = 0;		
		for(int k=myConnComp[n].z1; k<=myConnComp[n].z2; k++)		
		{
			for(int j=myConnComp[n].y1; j<=myConnComp[n].y2; j++)												
			{
				for(int i=myConnComp[n].x1; i<=myConnComp[n].x2; i++)				
				{
					val = subsegImg[ind];					
					if(val>0)
						segImagePtr[(k*numRows*numColumns)+(j*numColumns)+i] = val;
					
					ind++;
				}
			}
		}	
		std::cout<<"Done"<<std::endl;
		delete [] sublogImg;
		delete [] subsegImg;
		delete [] subclustImg;
		delete [] subDataImg;

		free(Dterms);
	}		

	//relabel the cells
	int numOfObjs = /*getConnCompImage*/getRelabeledImage(segImagePtr, 6, minObjSize, numRows, numColumns,numStacks, 1);			
	std::cout << "done with " << numOfObjs<<" found"<<std::endl;
	std::cout << "Creating Final Label Image" << std::endl;		
}

//Added by Yousef on 9/14/2009
//
//void yousef_nucleus_seg::runGradWeightedDistance()
//{	
//	if (!dataImagePtr)
//		return;
//
//	//computer the gradient image
//	float* gradIM = new float[numStacks*numRows*numColumns];
//	memset(gradIM/*destination*/,0.0/*value*/,numStacks*numRows*numColumns*sizeof(float)/*num bytes to move*/);
//	computeGradientImage(dataImagePtr, gradIM, numRows, numColumns, numStacks);
//	//pad with zeros	
//	float* gradIMpad;
//	if(numStacks == 1)
//	{
//		gradIMpad = new float[(numRows+2)*(numColumns+2)];
//		memset(gradIMpad/*destination*/,0.0/*value*/,(numRows+2)*(numColumns+2)*sizeof(float)/*num bytes to move*/);
//	}
//	else
//	{
//		gradIMpad = new float[(numStacks+2)*(numRows+2)*(numColumns+2)];
//		memset(gradIMpad/*destination*/,0.0/*value*/,(numStacks+2)*(numRows+2)*(numColumns+2)*sizeof(float)/*num bytes to move*/);
//	}
//
//	for(int k=0; k<numStacks; k++)
//	{
//		for(int i=0; i<numRows; i++)
//		{
//			for(int j=0; j<numColumns; j++)
//			{
//				int ind1 = k*numRows*numColumns + i*numColumns + j;
//				int ind2;
//				if(numStacks == 1)
//					ind2 = (i+1)*(numColumns+2)+j+1;
//				else
//					ind2 = (k+1)*(numRows+2)*(numColumns+2)+(i+1)*(numColumns+2)+j+1;
//				gradIMpad[ind2] = gradIM[ind1];
//			}
//		}
//	}
//	delete [] gradIM;
//
//	//Now prepare the input image to the gradient weighted distance transform
//	//In this image, each seed is set to zero, each background point is set to -inf and each foreground point is set to +inf
//	float* outImage = new float[numStacks*numRows*numColumns];
//	memset(outImage/*destination*/,0.0/*value*/,numStacks*numRows*numColumns*sizeof(float)/*num bytes to move*/);
//	prepareInputImage(binImagePtr, seedImagePtr,outImage, numRows, numColumns, numStacks);
//	
//	//finally, call the gradient-weighted distance transform function
//	gradient_enhanced_distance_map_2D( outImage, gradIMpad, numColumns, numRows);
//
//	delete [] gradIM;
//
//	//try this: copy the resulting image into the clustering image
//	for(int k=0; k<numStacks; k++)
//	{
//		for(int i=0; i<numRows; i++)
//		{
//			for(int j=0; j<numColumns; j++)
//			{
//				clustImagePtr[k*numColumns*numRows+i*numColumns+j] = (unsigned short) outImage[(k+1)*(numColumns+2)*(numRows+2)+(i+1)*(numColumns+2)+j+1];
//			}
//		}
//	}
//
//	delete [] outImage;
//}

unsigned char ***TriplePtr(int z, int r, int c)
{
	unsigned char ***ptr = new unsigned char**[z];
	for ( int j = 0; j<z; ++j )
	{
		ptr[j] = new unsigned char*[r];
		for ( int i = 0; i<r; ++i)
		{
			ptr[j][i] = new unsigned char[c];
		}
	}	
	return ptr;
}

void ucharToFloat(unsigned char* fromLoc, float* toLoc,int r, int c, int z, char invert)
{
	unsigned char val;
	size_t curNode;

	if ((toLoc != NULL) && (fromLoc != NULL))
	{
		for (size_t k=0; k<z; ++k)
		{
			for (size_t j=0; j<r; ++j)
			{
				for (size_t i=0; i<c; ++i)
				{
					curNode = (k*r*c)+(j*c)+i;
					val = fromLoc[curNode];
					if (invert == 1)
						val = 255-val;
					toLoc[curNode] = (float)val;	
				}
			}
		}
	}
	else
	{
		std::cerr << "POINTERS NOT INITIALIZED in ucharToFloat" << std::endl;
	}
}

void ushortToFloat(unsigned short* fromLoc, float* toLoc,int r, int c, int z, char invert)
{
	unsigned short val;
	int curNode;

	if ((toLoc != NULL) && (fromLoc != NULL))
	{
		for (int k=0; k<z; ++k)
		{
			for (int j=0; j<r; ++j)
			{
				for (int i=0; i<c; ++i)
				{
					curNode = (k*r*c)+(j*c)+i;
					val = fromLoc[curNode];
					if (invert == 1)
						val = 255-val;
					toLoc[curNode] = (float)val;	
				}
			}
		}
	}
	else
	{
		std::cerr << "POINTERS NOT INITIALIZED in ucharToFloat" << std::endl;
	}
}

void ucharToUshort(unsigned char* fromLoc, unsigned short* toLoc,int r, int c, int z, char invert)
{
	unsigned char val;
	int curNode;

	if ((toLoc != NULL) && (fromLoc != NULL))
	{
		for (int k=0; k<z; ++k)
		{
			for (int j=0; j<r; ++j)
			{
				for (int i=0; i<c; ++i)
				{
					curNode = (k*r*c)+(j*c)+i;
					val = fromLoc[curNode];
					if (invert == 1)
						val = 255-val;
					toLoc[curNode] = (unsigned short)val;	
				}
			}
		}
	}
	else
	{
		std::cerr << "POINTERS NOT INITIALIZED in ucharToFloat" << std::endl;
	}
}

//added by Yousef on 8-5-2008
//This function reads parameters from a .ini file following specific format
void yousef_nucleus_seg::readParametersFromFile(const char* pFname)
{
	char achBuffer[1024];
	char achBuffer2[1024];
	char* pchStr = NULL;
	int iCounter = 0;
	int m_iNumOfElements = 0;
	int params[12];//modified by yousef on 11/4/2008
	params[0]=0;	//sensitivity
	params[1]=0;	//adaptive binarization
	//added by vinay on 02/09/2012
	params[2]=30;	//loG size
	params[3]=5;	//min_scale
	params[4]=8;	//max_scale
	params[5]=5;	//xy_clustering
	params[6]=2;	//z_clust
	params[7]=0;	//finalize
	params[8]=2;	//sampling ratio
	params[9]=1;	//use_dist_map
	//added by yousef on 11/4/2008
	params[10]=6;	//refinement range
	//added by yousef on 12/5/2008
	params[11]=20;	//min_object_size
	std::ifstream inFile(pFname);
	if (! inFile)
	{
		/*cout << "Fatal Error: Could not open parameters file " << pFname
			<< " .... Terminating Program." << endl;
		exit(0);*/
		//if no file name is provided, then just use the default parameter values
		//and later replace them with the automatically estimated values.
		autoParamEstimation = true; 
		setParams(params);
		return;
	}

	while (inFile)
	{
		inFile.getline(achBuffer, 1024, '\n');
		if (achBuffer[0] != '\0' && achBuffer[0] != '!')
		{
			iCounter++;
		}
	}
	if (iCounter == 0)
	{
		cout << "Fatal Error: Empty Configuration File " << pFname
			<< " .... Terminating Program." << endl;
		exit(0);
	}

	m_pData = new TParamsEntry[iCounter];


	inFile.close();
	//inFile.open(pFname);

	std::ifstream inFile2(pFname);
	if (! inFile2)
	{
		cout << "Fatal Error: Could not open parameters file in the second time" << pFname
			<< " .... Terminating Program." << endl;
		exit(0);
	}
	while (inFile2)
	{
		inFile2.getline(achBuffer, 1024, '\n');
		if (achBuffer[0] != '\0' && achBuffer[0] != '!')
		{		
			strcpy(achBuffer2, achBuffer);
			pchStr = strtok(achBuffer, "\t ");

			if (pchStr)
			{
				strcpy(m_pData[m_iNumOfElements].m_pName, pchStr);
			}

			pchStr = strtok(NULL, "\t ");

			if (pchStr)
			{
				if (*pchStr == ':')
					pchStr = strtok(NULL, "\t " );

				m_pData[m_iNumOfElements].m_pValue= atoi(pchStr);
			}
			m_iNumOfElements++;
		}
	}
	inFile2.close();	

	for(int i=0; i<=iCounter; i++)
	{
		if(!strcmp(m_pData[i].m_pName,"high_sensitivity"))
			params[0] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"adaptive_binarization"))
			params[1] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"LoG_size"))
			params[2] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"min_scale"))
			params[3] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"max_scale"))
			params[4] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"xy_clustering_res"))
			params[5] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"z_clustering_res"))
			params[6] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"finalize_segmentation"))
			params[7] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"sampling_ratio_XY_to_Z"))
			params[8] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"Use_Distance_Map"))
			params[9] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"refinement_range"))
			params[10] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"min_object_size"))
			params[11] = m_pData[i].m_pValue;
		else
			continue;
	}

	setParams(params);

}


//added by Yousef on 9-2-2009
//This function writes the automatically estimated parameters into a text file
void yousef_nucleus_seg::writeParametersToFile()
{
	std::string pFname = dataFilename.substr(0,dataFilename.find_last_of("."))+"_parameters.ini";
	std::ofstream outFile(pFname.c_str());
	if (! outFile)
	{
		cout << "Fatal Error: Could not open parameters file " << pFname
			<< " for writing.... Terminating Program." << endl;
		exit(0);		
	}
	outFile << "! Segmentation parameters File"<<std::endl; 
	outFile << "! All parameters are case sensitive"<<std::endl;
	outFile << std::endl;
	outFile << "high_sensitivity\t:\t"<<shift<<std::endl;
	outFile << "adaptive_binarization\t:\t"<<adaptive_bin<<std::endl;
	outFile << "LoG_size\t:\t"<<sigma<<std::endl;
	outFile << "min_scale\t:\t"<<scaleMin<<std::endl;
	outFile << "max_scale\t:\t"<<scaleMax<<std::endl;
	outFile << "xy_clustering_res\t:\t"<<regionXY<<std::endl;
	outFile << "z_clustering_res\t:\t"<<regionZ<<std::endl;
	outFile << "finalize_segmentation\t:\t"<<finalizeSegmentation<<std::endl;
	outFile << "sampling_ratio_XY_to_Z\t:\t"<<sampling_ratio_XY_to_Z<<std::endl;
	outFile << "Use_Distance_Map\t:\t"<<useDistMap<<std::endl;
	outFile << "refinement_range\t:\t"<<refineRange<<std::endl;
	outFile << "min_object_size\t:\t"<<minObjSize<<std::endl;

	outFile.close();									
}

//Added by Yousef on 7-8-2008
int yousef_nucleus_seg::getRelabeledImage(unsigned short *IM, int connectivity, int minSize, int r, int c, int z, int runConnComp)
{
	typedef    unsigned short     InputPixelType;
	typedef    int     OutputPixelType;
	typedef itk::Image< InputPixelType,  3 >   InputImageType;
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType;
	

	InputImageType::Pointer im;
	im = InputImageType::New();
	InputImageType::PointType origin;
    origin[0] = 0; 
    origin[1] = 0;    	
	origin[2] = 0;    
    im->SetOrigin( origin );

    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    	
	start[2] =   0;  // first index on Z    
    InputImageType::SizeType  size;
    size[0]  = c;  // size along X
    size[1]  = r;  // size along Y
	size[2]  = z;  // size along Z

    InputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );

    im->SetRegions( region );
    im->Allocate();
    im->FillBuffer(0);
	im->Update();
	
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	for(int i=0; i<r*c*z; i++)
	{
		iterator1.Set(IM[i]);
		++iterator1;	
	}
	
	typedef itk::ScalarConnectedComponentImageFilter< InputImageType, OutputImageType > FilterType;	
	FilterType::Pointer filter = FilterType::New();
	typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelType;
	RelabelType::Pointer relabel = RelabelType::New();

	//if(runConnComp == 1)
	//{				
		//Compute the labeled connected component image		
		filter->SetInput (im);
		filter->SetFullyConnected( connectivity );					
		//use the connected component image as the input to the relabel component filter		
		relabel->SetInput( filter->GetOutput() );
	//}
	//else
	//{
		//use the input image as the input to the relabel component filter 		
		//relabel->SetInput( im );
	//}

    //set the minimum object size
		std::cout<<"Min obj size: "<<minSize<<std::endl;
	relabel->SetMinimumObjectSize( minSize );

	try
    {
		relabel->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
		std::cerr << "Relabel: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
    }
	
	//write the output of the labeling CC filter into our input image
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > IteratorType2;
	IteratorType2 iterator2(relabel->GetOutput(),relabel->GetOutput()->GetRequestedRegion());	
	for(int i=0; i<r*c*z; i++)
	{		
		if(binImagePtr[i] == 0)
			IM[i] = 0;
		else
		{
			int lbl = iterator2.Get()-1;		
			if(lbl < 0)
				IM[i]=0;
			else
				IM[i] = lbl;

		}		
		++iterator2;	
	}

	//return the number of CCs
	return relabel->GetNumberOfObjects()-1;    
}

void yousef_nucleus_seg::clearBinImagePtr(void)
{
	if(binImagePtr)
	{
		delete[] binImagePtr;
		binImagePtr = NULL;
	}
}
void yousef_nucleus_seg::clearSeedImagePtr(void)
{
	if(seedImagePtr)
	{
		delete[] seedImagePtr;
		seedImagePtr = NULL;
	}
}
void yousef_nucleus_seg::clearLogImagePtr(void)
{
	if(logImagePtr)
	{
		delete[] logImagePtr;
		logImagePtr = NULL;
	}
}
void yousef_nucleus_seg::clearSegImagePtr(void)
{
	if(segImagePtr)
	{
		delete[] segImagePtr;
		segImagePtr = NULL;
	}
}
void yousef_nucleus_seg::clearClustImagePtr(void)
{
	if(clustImagePtr)
	{
		delete[] clustImagePtr;
		clustImagePtr = NULL;
	}
}
void yousef_nucleus_seg::clearMyConnComp(void)
{
	if(myConnComp)
	{
		delete[] myConnComp;
		myConnComp = NULL;
	}
}

int yousef_nucleus_seg::saveIntoIDLFormat(std::string imageName)
{
  if(!segImagePtr && !clustImagePtr && !binImagePtr)
    return -1;

  std::string idl_name = imageName.substr(0,imageName.find_last_of("."))+"_seg_final.dat";
  std::cout<<"Save the image to "<<idl_name<<" in the IDL format ..."<<std::endl;
  //numStacks--;
  int array_size = numRows*numColumns*numStacks;
  unsigned short *memblock = new unsigned short[array_size];
  std::ofstream dat_file(idl_name.c_str(), std::ios::binary);

  // In the IDL format, the first direction is Z (fastest changing
  // index), the second is X and third is Y. The image is also flipped
  // in Y.
  int ind = 0;
  for(int y=numRows-1; y>=0; y--) {
  //for(int y=0; y<numRows; y++) {
    for(int x=0; x<numColumns; x++) {					       
      for(int z=0; z<numStacks; z++) {				
        if (segImagePtr)
          memblock[ind] = (unsigned short) segImagePtr[(z*numRows*numColumns)+(y*numColumns)+x];
		else if (clustImagePtr)
          memblock[ind] = (unsigned short) clustImagePtr[(z*numRows*numColumns)+(y*numColumns)+x];
		else
		  memblock[ind] = (unsigned short) binImagePtr[(z*numRows*numColumns)+(y*numColumns)+x];
        ind++;
      }
    }
  }

  dat_file.write(reinterpret_cast<char *>(memblock), sizeof(unsigned short)*array_size);
  dat_file.close();
  delete[] memblock;
  return 1;
}


int yousef_nucleus_seg::readFromIDLFormat(std::string fileName)
{
  if(!dataImagePtr)
    return -1;
  if(!segImagePtr)
	  segImagePtr = new unsigned short[numStacks*numRows*numColumns];

  std::string idl_name = fileName.substr(0,fileName.find_first_of("."))+"_seg_final.dat";
  std::cout<<"reading segmentation from "<<idl_name<<std::endl;

  int array_size = numRows*numColumns*numStacks;
  unsigned short *memblock = new unsigned short[array_size];
  std::ifstream dat_file(idl_name.c_str(), std::ios::binary);

  //read the file contents
  dat_file.read(reinterpret_cast<char *>(memblock), sizeof(unsigned short)*array_size);

  // In the IDL format, the first direction is Z (fastest changing
  // index), the second is X and third is Y. The image is also flipped
  // in Y.
  int ind = 0;
  for(int y=numRows-1; y>=0; y--) {
  //for(int y=0; y<numRows; y++) {
    for(int x=0; x<numColumns; x++) {					       
      for(int z=0; z<numStacks; z++) {				        
          segImagePtr[(z*numRows*numColumns)+(y*numColumns)+x] = memblock[ind];		
		  ind++;
      }
    }
  }

  dat_file.close();
  delete[] memblock;
  return 1;
}

ftk::Object::Point yousef_nucleus_seg::MergeInit(ftk::Object::Point P1, ftk::Object::Point P2, int* newID)
{
	//
	//if no label (segmentation) is available then return
	if(!clustImagePtr)
	{
		cerr<<"Run initial segmentation first"<<endl;				
		ftk::Object::Point err;
		err.t = err.x = err.y = err.z = 0;
		return err;
	}
	//Make sure the two points are inside two different cells
	int id1 = clustImagePtr[(P1.z*numRows*numColumns)+(P1.y*numColumns)+P1.x];
	int id2 = clustImagePtr[(P2.z*numRows*numColumns)+(P2.y*numColumns)+P2.x];
	if(id1 == id2)
	{
		std::cerr<<"Can't use two points inside the same cell for merging"<<std::endl;
		ftk::Object::Point err;
		err.t = err.x = err.y = err.z = 0;
		return err;
	}

	
	//Also, we have to make sure that we are trying to merge two adjacent nuclei
	//First, get the bounding boxes
	std::vector< ftk::Object::Point > BOX1 = getObjectBoundingBox(id1, 1);
	ftk::Object::Point bBox1 = BOX1.at(0);
	ftk::Object::Point bBox2 = BOX1.at(1);
	std::vector< ftk::Object::Point > BOX2 = getObjectBoundingBox(id2, 1);
	ftk::Object::Point bBox3 = BOX2.at(0);
	ftk::Object::Point bBox4 = BOX2.at(1);
	
	//then, go over all the points of the first cell and see if any of its neighbor points belongs to the other cell
	int adj = 0;
	if(bBox1.y == 0)
		bBox1.y = 1;
	if(bBox2.y == numRows-1)
		bBox2.y = numRows-2;
	if(bBox1.x == 0)
		bBox1.x = 1;
	if(bBox2.x == numColumns-1)
		bBox2.x = numColumns-2;
	if(bBox1.z == 0)
		bBox1.z = 1;
	if(bBox2.z == numStacks-1)
		bBox2.z = numStacks-2;
	for( int i=bBox1.y; i<=bBox2.y; i++)
	{
		for( int j=bBox1.x; j<=bBox2.x; j++)
		{
			for( int k=bBox1.z; k<=bBox2.z; k++)
			{
				int pix = clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j];
				if(pix != id1)
					continue;
				for( int si=i-1; si<=i+1; si++)
				{
					for( int sj=j-1; sj<=j+1; sj++)
					{
						for( int sk=k-1; sk<=k+1; sk++)
						{
							if(clustImagePtr[(sk*numRows*numColumns)+(si*numColumns)+sj]==id2)
							{
								adj = 1;
								break;
							}
						}
					}
				}
				
			}
		}
	}

	if(adj == 0)
	{
		std::cerr<<"The two cells need to be adjacent in order to merge them"<<std::endl;
	}
	//Update the initial segmentation image	
	int min_x = min(bBox1.x, bBox3.x);
	int min_y = min(bBox1.y, bBox3.y);
	int min_z = min(bBox1.z, bBox3.z);
	int max_x = max(bBox2.x, bBox4.x);
	int max_y = max(bBox2.y, bBox4.y);
	int max_z = max(bBox2.z, bBox4.z);

	std::vector <int> sz;
	sz.push_back(max_x - min_x + 1);
	sz.push_back(max_y - min_y + 1);
	sz.push_back(max_z - min_z + 1);

	//these will be needed to get the center of the new cell
	double med_x = 0.0;
	double med_y = 0.0;
	double med_z = 0.0;
	int id1_count = 0;
	int id2_count = 0;
	for( int i=min_y; i<=max_y; i++)
	{
		for( int j=min_x; j<=max_x; j++)
		{
			for( int k=min_z; k<=max_z; k++)
			{
				int pix = clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j];
				if(pix == id1)
				{
					id1_count++;
					med_x +=j;					
					med_y +=i;					
					med_z +=k;
				}
				else if(pix == id2)
				{
					id2_count++;
					med_x +=j;					
					med_y +=i;					
					med_z +=k;
					clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j] = id1;
				}
			}
		}
	}
	med_x /=(id1_count+id2_count);
	med_y /=(id1_count+id2_count);
	med_z /=(id1_count+id2_count);

	

	//get the indexes of the center with respect to the begining of the bounding box of the new cell
	int ind1 = ((med_z- min_z)*sz[0]*sz[1]) + ((med_y - min_y)*sz[0]) + (med_x - min_x);	

	
	//create a new itk images with the same size as the bounding box	 	 
	typedef    float     InputPixelType;
	typedef itk::Image< InputPixelType,  3 >   InputImageType;
	InputImageType::Pointer sub_im1;
	sub_im1 = InputImageType::New();	
	InputImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	origin[2] = 0.;    
    sub_im1->SetOrigin( origin );

    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    InputImageType::SizeType  size;
    size[0]  = sz[0];  // size along X
    size[1]  = sz[1];  // size along Y
	size[2]  = sz[2];  // size along Z
  
    InputImageType::RegionType rgn;
    rgn.SetSize( size );
    rgn.SetIndex( start );
    
    sub_im1->SetRegions( rgn );
    sub_im1->Allocate();
    sub_im1->FillBuffer(0);
	sub_im1->Update();	
	
	//set all the points in this image to zero except for the center point
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(sub_im1,sub_im1->GetRequestedRegion());
			
	for(int i=0; i<sz[0]*sz[1]*sz[2]; i++)
	{				
		if(i==ind1)
			iterator1.Set(255.0);
		else
			iterator1.Set(0.0);		
				
		++iterator1;			
	}
	

	//compute the distance transforms of that binary itk image
	typedef    float     OutputPixelType;
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType;
	typedef itk::DanielssonDistanceMapImageFilter<InputImageType, OutputImageType > DTFilter ;
	DTFilter::Pointer dt_obj1= DTFilter::New() ;
	dt_obj1->UseImageSpacingOn();
	dt_obj1->SetInput(sub_im1) ;
	
	try{
		dt_obj1->Update() ;
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "Error calculating distance transform: " << err << endl ;
	}
	
	//get the max distance so that you can invert the distance map
	IteratorType iterator3(dt_obj1->GetOutput(),dt_obj1->GetOutput()->GetRequestedRegion());
	int mx1 = 0;
	for(int k=min_z; k<=max_z; k++)
	{
		for(int i=min_y; i<=max_y; i++)
		{			
			for(int j=min_x; j<=max_x; j++)
			{
				int d1 = (int) fabs(iterator3.Get());	 

				if(d1>mx1)
					mx1 = d1;			
				++iterator3;
			}
		}
	}
	iterator3.GoToBegin();
	for(int k=min_z; k<=max_z; k++)
	{
		for(int i=min_y; i<=max_y; i++)
		{			
			for(int j=min_x; j<=max_x; j++)
			{
				int d1 = (int) fabs(iterator3.Get());	 					
				++iterator3;
				
				int pix = clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j];

				if(pix == id1)
					logImagePtr[(k*numRows*numColumns)+(i*numColumns)+j]= mx1-d1;				
			}
		}
	}	
	
	//send back the id of the new cell after merging
	newID[0] = id1;

	//return the coordinates of the seed (center) of the new cell
	ftk::Object::Point new_seed;
	new_seed.t = 0;
	new_seed.x = med_x;
	new_seed.y = med_y;
	new_seed.z = med_z;
	
	return new_seed;
}

vector< int > yousef_nucleus_seg::SplitInit(ftk::Object::Point P1, ftk::Object::Point P2)
{
	//
	//if no label (segmentation) or no data image is available then return
	if(!clustImagePtr)
	{
		cerr<<"Run initial segmentation first"<<endl;
		std::vector <int> ids_err;
		ids_err.push_back(0);
		ids_err.push_back(0);
		return ids_err;
	}
	//Check if the two points inside the same cell
	int id1 = clustImagePtr[(P1.z*numRows*numColumns)+(P1.y*numColumns)+P1.x];
	int id2 = clustImagePtr[(P2.z*numRows*numColumns)+(P2.y*numColumns)+P2.x];
	if(id1 != id2)
	{
		std::vector <int> ids_err;
		ids_err.push_back(0);
		ids_err.push_back(0);
		return ids_err;
	}

	//Update the initial segmentation image	
	//get the indexes of the two seeds with respect to the begining of the bounding box	
	std::vector< ftk::Object::Point > bBox = getObjectBoundingBox(id1, 1);
	ftk::Object::Point bBox1 = bBox.at(0);
	ftk::Object::Point bBox2 = bBox.at(1);
	std::vector <int> sz;
	sz.push_back(bBox2.x - bBox1.x + 1);
	sz.push_back(bBox2.y - bBox1.y + 1);
	sz.push_back(bBox2.z - bBox1.z + 1);
	
	int ind1 = ((P1.z- bBox1.z)*sz[0]*sz[1]) + ((P1.y - bBox1.y)*sz[0]) + (P1.x - bBox1.x);
	int ind2 = ((P2.z- bBox1.z)*sz[0]*sz[1]) + ((P2.y - bBox1.y)*sz[0]) + (P2.x - bBox1.x);

	//create two new itk images with the same size as the bounding box	 	 
	typedef    float     InputPixelType;
	typedef itk::Image< InputPixelType,  3 >   InputImageType;
	InputImageType::Pointer sub_im1;
	sub_im1 = InputImageType::New();
	InputImageType::Pointer sub_im2;
	sub_im2 = InputImageType::New();
	InputImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	origin[2] = 0.;    
    sub_im1->SetOrigin( origin );
	sub_im2->SetOrigin( origin );

    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    InputImageType::SizeType  size;
    size[0]  = sz[0];  // size along X
    size[1]  = sz[1];  // size along Y
	size[2]  = sz[2];  // size along Z
  
    InputImageType::RegionType rgn;
    rgn.SetSize( size );
    rgn.SetIndex( start );
    
   
    sub_im1->SetRegions( rgn );
    sub_im1->Allocate();
    sub_im1->FillBuffer(0);
	sub_im1->Update();	
	sub_im2->SetRegions( rgn );
    sub_im2->Allocate();
    sub_im2->FillBuffer(0);
	sub_im2->Update();

	//set all the points in those images to zeros except for the two points corresponding to the two new seeds
	//notice that one seed is set in each image
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(sub_im1,sub_im1->GetRequestedRegion());
	IteratorType iterator2(sub_im2,sub_im2->GetRequestedRegion());	
		
	for(int i=0; i<sz[0]*sz[1]*sz[2]; i++)
	{				
		if(i==ind1)
			iterator1.Set(255.0);
		else
			iterator1.Set(0.0);		
		if(i==ind2)
			iterator2.Set(255.0);
		else
		iterator2.Set(0.0);
		
		++iterator1;	
		++iterator2;	
	}
	

	//compute the distance transforms of those binary itk images
	typedef    float     OutputPixelType;
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType;
	typedef itk::DanielssonDistanceMapImageFilter<InputImageType, OutputImageType > DTFilter ;
	DTFilter::Pointer dt_obj1= DTFilter::New() ;
	DTFilter::Pointer dt_obj2= DTFilter::New() ;
	dt_obj1->UseImageSpacingOn();
	dt_obj1->SetInput(sub_im1) ;
	dt_obj2->UseImageSpacingOn();
	dt_obj2->SetInput(sub_im2) ;

	try{
		dt_obj1->Update() ;
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "Error calculating distance transform: " << err << endl ;
	}
	try{
		dt_obj2->Update() ;
	}
	catch( itk::ExceptionObject & err ){
		std::cerr << "Error calculating distance transform: " << err << endl ;
	}

	//Now, relabel the cell points into either the old id (id1) or a new id (newID) based on the distances to the seeds
	if(numObjects == 0)
		numObjects = getMaxID(1);

	++numObjects;		
	int newID1 = numObjects;
	++numObjects;		
	int newID2 = numObjects;
	IteratorType iterator3(dt_obj1->GetOutput(),dt_obj1->GetOutput()->GetRequestedRegion());
	IteratorType iterator4(dt_obj2->GetOutput(),dt_obj2->GetOutput()->GetRequestedRegion());	
	int mx1 = 0;
	int mx2 = 0;
	for(int k=bBox1.z; k<=bBox2.z; k++)
	{
		for(int i=bBox1.y; i<=bBox2.y; i++)
		{			
			for(int j=bBox1.x; j<=bBox2.x; j++)
			{
				int d1 = (int) fabs(iterator3.Get());	 
				int d2 = (int) fabs(iterator4.Get());	
				if(d1>mx1)
					mx1 = d1;
				if(d2>mx2)
					mx2 = d2;
				++iterator3;
				++iterator4;
			}
		}
	}
	iterator3.GoToBegin();
	iterator4.GoToBegin();
	for(int k=bBox1.z; k<=bBox2.z; k++)
	{
		for(int i=bBox1.y; i<=bBox2.y; i++)
		{			
			for(int j=bBox1.x; j<=bBox2.x; j++)
			{
				int d1 = (int) fabs(iterator3.Get());	 
				int d2 = (int) fabs(iterator4.Get());					
				++iterator3;
				++iterator4;				
				int pix = clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j];
				//if(pix != id1)
				//	continue; //maybe we need to update the LoG resp if pix = 0
				if(d1>d2)	
				{
					if(pix == id1)
						clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j] = newID1;
					//also update the LoG resp Image
					logImagePtr[(k*numRows*numColumns)+(i*numColumns)+j]= mx2-d2;
				}
				else
				{
					if(pix == id1)
						clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j] = newID2;
					//also update the LoG resp Image
					logImagePtr[(k*numRows*numColumns)+(i*numColumns)+j]= mx1-d1;
				}

			}
		}
	}	
	

	//return the ids of the two cells resulting from spliting
	std::vector <int> ids_ok;
	ids_ok.push_back(newID1);
	ids_ok.push_back(newID2);

	//also, add the old ID to the end of the list
	ids_ok.push_back(id1);

	return ids_ok;
}

//this is applied on the initial segmentatio, and it also updates the binary image
bool yousef_nucleus_seg::DeleteInit(ftk::Object::Point P1)
{
	if(!clustImagePtr)
	{
		cerr<<"Run initial segmentation first"<<endl;		
		return false;
	}
	//Check if the two points inside the same cell
	int id = clustImagePtr[(P1.z*numRows*numColumns)+(P1.y*numColumns)+P1.x];
	for (int i=0; i<numRows; i++)
	{
		for(int j=0; j<numColumns; j++)
		{
			for( int k=0; k<numStacks; k++)
			{
				int pix = clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j];
				if(pix == id)
				{
					clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j] = 0;
					binImagePtr[(k*numRows*numColumns)+(i*numColumns)+j] = 0;
				}
			}
		}
	}
	return true;
}

std::vector< ftk::Object::Point > yousef_nucleus_seg::getObjectBoundingBox(int id, int Int_Fin)
{
	int min_x, max_x, min_y, max_y, min_z, max_z;
	max_x = max_y = max_z = 0;
	min_x = numColumns;
	min_y = numRows;
	min_z = numStacks;
	
	unsigned short* imgPtr;
	if(Int_Fin == 1)
		imgPtr = clustImagePtr;
	else
		imgPtr = segImagePtr;

	int counter = 0;
	for( int i=0; i<numRows; i++)
	{
		for( int j=0; j<numColumns; j++)
		{
			for( int k=0; k<numStacks; k++)
			{
				int ID = imgPtr[(k*numRows*numColumns)+(i*numColumns)+j];
				if(ID!=id)
					continue;

				counter++;

				if(i<min_y)
					min_y = i;
				if(i>max_y)
					max_y = i;
				if(j<min_x)
					min_x = j;
				if(j>max_x)
					max_x = j;
				if(k<min_z)
					min_z = k;
				if(k>max_z)
					max_z = k;
			}
		}
	}
	std::vector< ftk::Object::Point > bBox;
	ftk::Object::Point P1;
	ftk::Object::Point P2;
	P1.t = 0;
	if(counter == 0)
		P1.x = P1.y = P1.z = P2.x = P2.y = P2.z = 0;
	else
	{
		P1.x = min_x;
		P1.y = min_y;
		P1.z = min_z;
		P2.t = 0;
		P2.x = max_x;
		P2.y = max_y;
		P2.z = max_z;
	}
	bBox.push_back(P1);
	bBox.push_back(P2);

	return bBox;
}

int yousef_nucleus_seg::getMaxID(int Int_Fin)
{
	int maxID = 0;
	int iid = 0;
	
		
	for( int i=0; i<numRows; i++)
	{
		for( int j=0; j<numColumns; j++)
		{
			for( int k=0; k<numStacks; k++)
			{				
				if(Int_Fin == 1)
					iid = clustImagePtr[(k*numRows*numColumns)+(i*numColumns)+j];
				else
					iid = segImagePtr[(k*numRows*numColumns)+(i*numColumns)+j];
				if(iid>maxID)
					maxID = iid;
			}
		}
	}
	return iid;
}

//int yousef_nucleus_seg::AddObject(ftk::Object::Point P1, ftk::Object::Point P2)
int yousef_nucleus_seg::AddObject(unsigned char* inImage, unsigned short* lbImage, std::vector<int> P1, std::vector<int> P2,
					std::vector<itk::SizeValueType> imSZ, int maxID)
{		
	//get the coordinates of the two points and the size of the box
	int x1 = P1[0];
	int x2 = P2[0];
	int y1 = P1[1];
	int y2 = P2[1];
	int z1 = P1[2];
	int z2 = P2[2];

	int sz_x = x2-x1+1;
	int sz_y = y2-y1+1;
	//if(z1==z2)
	//{
	//	//assume that the spacing ratio is 2
	//	int dz = sz_x/2;
	//	z1 -= dz;
	//	z2 += dz;
	//}
	
	int sz_z = z2-z1+1;
	/*if(sz_x<1 || sz_y<1 || sz_z<1)
		return 0;*/
	
	int cent_x = (x2+x1)/2;
	int cent_y = (y2+y1)/2;
	int cent_z = (z2+z1)/2;
	int max_dist = (int) sqrt((double)(x2-cent_x)*(x2-cent_x)+(y2-cent_y)*(y2-cent_y)+(z2-cent_z)*(z2-cent_z));
	//extract the region from the raw image	
	unsigned char* subDataImagePtr = new unsigned char[sz_x*sz_y*sz_z];
	int ind =0;	
	for(int k=z1; k<=z2; k++)
	{
		for(int i=y1; i<=y2; i++)
		{
			for(int j=x1; j<=x2; j++)
			{
				//try this: enhance the intensity at each point by multiplying by a weight
				//that is inversly proportional to its distance from the center of the box
				int d = (int) sqrt((double)(j-cent_x)*(j-cent_x)+(i-cent_y)*(i-cent_y)+(k-cent_z)*(k-cent_z));
				d = max_dist-d;				
				d/=2;

				subDataImagePtr[ind] = inImage[(k*imSZ[2]*imSZ[3])+i*imSZ[3]+j]*d;
				ind++;
			}
		}
	}

	//create an impty image to hold the binarization results
	unsigned short* subBinImagePtr = new unsigned short[sz_x*sz_y*sz_z];
	memset(subBinImagePtr/*destination*/,0/*value*/,sz_x*sz_y*sz_z*sizeof(unsigned short)/*num bytes to move*/);

	int ok = 0;
	//if (sz_z == 1)
	//{
	//	ok = Cell_Binarization_2D(subDataImagePtr,subBinImagePtr, sz_y, sz_x, 0); //Do Binarization		
	//}
	//else
	//{
		ok = Cell_Binarization_3D(subDataImagePtr,subBinImagePtr, sz_y, sz_x, sz_z, 0, 0, 128);		//Do Binarization		
	//}
	if(ok==0)
		return 0;
	

	//Added by Yousef on 9/26/2009
	//Apply binary morphological closing on the resulting binary object
	//First, create an ITK image to hold the binary sub-image that contain the new object
	typedef    unsigned short     TempPixelType;
	typedef itk::Image< TempPixelType,  3 >   TempImageType;
	TempImageType::Pointer sub_im;
	sub_im = TempImageType::New();	
	
	TempImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	origin[2] = 0.;    
    sub_im->SetOrigin( origin );

    TempImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    
	TempImageType::SizeType  size;
    size[0]  = sz_x;  // size along X
    size[1]  = sz_y;  // size along Y
	size[2]  = sz_z;  // size along Z
  
    TempImageType::RegionType rgn;
    rgn.SetSize( size );
    rgn.SetIndex( start );
    
   
    sub_im->SetRegions( rgn );
    sub_im->Allocate();
    sub_im->FillBuffer(0);
	sub_im->Update();	
	

	//Copy the sub-binary image into the ITK image
	//notice that one seed is set in each image
	typedef itk::ImageRegionIteratorWithIndex< TempImageType > IteratorType;
	IteratorType iterator(sub_im,sub_im->GetRequestedRegion());			
	for(int i=0; i<sz_x*sz_y*sz_z; i++)
	{						
		iterator.Set(subBinImagePtr[i]);				
		++iterator;	
	}

	//Binary Morphological closing....
	//1-Create the structuring element
	typedef itk::BinaryBallStructuringElement< TempPixelType, 3 > KernelType;
	KernelType ball;
	KernelType::SizeType ballSize;
	ballSize.Fill( 4 ); //for now, set the radius to 4
	ball.SetRadius( ballSize );
	ball.CreateStructuringElement();
	
	typedef itk::BinaryMorphologicalClosingImageFilter< TempImageType, TempImageType, KernelType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( sub_im );
	filter->SetKernel( ball );
	// test the default attribute values, and exercise the accesors
	if( !filter->GetSafeBorder() )
	{
		std::cerr << "Wrong SafeBorder default value" << std::endl;
		return EXIT_FAILURE;
    }
	filter->SafeBorderOff();
	filter->SafeBorderOn();
	filter->SetSafeBorder( true );

	/*if( filter->GetForegroundValue() != 255 )
    {
    std::cerr << "Wrong Foreground default value" << std::endl;
    return EXIT_FAILURE;
    }*/
	filter->SetForegroundValue( 255 ); 

    try
    {
		filter->Update();
    } 
    catch ( itk::ExceptionObject & excp )
    {
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
    }
		
	IteratorType iterator2(filter->GetOutput(),filter->GetOutput()->GetRequestedRegion());				

	//typedef itk::ImageFileWriter< TempImageType > WriterType;
	//WriterType::Pointer writer = WriterType::New();
	//writer->SetInput( filter->GetOutput() );
	//writer->SetFileName( "anothertest.tif" );
	//writer->Update();
	////
	maxID++;

	//ind = -1;
	for(int k=z1; k<=z2; k++)
	{
		for(int i=y1; i<=y2; i++)
		{
			for(int j=x1; j<=x2; j++)
			{												
				//ind++;
				//if(subBinImagePtr[ind] == 0 || lbImage[(k*imSZ[2]*imSZ[3])+i*imSZ[3]+j]!=0)
				//	continue;

				//This following code replaces the three commented lines above
				unsigned short tmp = iterator2.Get();
				++iterator2;
				if(tmp == 0 || lbImage[(k*imSZ[2]*imSZ[3])+i*imSZ[3]+j]!=0)
					continue;				
				////////////////////////////////////////////////////////////

				if(k==z1 || k==z2 || i==y1 || i==y2 || j==x1 || j==x2) //remove border points
					continue;

				lbImage[(k*imSZ[2]*imSZ[3])+i*imSZ[3]+j] = maxID;				
			}
		}
	}
	//free memory
	delete [] subDataImagePtr;
	delete [] subBinImagePtr;

	//return the new object ID
	return maxID;					
}


int yousef_nucleus_seg::AddObject2D(unsigned char* inImage, unsigned short* lbImage, std::vector<int> P1, std::vector<int> P2,
					std::vector<itk::SizeValueType> imSZ, int maxID)
{		
	//get the coordinates of the two points and the size of the box
	int x1 = P1[0];
	int x2 = P2[0];
	int y1 = P1[1];
	int y2 = P2[1];

	int sz_x = x2-x1+1;
	int sz_y = y2-y1+1;		
	
	int cent_x = (x2+x1)/2;
	int cent_y = (y2+y1)/2;	
	int max_dist = (int) sqrt((double)(x2-cent_x)*(x2-cent_x)+(y2-cent_y)*(y2-cent_y));
	//extract the region from the raw image	
	unsigned char* subDataImagePtr = new unsigned char[sz_x*sz_y];
	int ind =0;		
	for(int i=y1; i<=y2; i++)
	{
		for(int j=x1; j<=x2; j++)
		{
			//try this: enhance the intensity at each point by multiplying by a weight
			//that is inversly proportional to its distance from the center of the box
			int d = (int) sqrt((double)(j-cent_x)*(j-cent_x)+(i-cent_y)*(i-cent_y));
			d = max_dist-d;				
			d/=2;

			subDataImagePtr[ind] = inImage[i*imSZ[3]+j]*d;
			ind++;
		}
	}
	
	//create an impty image to hold the binarization results
	unsigned short* subBinImagePtr = new unsigned short[sz_x*sz_y];
	memset(subBinImagePtr/*destination*/,0/*value*/,sz_x*sz_y*sizeof(unsigned short)/*num bytes to move*/);
	
	int ok = Cell_Binarization_2D(subDataImagePtr,subBinImagePtr, sz_y, sz_x, 0); //Do Binarization		
	if(ok==0)
		return 0;
	

	//Added by Yousef on 9/26/2009
	//Apply binary morphological closing on the resulting binary object
	//First, create an ITK image to hold the binary sub-image that contain the new object
	typedef    unsigned short     TempPixelType;
	typedef itk::Image< TempPixelType,  2 >   TempImageType;
	TempImageType::Pointer sub_im;
	sub_im = TempImageType::New();	
	
	TempImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	
    sub_im->SetOrigin( origin );

    TempImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	    
	TempImageType::SizeType  size;
    size[0]  = sz_x;  // size along X
    size[1]  = sz_y;  // size along Y	
  
    TempImageType::RegionType rgn;
    rgn.SetSize( size );
    rgn.SetIndex( start );
    
   
    sub_im->SetRegions( rgn );
    sub_im->Allocate();
    sub_im->FillBuffer(0);
	sub_im->Update();	
	

	//Copy the sub-binary image into the ITK image
	//notice that one seed is set in each image
	typedef itk::ImageRegionIteratorWithIndex< TempImageType > IteratorType;
	IteratorType iterator(sub_im,sub_im->GetRequestedRegion());			
	for(int i=0; i<sz_x*sz_y; i++)
	{						
		iterator.Set(subBinImagePtr[i]);				
		++iterator;	
	}

	//Binary Morphological closing....
	//1-Create the structuring element
	typedef itk::BinaryBallStructuringElement< TempPixelType, 2 > KernelType;
	KernelType ball;
	KernelType::SizeType ballSize;
	ballSize.Fill( 4 ); //for now, set the radius to 4
	ball.SetRadius( ballSize );
	ball.CreateStructuringElement();
	
	typedef itk::BinaryMorphologicalClosingImageFilter< TempImageType, TempImageType, KernelType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( sub_im );
	filter->SetKernel( ball );
	// test the default attribute values, and exercise the accesors
	if( !filter->GetSafeBorder() )
	{
		std::cerr << "Wrong SafeBorder default value" << std::endl;
		return EXIT_FAILURE;
    }
	filter->SafeBorderOff();
	filter->SafeBorderOn();
	filter->SetSafeBorder( true );
	
	filter->SetForegroundValue( 255 ); 

    try
    {
		filter->Update();
    } 
    catch ( itk::ExceptionObject & excp )
    {
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
    }
		
	IteratorType iterator2(filter->GetOutput(),filter->GetOutput()->GetRequestedRegion());				

	//typedef itk::ImageFileWriter< TempImageType > WriterType;
	//WriterType::Pointer writer = WriterType::New();
	//writer->SetInput( filter->GetOutput() );
	//writer->SetFileName( "anothertest.tif" );
	//writer->Update();
	////
	maxID++;

	//ind = -1;
	for(int i=y1; i<=y2; i++)
	{
		for(int j=x1; j<=x2; j++)
		{												
			//ind++;
			//if(subBinImagePtr[ind] == 0 || lbImage[(k*imSZ[2]*imSZ[3])+i*imSZ[3]+j]!=0)
			//	continue;

			//This following code replaces the three commented lines above
			unsigned short tmp = iterator2.Get();
			++iterator2;
			if(tmp == 0 || lbImage[i*imSZ[3]+j]!=0)
				continue;				
			////////////////////////////////////////////////////////////

			if(i==y1 || i==y2 || j==x1 || j==x2) //remove border points
				continue;

			lbImage[i*imSZ[3]+j] = maxID;					
		}
	}
	//free memory
	delete [] subDataImagePtr;
	delete [] subBinImagePtr;

	//return the new object ID
	return maxID;					
}

/*
int yousef_nucleus_seg::saveIntoIDLFormat()
{
	if(!segImagePtr || !clustImagePtr)
		return -1;

	//added by Yousef on 8-4-2008: Write the segmentation output into a dat file readable by the IDL farsight
	unsigned short *IDL_seg = new unsigned short[numRows*numColumns*numStacks];
  #ifdef WIN32
	  FILE* fid = fopen("c:\w1_seg_final.dat","w");
  #else
	  FILE* fid = fopen("/tmp/w1_seg_final.dat","w");
  //we should probably also have a case for MacOs...
  #endif
	int ind = -1;
	int mx = 1;
	for(int k=0; k<numStacks; k++)										
	{
		for(int i=0; i<numColumns; i++)										
		{						
			for(int j=0; j<numRows; j++)
			{				
				ind++;		
				if(segImagePtr)
					IDL_seg[ind] = (unsigned short) segImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];							else
					IDL_seg[ind] = (unsigned short) clustImagePtr[(k*numRows*numColumns)+(j*numColumns)+i];			
				//fputc(IDL_seg[ind],fid);
			}
		}
	}
	int siz = sizeof(IDL_seg[0]);
	fwrite(IDL_seg,siz,numRows*numColumns*numStacks,fid);
-	fclose(fid);

	return 1;
}
*/

