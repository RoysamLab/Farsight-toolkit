/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#include "yousef_seg.h"
#include <fstream>

//Constructor
yousef_nucleus_seg::yousef_nucleus_seg()
{
	dataImagePtr = NULL;
	binImagePtr = NULL;
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

void yousef_nucleus_seg::setDataImage( unsigned char *imgPtr,  int x, int y, int z, const char *filename )
{
	numStacks = z;
	numRows = y;//y();			//y-direction
	numColumns = x; 		//x-direction

	dataFilename = filename;

	//delete[] dataImagePtr; I will never delete the dataImagePtr because it is created outside this class
	dataImagePtr = imgPtr;
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
void yousef_nucleus_seg::runBinarization()
{
	//First check to be sure that we have a dataImage to use
	if (!dataImagePtr)
		return;
	
	//Now clear all subsequent variables (dependent upon this binary image
	numConnComp = 0;
	clearBinImagePtr();
	clearSeedImagePtr();
	clearLogImagePtr();
	clearSegImagePtr();
	clearClustImagePtr();
	clearMyConnComp();
	mySeeds.clear();

	//allocate space for the binary image
	binImagePtr = new int[numStacks*numRows*numColumns];

	int ok = 0;
	if (numStacks == 0)
	{}
	else if (numStacks == 1)
	{
		ok = Cell_Binarization_2D(dataImagePtr,binImagePtr, numRows, numColumns, shift);
		numConnComp = getConnCompImage(binImagePtr, 26, minObjSize, numRows, numColumns, numStacks,1);
		getConnCompInfo3D();
	}
	else
	{
		ok = Cell_Binarization_3D(dataImagePtr,binImagePtr, numRows, numColumns, numStacks, shift);		//Do Binarization
		numConnComp = getConnCompImage(binImagePtr, 26, minObjSize, numRows, numColumns, numStacks,1);			//Find connected components
		getConnCompInfo3D();																			//Populate myConnComp
	}

	if(ok)
	{
		cerr << "Cell Binarized.. with " << numConnComp << " connected components" << endl;	
	}
	else
	{
		cerr << "Binarization Failed!!" << endl;
	}
}

void yousef_nucleus_seg::runSeedDetection()
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

	//allocate space for the binary image of seed points
	seedImagePtr = new int[numStacks*numRows*numColumns];
	
	//copy the binary image into the seeds image for now
	//memcpy(seedImagePtr/*destination*/, binImagePtr/*source*/, numStacks*numRows*numColumns*sizeof(int)/*num bytes to move*/);

	//need to pass a float pointer with input image in it, so create it here
	float *imgPtr = new float[numStacks*numRows*numColumns];
	ucharToFloat(dataImagePtr /*from*/, imgPtr /*to*/, numRows, numColumns, numStacks, 1 /*invert*/);

	//allocate space for the laplacian of gaussian
	logImagePtr = new float[numStacks*numRows*numColumns];
	
	//Now do seed detection
	int ok = 0;
	if (numStacks == 1)
	{
		ok = detectSeeds2D( imgPtr, logImagePtr, seedImagePtr, numRows, numColumns, scaleMin, scaleMax, regionXY, binImagePtr );
	}
	else
	{
		minLoGImg = Seeds_Detection_3D( imgPtr, logImagePtr, seedImagePtr, numRows, numColumns, numStacks, scaleMin, scaleMax, regionXY, regionZ, getSamplingRatio(), binImagePtr, useDistMap );
		ok = 1;
	}
	std::cout << "Seeds Detected? " << ok << std::endl;

	delete [] imgPtr;	//cleanup

	//Make sure all seeds are in foreground and extract vector of seeds
	ExtractSeeds();

}

void yousef_nucleus_seg::runClustering()
{
	//Check for required images
	if( !dataImagePtr || !logImagePtr || !seedImagePtr || !binImagePtr )
		return;

	//Now clear all subsequent variables				
	clearSegImagePtr();
	clearClustImagePtr();

	//Allocate space
	clustImagePtr = new int[numStacks*numRows*numColumns];

	if (numStacks == 1)
	{
		std::cout << "Starting Initial Clustering" << std::endl;
		//ExtractSeeds();
		int *seed_xmclust, *seed_ymclust;
		int numseedsmclust = mySeeds.size();
		seed_xmclust = (int *) malloc(mySeeds.size()*sizeof(int));
		seed_ymclust = (int *) malloc(mySeeds.size()*sizeof(int));
		for (int i=0; i<((int)mySeeds.size()); ++i)
		{
			seed_ymclust[i] = mySeeds[i].y();
			seed_xmclust[i] = mySeeds[i].x();
		}
		local_max_clust_2D(logImagePtr, numRows, numColumns, 5.0, clustImagePtr, seed_xmclust, seed_ymclust, numseedsmclust, binImagePtr);
		free( seed_xmclust );
		free( seed_ymclust );
	}
	else
	{
		std::cerr << "Starting Initial Clustering" << std::endl;
		local_max_clust_3D(logImagePtr/*LoG*/, seedImagePtr/*local max vals*/, binImagePtr/*binary mask*/,clustImagePtr/*output*/,\
			numRows, numColumns, numStacks, regionXY, regionZ);		
	}
}

void yousef_nucleus_seg::ExtractSeeds()
{
	mySeeds.clear();

	int seedVal;
	int binVal;
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
						seedImagePtr[curNode] = -1;
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
	std::cerr << id-1 << " seeds were detected"<<std::endl;
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


int yousef_nucleus_seg::getConnCompImage(int *IM, int connectivity, int minSize, int r, int c, int z, int runConnComp)
{
	typedef    int     InputPixelType;
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
				
	typedef itk::ConnectedComponentImageFilter< InputImageType, OutputImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelType;
	RelabelType::Pointer relabel = RelabelType::New();

	if(runConnComp == 1)
	{
		//Compute the labeled connected component image		
		filter->SetInput (im);
		filter->SetFullyConnected( connectivity );		
		//use the connected component image as the input to the relabel component filter		
		relabel->SetInput( filter->GetOutput() );
	}
	else
	{
		//use the input image as the input to the relabel component filter 		
		relabel->SetInput( im );
	}

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
    }

	
    //write the output of the labeling CC filter into our input image
	IteratorType iterator2(relabel->GetOutput(),relabel->GetOutput()->GetRequestedRegion());
	for(int i=0; i<r*c*z; i++)
	{		
		IM[i] = iterator2.Get();		
		++iterator2;	
	}

	//return the number of CCs
	return relabel->GetNumberOfObjects();;
}

//By Yousef: 5-21-2008
//Get the connected component information of the binary (conncomp) image
void yousef_nucleus_seg::getConnCompInfo3D()
{
	int val, curNode;
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

void yousef_nucleus_seg::runAlphaExpansion2D(){
	//First check for necessary prerequisites
	if( !dataImagePtr || !logImagePtr || !seedImagePtr || !binImagePtr || !myConnComp ){
		return;
	}

	//Now clear all subsequent variables
	clearSegImagePtr();

	std::cerr<<"Finalizing Segmentation"<<std::endl;

	//Now, we apply the next steps into the connected components one by one
	int ind, x_len, y_len, val;

	segImagePtr = new int[numRows*numColumns];
	memset(segImagePtr/*destination*/,0/*value*/,numStacks*numRows*numColumns*sizeof(int)/*num bytes to move*/);

	for( int n=0; n<numConnComp; n++ ){
		std::cerr<<"Processing Connected Component #"<<n+1<<"...";
		//Now, get the subimages (the bounding box) for the current connected component
		ind = 0;
		x_len = myConnComp[n].x2 - myConnComp[n].x1 + 1;
		y_len = myConnComp[n].y2 - myConnComp[n].y1 + 1;
		float* sublogImg = new float[x_len*y_len];
		int* subclustImg = new int[x_len*y_len];	
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
			std::cerr<<"Done with only one object"<<std::endl;
			delete[] sublogImg;
			delete[] subclustImg;
			continue;
		}
		std::cerr<<std::endl<<"    "<<labelsList.size()<<" objects found"<<std::endl;
		//If you reach here, it means that the current connected component contains two or more cells
		//First, sort the labels list
		std::cerr<<"    "<<"sorting labels"<<std::endl;
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
		std::cerr<<"Done"<<std::endl;
		delete [] sublogImg;
		delete [] subclustImg;
		//delete [] subDataImg;
	}
	//relabel the cells
	int numOfObjs = /*getConnCompImage*/getRelabeledImage(segImagePtr, 8, 25, numRows, numColumns,numStacks, 1);	
    numOfObjs--;
	std::cerr << "done with " << numOfObjs<<" found"<<std::endl;
	std::cerr << "Creating Final Label Image" << std::endl;	
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

	std::cerr<<"Finalizing Segmentation"<<std::endl;

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
	segImagePtr = new int[numStacks*numRows*numColumns];
	//memcpy(segImagePtr/*destination*/, clustImagePtr/*source*/, numStacks*numRows*numColumns*sizeof(int)/*num bytes to move*/);
	memset(segImagePtr/*destination*/,0/*value*/,numStacks*numRows*numColumns*sizeof(int)/*num bytes to move*/);

	for(int n=0; n<numConnComp; n++)
	{
		std::cerr<<"Processing Connected Component #"<<n+1<<"...";
		//Now, get the subimages (the bounding box) for the current connected component
		ind = 0;
		x_len = myConnComp[n].x2 - myConnComp[n].x1 + 1;
		y_len = myConnComp[n].y2 - myConnComp[n].y1 + 1;
		z_len = myConnComp[n].z2 - myConnComp[n].z1 + 1;
		float* sublogImg = new float[x_len*y_len*z_len];
		int* subclustImg = new int[x_len*y_len*z_len];	
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
			std::cerr<<"Done with only one object"<<std::endl;
			delete[] sublogImg;
			delete[] subclustImg;
			continue;
		}
		
		std::cerr<<std::endl<<"    "<<labelsList.size()<<" objects found"<<std::endl;
		//If you reach here, it means that the current connected component contains two or more cells
		//First, sort the labels list
		std::cerr<<"    "<<"sorting labels"<<std::endl;
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
		int* subsegImg = new int[x_len*y_len*z_len];			
		float* Dterms  = multiColGraphLearning(sublogImg, subclustImg, subsegImg, y_len, x_len, z_len, &NC,refineRange);
		std::cerr<<"    Graph Coloring done with "<<NC<<" colors"<<std::endl;

		std::cerr<<"    Starting alpha-expansion..";
		//Call the alpha expansion module		
		//memset(subsegImg/*destination*/,0/*value*/,x_len*y_len*z_len*sizeof(int)/*num bytes to move*/);
		//memcpy(subsegImg/*destination*/,subclustImg/*source*/,x_len*y_len*z_len*sizeof(int)/*num bytes to move*/);
		start_alpha_expansion(subDataImg, subsegImg, Dterms, y_len, x_len, z_len, NC+1);	
		
		//Try this for now: Set the label at each pixel based on the maximum probability (min D-trerm)
		//ind = 0;
		//float tmp_mn=100.0;		
		//int tmp_lbl;
		//int H[8];
		//for(int jk=0; jk<8; jk++)
		//	H[jk] = 0;
		//for(int k=0; k<z_len; k++)		
		//{		
		//	for(int j=0; j<y_len; j++)			
		//	{															
		//		for(int i=0; i<x_len; i++)		
		//		{					
		//			tmp_mn = Dterms/*[(k*x_len*y_len)+(i*y_len)+j]*/[(i+j*x_len+k*x_len*y_len)*(NC+1)];
		//			tmp_lbl = 0;
		//			for(int h=1; h<=NC; h++)
		//			{																	
		//				if(Dterms/*[(h*x_len*y_len*z_len)+(k*x_len*y_len)+(i*y_len)+j]*/[(i+j*x_len+k*x_len*y_len)*(NC+1)+h]<tmp_mn)
		//				{
		//					tmp_lbl = h;
		//					tmp_mn = Dterms[(i+j*x_len+k*x_len*y_len)*(NC+1)+h];
		//				}
		//			}
		//			H[tmp_lbl]++;
		//			subsegImg[ind] = tmp_lbl*20;					
		//			ind++;
		//		}
		//	}			
		//}
					

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
		std::cerr<<"Done"<<std::endl;
		delete [] sublogImg;
		delete [] subsegImg;
		delete [] subclustImg;
		delete [] subDataImg;
	}		

	//relabel the cells
	int numOfObjs = /*getConnCompImage*/getRelabeledImage(segImagePtr, 26, 25, numRows, numColumns,numStacks, 1);	
	std::cerr << "done with " << numOfObjs<<" found"<<std::endl;
	std::cerr << "Creating Final Label Image" << std::endl;		
}

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

//added by Yousef on 8-5-2008
//This function reads parameters from a .ini file following specific format
void yousef_nucleus_seg::readParametersFromFile(const char* pFname)
{
	char achBuffer[1024];
	char achBuffer2[1024];
	char* pchStr = NULL;
	int iCounter = 0;
	int m_iNumOfElements = 0;
	int params[11];//modified by yousef on 11/4/2008
	params[0]=0;
	params[1]=30;
	params[2]=5;
	params[3]=8;
	params[4]=5;
	params[5]=2;
	params[6]=1;
	params[7]=3;
	params[8]=1;
	//added by yousef on 11/4/2008
	params[9]=6;
	//added by yousef on 12/5/2008
	params[10]=100;
	std::ifstream inFile(pFname);
	if (! inFile)
	{
		cout << "Fatal Error: Could not open parameters file " << pFname
			<< " .... Terminating Program." << endl;
		exit(0);
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
		else if(!strcmp(m_pData[i].m_pName,"LoG_size"))
			params[1] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"min_scale"))
			params[2] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"max_scale"))
			params[3] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"xy_clustering_res"))
			params[4] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"z_clustering_res"))
			params[5] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"finalize_segmentation"))
			params[6] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"sampling_ratio_XY_to_Z"))
			params[7] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"Use_Distance_Map"))
			params[8] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"refinement_range"))
			params[9] = m_pData[i].m_pValue;
		else if(!strcmp(m_pData[i].m_pName,"min_object_size"))
			params[10] = m_pData[i].m_pValue;
		else
			continue;
	}

	setParams(params);

}

//Added by Yousef on 7-8-2008
int yousef_nucleus_seg::getRelabeledImage(int *IM, int connectivity, int minSize, int r, int c, int z, int runConnComp)
{
	typedef    int     InputPixelType;
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
	//typedef itk::ConnectedComponentImageFilter< InputImageType, OutputImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelType;
	RelabelType::Pointer relabel = RelabelType::New();

	if(runConnComp == 1)
	{				
		//Compute the labeled connected component image		
		filter->SetInput (im);
		filter->SetFullyConnected( connectivity );	
		//use the connected component image as the input to the relabel component filter		
		relabel->SetInput( filter->GetOutput() );
	}
	else
	{
		//use the input image as the input to the relabel component filter 		
		relabel->SetInput( im );
	}

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
    }
	
	//write the output of the labeling CC filter into our input image
	IteratorType iterator2(relabel->GetOutput(),relabel->GetOutput()->GetRequestedRegion());
	for(int i=0; i<r*c*z; i++)
	{		
		if(binImagePtr[i] == 0)
			IM[i] = 0;
		else
			IM[i] = iterator2.Get()-1;		
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
	  segImagePtr = new int[numStacks*numRows*numColumns];

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

