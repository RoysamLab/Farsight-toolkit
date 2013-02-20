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

//seedsDetection_2D.cxx

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "seedsdetection.h"
//added by Yousef on 9/26/2009
#include "itkExtractImageFilter.h"

#include "itkSignedMaurerDistanceMapImageFilter.h"
//#include "conio.h"
//,.,.

typedef    float     InputPixelType;
typedef itk::Image< InputPixelType,  2 >   InputImageType;

float get_maximum(float** A, int r1, int r2, int c1, int c2);

void Detect_Local_MaximaPoints(float* im_vals, int r, int c, double scale, unsigned short* im_bin);

int distMap(itk::SmartPointer<InputImageType> im, int r, int c, float* IMG);

int detect_seeds(itk::SmartPointer<InputImageType>, int , int , const double, float*);

void estimateMinMaxScales2D(itk::SmartPointer<InputImageType> im, float* distIm, double* minScale, double* maxScale, int r, int c);

unsigned short get_maximumV2(unsigned short* A, int r1, int r2, int c1, int c2, int r, int c);

int computeMedian2D(std::vector< std::vector<unsigned short> > scales, int cntr);

int computeWeightedMedian2D(std::vector< std::vector<float> > scales, int cntr);

InputImageType::Pointer extract2DImageRegion(itk::SmartPointer<InputImageType> im, int sz_x, int sz_y, int start_x, int start_y);

// detectSeeds2D( imgPtr, logImagePtr, seedImagePtr, numRows, numColumns, scaleMin, scaleMax, regionXY, binImagePtr );
int detectSeeds2D( float* IM, float* IM_out, unsigned short* IM_bin, int r, int c, double* sigma_min_in, double* sigma_max_in, double* scale_in, unsigned short* bImg, bool paramEstimation)
{
    //copy the the input parameters into local variables
	double sigma_min = sigma_min_in[0];
	double sigma_max = sigma_max_in[0];
	double scale = scale_in[0];

	//Create an itk image from the input image
	InputImageType::Pointer im;
	im = InputImageType::New();
	InputImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
    im->SetOrigin( origin );

    InputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
    InputImageType::SizeType  size;
    size[0]  = c;  // size along X
    size[1]  = r;  // size along Y
  
    InputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
    im->SetRegions( region );
    im->Allocate();
    im->FillBuffer(0);
	im->Update();
	
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	
	// Copy the binary image array to an itk image:
	for(int i=0; i<r*c; i++)
	{		
		if(bImg[i]>0)
		{
			//iterator1.Set(255.0);
			iterator1.Set(65535.0);
			//std::cout<<"1\n";
		}
		else
			iterator1.Set(0.0);
		++iterator1;	
	}
	


	//Compute Distance Map
	float* dImg = (float *) malloc(r*c*sizeof(float));
	distMap(im, r, c, dImg);
	iterator1.GoToBegin();	
	
	//Copy the input image into an ITK image
	for(int i=0; i<r*c; i++)
	{				
		iterator1.Set(IM[i]);		
		++iterator1;	

		//try this
		if(bImg[i] == 0)
			dImg[i] = 0;
	}

	//By Yousef (9/26/2009)
	//Estimate the segmentation parameters
	if(paramEstimation)
	{
		std::cout<<"Estimating parameters..."<<std::endl;
		//estimateMinMaxScales2D(im, dImg, &sigma_min, &sigma_max, r, c);		
		scale = sigma_min;

		//sigma_min = 5;
		//sigma_max = 10;

		if(scale<3)
			scale = 3; //just avoid very small search boxes		
		std::cout<<"    Minimum scale = "<<sigma_min<<std::endl;
		std::cout<<"    Maximum scale = "<<sigma_max<<std::endl;
		std::cout<<"    Clustering Resolution = "<<scale<<std::endl;
	

		//write out the parameters
		sigma_min_in[0] = sigma_min;
		sigma_max_in[0] = sigma_max;
		scale_in[0] =  scale;		
	}
	

	else
	{
		// FOR THE MODEL BASED NUCLEUS MERGING,I SHATTER SMALL NUCLEI
		// TYPICALLY, THE IMAGE ROWS AND COLUMNS ARE LESS THAN 100 PIXELS
		// USE THIS HEURISTIC TO SET MAX AND MAX SCALES TO 3 AND 4
		// INSERTED BY RAGHAV 
		if(r <100 || c<100)
		{	
		
		//write out the parameters
		sigma_min_in[0] = sigma_min;
		//sigma_max_in[0] = ((sigma_max -scaleDiff) > sigma_min)?(sigma_max -2):sigma_min;
		sigma_max_in[0] =	sigma_max - floor((sigma_max - sigma_min)/2);
		scale_in[0] =  scale;	
		sigma_min = sigma_min_in[0];
		sigma_max = sigma_max_in[0];	
		}

		std::cout<<"    Minimum scale = "<<sigma_min<<std::endl;
		std::cout<<"    Maximum scale = "<<sigma_max<<std::endl;
		std::cout<<"    Clustering Resolution = "<<scale<<std::endl;
	}

	
	//Start from sigma_min to sigma sigma_max	
	double conv = 0;	
	double sigma = sigma_min;
	float *IMG_tmp = (float *) malloc(r*c*sizeof(float));	



	while(!conv)
	{		
		sigma = sigma+1;		
		if(sigma>sigma_max)
		{
			conv=1;
			break;
		}

		//Detecting seeds at next scale
		detect_seeds(im,r,c,sigma,IMG_tmp);
	//	detect_seeds(original_image,r,c,sigma,IM_out);	// change this back


		for(int i=0; i<r*c; i++)
		{
			//if(sigma<=dImg[i]*2)
				IM_out[i] = (IM_out[i]>=IMG_tmp[i])? IM_out[i] : IMG_tmp[i];				
		}		
	}
	free(IMG_tmp);
	free(dImg);
	

    //get seed points (local maxima points in the LoG response image
	//IM_bin = new unsigned short[r*c];
	Detect_Local_MaximaPoints(IM_out, r, c, scale, IM_bin);
			
	return 1;
}

int detect_seeds(itk::SmartPointer<InputImageType> im, int r, int c, const double sigma, float* IMG)
{
  
  //  Software Guide : BeginLatex
  //
  //  Types should be selected on the desired input and output pixel types.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef    float     InputPixelType;
  typedef    float     OutputPixelType;
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  The input and output image types are instantiated using the pixel types.
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::Image< InputPixelType,  2 >   InputImageType;
  typedef itk::Image< OutputPixelType, 2 >   OutputImageType;
  // Software Guide : EndCodeSnippet


  //typedef itk::ImageFileReader< InputImageType >  ReaderType;


  //  Software Guide : BeginLatex
  //
  //  The filter type is now instantiated using both the input image and the
  //  output image types.
  //
  //  \index{itk::RecursiveGaussianImageFilter!Instantiation}
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  typedef itk::LaplacianRecursiveGaussianImageFilter<
                        InputImageType, OutputImageType >  FilterType;
  // Software Guide : EndCodeSnippet


  //ReaderType::Pointer reader = ReaderType::New();
  //reader->SetFileName( in_image_name );


  
  // Software Guide : BeginCodeSnippet
  FilterType::Pointer laplacian = FilterType::New();
  // Software Guide : EndCodeSnippet



  //  Software Guide : BeginLatex
  //  
  //  The option for normalizing across scale space can also be selected in this filter.
  //
  //  \index{LaplacianRecursiveGaussianImageFilter!SetNormalizeAcrossScale()}
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  laplacian->SetNormalizeAcrossScale( true );
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  The input image can be obtained from the output of another
  //  filter. Here, an image reader is used as the source. 
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  laplacian->SetInput( im);//reader->GetOutput() );
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  It is now time to select the $\sigma$ of the Gaussian used to smooth the
  //  data.  Note that $\sigma$ must be passed to both filters and that sigma
  //  is considered to be in millimeters. That is, at the moment of applying
  //  the smoothing process, the filter will take into account the spacing
  //  values defined in the image.
  //
  //  \index{itk::LaplacianRecursiveGaussianImageFilter!SetSigma()}
  //  \index{SetSigma()!itk::LaplacianRecursiveGaussianImageFilter}
  //
  //  Software Guide : EndLatex 

  //const double sigma = atof( argv[3] );

  // Software Guide : BeginCodeSnippet
  laplacian->SetSigma( sigma );
 
  // Software Guide : EndCodeSnippet


  //  Software Guide : BeginLatex
  //
  //  Finally the pipeline is executed by invoking the \code{Update()} method.
  //
  //  \index{itk::LaplacianRecursiveGaussianImageFilter!Update()}
  //
  //  Software Guide : EndLatex 

  // Software Guide : BeginCodeSnippet
  try
    {
    laplacian->Update();
    }
  catch( itk::ExceptionObject & err ) 
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return EXIT_FAILURE;
    } 
  // Software Guide : EndCodeSnippet




  //By Yousef: Now, instead of writing the output, let's copy the resulting image into an array and return it  
  long int i = 0;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
  IteratorType iterate(laplacian->GetOutput(),laplacian->GetOutput()->GetRequestedRegion());
  while ( i<r*c)
  {
    //IMG[i] = sigma*iterate.Get();
	IMG[i] = iterate.Get();
    ++i;
	++iterate;
  }
	

    return EXIT_SUCCESS;

}

float get_maximum(float** A, int r1, int r2, int c1, int c2)
{
    float mx = A[r1][c1];
    for(int i=r1; i<=r2; i++)
    {
        for(int j=c1; j<=c2; j++)
        {
            if(A[i][j]>mx)
                mx = A[i][j];
        }
    }
    return mx;
}

unsigned short get_maximumV2(float* A, int r1, int r2, int c1, int c2, int r, int c)
{
    unsigned short mx = A[0];
    for(int i=r1; i<=r2; i++)
    {
        for(int j=c1; j<=c2; j++)
        {
            if(A[i*c+j]>mx)
                mx = (unsigned short) A[i*c+j];
        }
    }
    return mx;
}


void Detect_Local_MaximaPoints(float* im_vals, int r, int c, double scale, unsigned short* out1)
{  
    float** im;
    int min_r, min_c, max_r, max_c;    
   
    
    
    im = (float **) malloc(r*sizeof(float*));       
    for(int i=0; i<r; i++)
    {
        im[i] = (float *) malloc(c*sizeof(float));
        for(int j=0; j<c; j++)
        {
            im[i][j] = im_vals[(i*c)+j];            
        }
    }    
        

    //start by getting local maxima points
    //if a point is a local maximam give it a local maximum ID
        
	int IND = 0;
    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {					
			min_r = (int) std::max(0.0,i-scale);
            min_c = (int) std::max(0.0,j-scale);
            max_r = (int) std::min((double)r-1,i+scale);
            max_c = (int) std::min((double)c-1,j+scale);                         
            float mx = get_maximum(im, min_r, max_r, min_c, max_c);
            if(im[i][j] == mx)  
			{
				IND = IND+1;
                out1[(i*c)+j]=255;//(unsigned char)IND; //same as sub2ind(size(im),i,j);                                         
			}
			else
				out1[(i*c)+j]=0;
        }
    }  
    	        
}

int distMap(itk::SmartPointer<InputImageType> im, int r, int c, float* IMG)
{
  
  //  Types should be selected on the desired input and output pixel types.
  typedef    float     InputPixelType;
  typedef    float     OutputPixelType;


  //  The input and output image types are instantiated using the pixel types.
  typedef itk::Image< InputPixelType,  2 >   InputImageType;
  typedef itk::Image< OutputPixelType, 2 >   OutputImageType;


  //  The filter type is now instantiated using both the input image and the
  //  output image types.
  /*typedef itk::DanielssonDistanceMapImageFilter <InputImageType, OutputImageType > DTFilter ;
  DTFilter::Pointer dt_obj= DTFilter::New() ;
  dt_obj->SetInput(im) ;*/
  typedef itk::SignedMaurerDistanceMapImageFilter<InputImageType, OutputImageType>  DTFilter;
  DTFilter::Pointer dt_obj= DTFilter::New() ;
  dt_obj->SetInput(im) ;
  dt_obj->SetSquaredDistance( false );      
  dt_obj->SetInsideIsPositive( false );
  

  try{
	 dt_obj->Update() ;
  }
  catch( itk::ExceptionObject & err ){
	  std::cerr << "Error calculating distance transform: " << err << std::endl ;
    return -1;
  }
 
  //   Copy the resulting image into the input array
  long int i = 0;
  typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
  IteratorType iterate(dt_obj->GetOutput(),dt_obj->GetOutput()->GetRequestedRegion());
  while ( i<r*c)
  {
	  IMG[i] = fabs(iterate.Get());	  
      ++i;
 	  ++iterate;
  }
	
 
  return EXIT_SUCCESS;
}

//added by Yousef on 9/26/2009
//Estimate the min and max scales based on the local maxima points of the distance map
void estimateMinMaxScales2D(itk::SmartPointer<InputImageType> im, float* distIm, double* minScale, double* maxScale, int r, int c)
{
	int min_r, min_c, max_r, max_c;
	int II = 0;
	minScale[0] = 1000.0;
	maxScale[0] = 0.0;
	
	std::vector< std::vector<unsigned short> > scales;
	int cnt = 0;
	//ofstream p;
	//p.open("checkme.txt");
	//unsigned short mmx = 0;
	for(int i=1; i<r-1; i++)
    {
        for(int j=1; j<c-1; j++)
        {								
			min_r = (int) std::max(0.0,(double)i-2);
			min_c = (int) std::max(0.0,(double)j-2);
			max_r = (int) std::min((double)r-1,(double)i+2);
			max_c = (int) std::min((double)c-1,(double)j+2);                         			
			unsigned short mx = get_maximumV2(distIm, min_r, max_r, min_c, max_c, r, c);
			
			if(mx <= 2)
				continue; //background or edge point
			II = (i*c)+j;
			if(distIm[II] == mx)    
			{		
				//since we have scaled by 100 earlier, scale back to the original value
				//also, we want to dived by squre root of 2 or approximately 1.4
				//mx = mx/140;							
				
				//add the selected scale to the list of scales
				std::vector <unsigned short> lst;
				lst.push_back(mx);
				lst.push_back(i);
				lst.push_back(j);				
				//p<<j<<" "<<i<<" "<<k<<" "<<mx<<std::endl;
				scales.push_back(lst);
				
				cnt++;														
			}			
        }
    } 
	//p.close();
	
	//get the median of the scales(distances)	
	int medianS = computeMedian2D(scales, cnt); 
	//ofstream p2;
	//p2.open("checkme2.txt");
	//p2<<"med = "<<medianS<<std::endl;				

	//ofstream p3;
	//p3.open("checkme3.txt");	
	//For each local maximum point,try to find the best LoG scale
	//To do that, suppose the distance at a given local maximum point is d, 
	//then compute the its LoG responses at scales from d/2 to d
	//Then, select the scale the gave us the maximum LoG response
	int mnScl = 10000;
	int mxScl = 0;
	int cnt2 = 0;
	std::vector<std::vector<float> > smallScales;
	std::vector<std::vector<float> > largeScales;
	int numSmall = 0;
	int numLarge = 0;

	for(int ind=0; ind<cnt; ind++)
	{
		int mx = scales[ind][0];
		int i = scales[ind][1];
		int j = scales[ind][2];		

		int smin = (int) std::ceil(mx/2.0);
		if(smin == mx)
			continue;
		if(smin == 1)
			smin++;
		cnt2++;
		min_r = (int) std::max(0.0,(double)i-mx);
		min_c = (int) std::max(0.0,(double)j-mx);		
		max_r = (int) std::min((double)r-1,(double)i+mx);
		max_c = (int) std::min((double)c-1,(double)j+mx);                         
					
		int sub_r = i-min_r;
		int sub_c = j-min_c;		
		int sz_r = (max_r-min_r+1);
		int sz_c = (max_c-min_c+1);		
		int ind_i = sub_r*sz_c+sub_c;
													
		InputImageType::Pointer im_Small = extract2DImageRegion(im, sz_c, sz_r, min_c, min_r);															
		float* IMG = new float[sz_c*sz_r];
		float max_resp = -100000.0;	
		int best_scale = 0.0;
		double sigma;	
	
		for(int kk=smin; kk<=mx; kk++)
		{						
			sigma = kk;
			detect_seeds(im_Small, sz_r, sz_c, sigma, IMG);
			
			//Get the scale at which the LoG response at our point of interest is maximum
			if(IMG[ind_i]>=max_resp)
			{
				max_resp = IMG[ind_i];								
				best_scale = kk;				
			}			
		}
		std::vector<float> pp;
		pp.push_back(best_scale);
		pp.push_back(max_resp);

		if(mx<=medianS)
		{	
			numSmall++;		
			smallScales.push_back(pp);
		}
		else
		{
			numLarge++;		
			largeScales.push_back(pp);
		}

		//p3<<j<<" "<<i<<" "<<k<<" "<<mx<<" "<<best_scale<<" "<<max_resp<<std::endl;
		/*mx = best_scale;	
		if(mx<mnScl)
			mnScl = mx;
		if(mx>mxScl)
			mxScl = mx;*/

		delete [] IMG;
	}
	//p3.close();
	//set the min and max scales to the LoG-weighted medians of the small and large scale sets
	mnScl =  computeWeightedMedian2D(smallScales, numSmall); //work on this
	mxScl =  computeWeightedMedian2D(largeScales, numLarge); //work on this

	scales.clear();
	smallScales.clear();
	largeScales.clear();

	//I assume at least 4 scales must be used (will be relaxed later)
	if(mxScl<mnScl+3)
		mxScl = mnScl+3;
	minScale[0] = mnScl;
	maxScale[0] = mxScl;	
	//p2<<"min_scale="<<mnScl<<std::endl;
	//p2<<"max_scale="<<mxScl<<std::endl;
	//p2.close();	
}

int computeMedian2D(std::vector< std::vector<unsigned short> > scales, int cntr)
{
	if(cntr == 1)
		return scales[0][0];

	unsigned short* srtList = new unsigned short[cntr];
	for(int i=0; i<cntr-1; i++)
		srtList[i] = scales[i][0];

	//sort the distances
	for(int i=0; i<cntr-1; i++)
	{
		for(int j=i+1; j<cntr; j++)
		{
			if(srtList[j]<srtList[i])
			{
				unsigned short tmp = srtList[i];
				srtList[i] = srtList[j];
				srtList[j] = tmp;
			}
		}
	}

	//if the number of points is odd then the median is the mid point
	//else, it is in between the two mid points
	int res = cntr % 2;
	int mdn;
	if(res!=0)
	{
		mdn = (cntr+1)/2;
		mdn = (int) srtList[mdn-1];
	}
	else
	{
		int mdn1 = cntr/2;
		mdn = (int) (srtList[mdn1]+srtList[mdn1+1])/2;
	}
		
	return mdn;
}

int computeWeightedMedian2D(std::vector< std::vector<float> > scales, int cntr)
{
	if(cntr == 1)
		return scales[0][0];

	float* srtList = new float[cntr];
	float* wgtList = new float[cntr];
	

	float sumWeights = 0;
	for(int i=0; i<cntr; i++)
	{
		srtList[i] = scales[i][0];
		wgtList[i] = scales[i][1];		
		sumWeights+= wgtList[i];
	}

	//normalize the list of weights
	for(int i=0; i<cntr; i++)
		wgtList[i] /= sumWeights;

	//sort the distances
	for(int i=0; i<cntr-1; i++)
	{
		for(int j=i+1; j<cntr; j++)
		{
			if(srtList[j]<srtList[i])
			{
				float tmp = srtList[i];
				srtList[i] = srtList[j];
				srtList[j] = tmp;
				tmp = wgtList[i];
				wgtList[i] = wgtList[j];
				wgtList[j] = tmp;

			}
		}
	}

	//Find the point at which the cummulative sum exceeds .5
	float cumSum = 0;
	int mdn = 0;

	for(int i=0; i<cntr; i++)
	{
		cumSum+=wgtList[i];
		if(cumSum > .5)
		{
			mdn = srtList[i];
			break;
		}

	}	
		
	delete [] srtList;
	delete [] wgtList;
	return mdn;
}

InputImageType::Pointer extract2DImageRegion(itk::SmartPointer<InputImageType> im, int sz_x, int sz_y, int start_x, int start_y)
{
    typedef itk::ExtractImageFilter< InputImageType, InputImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();   
    
    InputImageType::SizeType size;
	size[0] = sz_x;
	size[1] = sz_y;
	//size[2] = 0.;
    
    InputImageType::IndexType start;
    start[0] = start_x;
	start[1] = start_y;
	//start[2] = 0.;
    
    InputImageType::RegionType desiredRegion;
    desiredRegion.SetSize(  size  );
    desiredRegion.SetIndex( start );
	
    
    filter->SetExtractionRegion( desiredRegion );
    
    filter->SetInput( im );

    InputImageType::Pointer img = filter->GetOutput();
    img->Update();

    return img;
}
