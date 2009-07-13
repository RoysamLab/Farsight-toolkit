//seedsDetection_2D.cxx

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "seedsdetection.h"

//,.,.

// detectSeeds2D( imgPtr, logImagePtr, seedImagePtr, numRows, numColumns, scaleMin, scaleMax, regionXY, binImagePtr );
int detectSeeds2D( float* IM, float* IM_out, int* IM_bin, int r, int c, double sigma_min, double sigma_max, double scale, int* bImg)
{
    
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
	
	for(int i=0; i<r*c; i++)
	{		
		if(bImg[i]>0)
			iterator1.Set(255.0);
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
	}
    
	
	
	//Start from sigma_min to sigma sigma_max	
	double conv = 0;	
	double sigma = sigma_min;
	float *IMG_tmp = (float *) malloc(r*c*sizeof(float));	

	//Detecting seeds at min scale
	detect_seeds(im,r,c,sigma_min,IM_out);
	
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

		for(int i=0; i<r*c; i++)
		{
			if(sigma<=dImg[i]*2)
				IM_out[i] = (IM_out[i]>=IMG_tmp[i])? IM_out[i] : IMG_tmp[i];				
		}		
	}
	free(IMG_tmp);
	free(dImg);
	
    //get seed points (local maxima points in the LoG response image
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
    IMG[i] = sigma*iterate.Get();
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


void Detect_Local_MaximaPoints(float* im_vals, int r, int c, double scale, int* out1)
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
            min_r = (int) max(0.0,i-scale);
            min_c = (int) max(0.0,j-scale);
            max_r = (int)min((double)r-1,i+scale);
            max_c = (int)min((double)c-1,j+scale);                         
            float mx = get_maximum(im, min_r, max_r, min_c, max_c);
            if(im[i][j] == mx)  
			{
				IND = IND+1;
                out1[(i*c)+j]=255;//(unsigned char)IND; //same as sub2ind(size(im),i,j);                                         
			}
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
  //typedef itk::ApproximateSignedDistanceMapImageFilter<InputImageType, OutputImageType > DTFilter ;
  typedef itk::DanielssonDistanceMapImageFilter <InputImageType, OutputImageType > DTFilter ;
  DTFilter::Pointer dt_obj= DTFilter::New() ;
  dt_obj->SetInput(im) ;
  //dt_obj->SetInsideValue(0.0);
  //dt_obj->SetOutsideValue(255.0);

  try{
	 dt_obj->Update() ;
  }
  catch( itk::ExceptionObject & err ){
	//std::cerr << "Error calculating distance transform: " << err << endl ;
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
