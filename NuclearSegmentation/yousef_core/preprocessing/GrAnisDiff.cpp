#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "GrAnisDiff.h"

int runGrAnisDiff(unsigned char* imgIn, int r, int c, int z, int iter, int timeStep, int comductance)
{
	
	//pixel types
	typedef unsigned char InputPixelType;
	typedef float OutputPixelType;
	typedef itk::Image< InputPixelType, 3 > InputImageType;
	typedef itk::Image< OutputPixelType, 3 > OutputImageType;
	
	//create an ITK image from the input image
	OutputImageType::Pointer im;
	im = OutputImageType::New();
	OutputImageType::PointType origin;
    origin[0] = 0.; 
    origin[1] = 0.;    
	origin[2] = 0.;    
    im->SetOrigin( origin );

    OutputImageType::IndexType start;
    start[0] =   0;  // first index on X
    start[1] =   0;  // first index on Y    
	start[2] =   0;  // first index on Z    
    OutputImageType::SizeType  size;
    size[0]  = c;  // size along X
    size[1]  = r;  // size along Y
	size[2]  = z;  // size along Z
  
    OutputImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    
 //   double spacing[3];
	//spacing[0] = 1; //spacing along x
	//spacing[1] = 1; //spacing along y
	//spacing[2] = sampl_ratio; //spacing along z

    im->SetRegions( region );
	//im->SetSpacing(spacing);
    im->Allocate();
    im->FillBuffer(0);
	im->Update();	
	
	//copy the input image into the ITK image
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > IteratorType;
	IteratorType iterator1(im,im->GetRequestedRegion());
	
	for(int i=0; i<r*c*z; i++)
	{					
		iterator1.Set((float)imgIn[i]);
		++iterator1;		
	}

	//apply gradient anisotropic diffusion
	typedef itk::GradientAnisotropicDiffusionImageFilter< OutputImageType, OutputImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( im );

	filter->SetNumberOfIterations( iter );
	filter->SetTimeStep( timeStep );
	filter->SetConductanceParameter( comductance );

	try
	{
		filter->Update();
	}
	catch( itk::ExceptionObject & err ) 
    { 
		std::cerr << err << std::endl;
		return 0;
	}

	typedef itk::RescaleIntensityImageFilter< OutputImageType, InputImageType > RescaleFilterType;

	RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
	rescaler->SetOutputMinimum(   0 );
	rescaler->SetOutputMaximum( 255 );
	rescaler->SetInput( filter->GetOutput() );
	rescaler->Update();

	//overwrite the input image for now in order to save space
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType2;
	IteratorType2 iterator2(rescaler->GetOutput(),rescaler->GetOutput()->GetRequestedRegion());
	
	for(int i=0; i<r*c*z; i++)
	{					
		imgIn[i] = iterator2.Get();
		++iterator2;		
	}

	return 1;
}