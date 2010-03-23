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
#include "ftkPreprocess2.h"

namespace ftk
{

Preprocess::Preprocess(ImageType3D::Pointer img)
{
	myImg = img;
}

Preprocess::Preprocess(RGBImageType3D::Pointer img)
{
	typedef itk::RGBToLuminanceImageFilter< RGBImageType3D, ImageType3D > ConvertFilterType;
	ConvertFilterType::Pointer convert = ConvertFilterType::New();
	convert->SetInput( img );

	try
	{
		convert->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
	}

	myImg = convert->GetOutput();
}

Preprocess::Preprocess(RGBImageType3D::Pointer img, const char color)
{
	int size1 = (int)img->GetLargestPossibleRegion().GetSize()[0];
	int size2 = (int)img->GetLargestPossibleRegion().GetSize()[1];
	int size3 = (int)img->GetLargestPossibleRegion().GetSize()[2];

	int component = 0;
	switch(color)
	{
	case 'r':
		component = 0;
		break;
	case 'g':
		component = 1;
		break;
	case 'b':
		component = 2;
		break;
	}

	myImg = ImageType3D::New();
	ImageType3D::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	myImg->SetOrigin( origin );

	ImageType3D::IndexType start = { 0,0,0 };
	ImageType3D::SizeType  size = { size1, size2, size3 };
	ImageType3D::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	myImg->SetRegions( region ); 
	myImg->Allocate();

	typedef itk::ImageRegionIterator< RGBImageType3D > rgbIteratorType;
	rgbIteratorType iteratorRGB ( img, img->GetLargestPossibleRegion() );

	typedef itk::ImageRegionIterator< ImageType3D > myIteratorType;
	myIteratorType iteratorIn( myImg, myImg->GetLargestPossibleRegion() );

	for( iteratorRGB.GoToBegin(), iteratorIn.GoToBegin();
		!iteratorRGB.IsAtEnd(), !iteratorIn.IsAtEnd();
		++iteratorRGB, ++iteratorIn )
	{ 
		RGBPixelType p_rgb = iteratorRGB.Get();
		iteratorIn.Set( p_rgb[component] );
	}
}

void Preprocess::RescaleIntensities(int min, int max)
{
	//Rescale the pixel values
    typedef itk::RescaleIntensityImageFilter< ImageType3D, ImageType3D > RescaleFilterType;
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput( myImg );
	rescaleFilter->InPlaceOn();
	rescaleFilter->SetOutputMinimum( min );
    rescaleFilter->SetOutputMaximum( max );

	try
	{
		rescaleFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
	}

	myImg = rescaleFilter->GetOutput();
}

void Preprocess::LaplacianOfGaussian(int sigma, int min)
{
	int size3 = myImg->GetLargestPossibleRegion().GetSize()[2];

	if(size3 == 1)
	{
		typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType2D, FloatImageType2D >  FilterType;
		FilterType::Pointer laplacian = FilterType::New();
		laplacian->SetNormalizeAcrossScale( true );
		laplacian->SetSigma( sigma );
		laplacian->SetInput( this->ExtractSlice( myImg, 0) );

		try
		{
			laplacian->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
			return;
		}

		ImageType2D::Pointer temp = RescaleFloatToImageType( laplacian->GetOutput(), min );
		myImg = SliceTo3D( temp );
	}
	else
	{
		typedef itk::LaplacianRecursiveGaussianImageFilter< ImageType3D, FloatImageType3D >  FilterType;
		FilterType::Pointer laplacian = FilterType::New();
		laplacian->SetNormalizeAcrossScale( true );
		laplacian->SetSigma( sigma );
		laplacian->SetInput( myImg );

		try
		{
			laplacian->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
			return;
		}

		myImg = RescaleFloatToImageType( laplacian->GetOutput(), min );
	}
}

void Preprocess::InvertIntensity(void)
{
	typedef itk::Image< float, 3 > FloatImageType;
	typedef itk::InvertIntensityImageFilter< ImageType3D, ImageType3D > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( myImg );
	filter->InPlaceOn();

	try
	{
		filter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
	}

	myImg = filter->GetOutput();
	
}

void Preprocess::CurvatureAnisotropicDiffusion( double timestep, double conductance, int iterations )
{
	int size3 = myImg->GetLargestPossibleRegion().GetSize()[2];

	if(size3 == 1)
	{
		typedef itk::CurvatureAnisotropicDiffusionImageFilter< ImageType2D, FloatImageType2D >  FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( this->ExtractSlice( myImg, 0) );
		filter->SetTimeStep(timestep);
		filter->SetConductanceParameter(conductance);
		filter->SetNumberOfIterations(iterations);

		try
		{
			filter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
			return;
		}

		ImageType2D::Pointer temp = RescaleFloatToImageType( filter->GetOutput() );
		myImg = SliceTo3D( temp );
	}
	else
	{
		typedef itk::CurvatureAnisotropicDiffusionImageFilter< ImageType3D, FloatImageType3D > FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(myImg);
		filter->SetTimeStep(timestep);
		filter->SetConductanceParameter(conductance);
		filter->SetNumberOfIterations(iterations);

		try
		{
			filter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
			return;
		}

		myImg = RescaleFloatToImageType( filter->GetOutput() );
	}
}
	
void Preprocess::GradientAnisotropicDiffusion( double timestep, double conductance, int iterations )
{
	int size3 = myImg->GetLargestPossibleRegion().GetSize()[2];

	if(size3 == 1)
	{
		typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType2D, FloatImageType2D >  FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( this->ExtractSlice( myImg, 0) );
		filter->SetTimeStep(timestep);
		filter->SetConductanceParameter(conductance);
		filter->SetNumberOfIterations(iterations);

		try
		{
			filter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
			return;
		}

		ImageType2D::Pointer temp = RescaleFloatToImageType( filter->GetOutput() );
		myImg = SliceTo3D( temp );
	}
	else
	{
		typedef itk::GradientAnisotropicDiffusionImageFilter< ImageType3D, FloatImageType3D > FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(myImg);
		filter->SetTimeStep(timestep);
		filter->SetConductanceParameter(conductance);
		filter->SetNumberOfIterations(iterations);

		try
		{
			filter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
			return;
		}

		myImg = RescaleFloatToImageType( filter->GetOutput() );
	}
}

void Preprocess::MedianFilter( int radiusX, int radiusY, int radiusZ)
{
	int size3 = myImg->GetLargestPossibleRegion().GetSize()[2];
	if(size3==1)
		radiusZ = 0;

	ImageType3D::SizeType radius; 
	radius[0] = radiusX; // radius along x 
	radius[1] = radiusY; // radius along y 
	radius[2] = radiusZ; // radius along y 

	typedef itk::MedianImageFilter<ImageType3D,ImageType3D> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetRadius( radius );  
	filter->SetInput(myImg);

	try
    {
		filter->Update();
    }
	catch( itk::ExceptionObject & err )
    {
		std::cerr << "Exception caught: " << err << std::endl;
    }

	myImg = filter->GetOutput();
}

void Preprocess::ManualThreshold( PixelType threshold, bool binary )
{
	//Apply threshold (any thing below threshold is set to zero)
    typedef itk::ImageRegionIterator<ImageType3D> IteratorType;
	IteratorType itr( myImg, myImg->GetLargestPossibleRegion() );
	for(itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
	{
		PixelType val = itr.Get();
		if(val < threshold) 
        {
			itr.Set(0);
        }
		else
		{
			if(binary)
				itr.Set( itk::NumericTraits<PixelType>::max() );
		}
	}
}

//This function removes all objects that are <= minObjSize from the foreground.
//The foreground remains grayscale after this filter
void Preprocess::RemoveConnectedComponents(int minObjSize)
{
	typedef itk::Image< unsigned short, 3 > ShortImageType;
	typedef itk::ConnectedComponentImageFilter< ImageType3D, ShortImageType > CCFilterType;
	typedef itk::RelabelComponentImageFilter< ShortImageType, ShortImageType > RelabelType;

	CCFilterType::Pointer ccfilter = CCFilterType::New();
	RelabelType::Pointer relabel = RelabelType::New();
	
	ccfilter->SetInput( myImg );
	ccfilter->FullyConnectedOn();

	relabel->SetInput( ccfilter->GetOutput() );
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
	
	unsigned short numObjects = relabel->GetNumberOfObjects();

	ShortImageType::Pointer ccImage = relabel->GetOutput();

	//Use connected component image as a mask:
	itk::ImageRegionIterator< ShortImageType > itr1( ccImage, ccImage->GetLargestPossibleRegion() );
	itk::ImageRegionIterator< ImageType3D > itr2( myImg, myImg->GetLargestPossibleRegion() );
	for(itr1.GoToBegin(), itr2.GoToBegin() ; !itr1.IsAtEnd(); ++itr1, ++itr2)
	{
		if(itr1.Get() == 0)
		{
			itr2.Set( 0 );
		}
	}
}

void Preprocess::OtsuBinarize(int num_thresholds, int num_in_foreground, bool fgrnd_dark)
{
	//Create histogram:
	typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType3D > HistogramGeneratorType;
	typedef HistogramGeneratorType::HistogramType HistogramType;
	HistogramGeneratorType::Pointer histoGenerator = HistogramGeneratorType::New();
	histoGenerator->SetNumberOfBins( 256 );
	histoGenerator->SetInput( myImg );
	try
	{
		histoGenerator->Compute();
	}
	catch( itk::ExceptionObject & excep )
    {
		std::cerr << "    Histogram Computation: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return;
    }

	//Estimate thresholds:
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;
	CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetNumberOfThresholds( num_thresholds );
	calculator->SetInputHistogram( histoGenerator->GetOutput() );
	try
	{
		calculator->Update();
	}
	catch( itk::ExceptionObject & excep )
    {
		std::cerr << "    Otsu Computation: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return;
    }

	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput();

	//Estimate Binarization thresholds:
	PixelType lowerThreshold,upperThreshold;
	if( fgrnd_dark )	//Do I want to make the foregound of the binary dark:
	{
		CalculatorType::OutputType::const_iterator itNum = thresholdVector.end();
 		for( int i=0; i<(num_thresholds-num_in_foreground+1); ++i ) 
			--itNum;
		upperThreshold = static_cast<PixelType>(*itNum);
		lowerThreshold = itk::NumericTraits<PixelType>::min();
	} 
	else
	{
		CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();
		for( int i=0; i<(num_thresholds-num_in_foreground); ++i ) 
			++itNum;
		lowerThreshold = static_cast<PixelType>(*itNum);
		upperThreshold = itk::NumericTraits<PixelType>::max();
	}

	//std::cerr << "Binarization Thresholds: " << lowerThreshold << "  " << upperThreshold << std::endl;

	typedef itk::BinaryThresholdImageFilter< ImageType3D, ImageType3D >  ThreshFilterType;
	ThreshFilterType::Pointer threshfilter = ThreshFilterType::New();
	threshfilter->SetOutsideValue( 0 );
	threshfilter->SetInsideValue( (int)itk::NumericTraits<PixelType>::max() );
	threshfilter->SetInput( myImg );
	threshfilter->SetLowerThreshold( lowerThreshold );
	threshfilter->SetUpperThreshold( upperThreshold );

	try
	{
		threshfilter->Update();
	}
	catch( itk::ExceptionObject & excep )
    {
		std::cerr << "    Threshold: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return;
    }

	myImg = threshfilter->GetOutput();
}

void Preprocess::VotingHoleFilling(int radiusX, int radiusY, int radiusZ, int iterations)
{
	int size3 = myImg->GetLargestPossibleRegion().GetSize()[2];
	if(size3==1)
		radiusZ = 0;

	ImageType3D::SizeType radius;
	radius[0] = radiusX; 
	radius[1] = radiusY;
    radius[2] = radiusZ;

	typedef itk::VotingBinaryIterativeHoleFillingImageFilter< ImageType3D > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetRadius(radius);
	filter->SetMajorityThreshold(2);
	filter->SetMaximumNumberOfIterations(iterations);
	filter->SetInput(myImg);
    try
	{
		filter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return;
	}
	
	myImg = filter->GetOutput();
}

void Preprocess::MedianHoleFilling(int radiusX, int radiusY, int radiusZ)
{
	int size3 = myImg->GetLargestPossibleRegion().GetSize()[2];
	if(size3==1)
		radiusZ = 0;

	ImageType3D::SizeType radius;
	radius[0] = radiusX; 
	radius[1] = radiusY;
    radius[2] = radiusZ;

	typedef itk::BinaryMedianImageFilter< ImageType3D, ImageType3D > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetRadius(radius);
	filter->SetInput(myImg);
    try
	{
		filter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
		return;
	}
	
	myImg = filter->GetOutput();
}

void Preprocess::DanielssonDistanceMap(void)
{
    typedef itk::DanielssonDistanceMapImageFilter<ImageType3D, FloatImageType3D>  DT_Type;
	DT_Type::Pointer filter = DT_Type::New();
	filter->SetInput( myImg );
	try
	{
		filter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl;
	}
	
	myImg = RescaleFloatToImageType( filter->GetDistanceMap() );
}

void Preprocess::BinaryThinning()
{
	typedef itk::BinaryThinningImageFilter< ImageType3D, ImageType3D > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( myImg );

	//Rescale Output:
	typedef itk::RescaleIntensityImageFilter< ImageType3D, ImageType3D > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetOutputMaximum( 255 );
	rescale->SetOutputMinimum( 0 );
	rescale->SetInput( filter->GetOutput() );
	try
	{
		rescale->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl;
	}

	myImg = rescale->GetOutput();
}

void Preprocess::GradientVectorFlow()
{
	//typedef itk::GradientVectorFlowImageFilter<ImageType3D, ImageType3D> FilterType;
	//FilterType::Pointer filter = FilterType::New();

}

Preprocess::ImageType3D::Pointer Preprocess::RescaleFloatToImageType(FloatImageType3D::Pointer img, int inMin)
{
	//Clamp at min:
	typedef itk::ImageRegionIterator< FloatImageType3D > IteratorType;
	IteratorType it( img, img->GetLargestPossibleRegion() );

	for( it.GoToBegin(); !it.IsAtEnd(); ++it )
	{
		if(it.Get() < inMin)
			it.Set(0.0);
	}

	typedef itk::RescaleIntensityImageFilter< FloatImageType3D, ImageType3D > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetOutputMaximum( 255 );
	rescale->SetOutputMinimum( 0 );
	rescale->SetInput( img );
	rescale->Update();
	return rescale->GetOutput();
}

Preprocess::ImageType2D::Pointer Preprocess::RescaleFloatToImageType(FloatImageType2D::Pointer img, int inMin)
{
	typedef itk::ImageRegionIterator< FloatImageType2D > IteratorType;
	IteratorType it( img, img->GetLargestPossibleRegion() );

	for( it.GoToBegin(); !it.IsAtEnd(); ++it )
	{
		if(it.Get() < inMin)
			it.Set(0.0);
	}

	//Rescale weights:
	typedef itk::RescaleIntensityImageFilter< FloatImageType2D, ImageType2D > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetOutputMaximum( 255 );
	rescale->SetOutputMinimum( 0 );
	rescale->SetInput( img );
	rescale->Update();
	return rescale->GetOutput();
}

Preprocess::ImageType3D::Pointer Preprocess::RescaleFloatToImageType(FloatImageType3D::Pointer img)
{
	//Rescale weights:
	typedef itk::RescaleIntensityImageFilter< FloatImageType3D, ImageType3D > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetOutputMaximum( 255 );
	rescale->SetOutputMinimum( 0 );
	rescale->SetInput( img );
	rescale->Update();
	return rescale->GetOutput();
}

Preprocess::ImageType2D::Pointer Preprocess::RescaleFloatToImageType(FloatImageType2D::Pointer img)
{
	//Rescale weights:
	typedef itk::RescaleIntensityImageFilter< FloatImageType2D, ImageType2D > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetOutputMaximum( 255 );
	rescale->SetOutputMinimum( 0 );
	rescale->SetInput( img );
	rescale->Update();
	return rescale->GetOutput();
}


Preprocess::ImageType2D::Pointer Preprocess::ExtractSlice(ImageType3D::Pointer img, int slice)
{
	int size1 = img->GetLargestPossibleRegion().GetSize()[0];
	int size2 = img->GetLargestPossibleRegion().GetSize()[1];
	int size3 = img->GetLargestPossibleRegion().GetSize()[2];

	if(slice >= size3)
		return NULL;

	ImageType3D::IndexType extractStart = { 0,0,slice };
	ImageType3D::SizeType extractSize = { size1, size2, 0 };
	ImageType3D::RegionType extractRegion;
	extractRegion.SetSize( extractSize );
	extractRegion.SetIndex( extractStart );

	typedef itk::ExtractImageFilter< ImageType3D, ImageType2D > ExtractFilterType;
	ExtractFilterType::Pointer extract = ExtractFilterType::New();
	extract->SetInput( img );
	extract->SetExtractionRegion( extractRegion );

	try
    {
		extract->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl ;
	}
	return extract->GetOutput();
}

Preprocess::ImageType3D::Pointer Preprocess::SliceTo3D(ImageType2D::Pointer img)
{
	int size1 = img->GetLargestPossibleRegion().GetSize()[0];
	int size2 = img->GetLargestPossibleRegion().GetSize()[1];

	ImageType3D::Pointer nImg = ImageType3D::New();

	ImageType3D::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	nImg->SetOrigin( origin );

	ImageType3D::IndexType index = { 0,0,0 };
	ImageType3D::SizeType size = { size1, size2, 1 };
	ImageType3D::RegionType regions;
	regions.SetSize( size );
	regions.SetIndex( index );
	nImg->SetRegions(regions);
	nImg->Allocate();

	typedef itk::ImageRegionIterator< ImageType2D > IteratorType2D;
	typedef itk::ImageRegionIterator< ImageType3D > IteratorType3D;
	IteratorType2D iteratorIn( img, img->GetRequestedRegion() );
	IteratorType3D iteratorOut( nImg, nImg->GetRequestedRegion() );
	for( iteratorIn.GoToBegin(), iteratorOut.GoToBegin();
		!iteratorIn.IsAtEnd(), !iteratorOut.IsAtEnd();
		++iteratorIn, ++iteratorOut )
	{ 
		iteratorOut.Set( iteratorIn.Get() );
	}

	return nImg;
}

} // end namespace ftk








