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
	//typedef itk::ImageDuplicator< ImageType3D > DuplicatorType; 
	//DuplicatorType::Pointer duplicator = DuplicatorType::New(); 
	//duplicator->SetInputImage(img); 
	//duplicator->Update(); 
	//myImg = duplicator->GetOutput();
	myImg=img;
	SetImage(myImg);
}

Preprocess::Preprocess()
{		
	//will need SetImage called to function correctly
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
	SetImage(myImg);
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

	ImageType3D::IndexType start = {{ 0,0,0 }};
	ImageType3D::SizeType  size = {{ size1, size2, size3 }};
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
	SetImage(myImg);
}

void Preprocess::SetImage(ImageType3D::Pointer img){	
	//resets the beginning of the pipeline to the image we wish to run it on right now
	myImg=img;
}

std::map<std::string, std::string> Preprocess::filterMap = Preprocess::CreateFilterMap();
std::map<std::string, std::string> Preprocess::CreateFilterMap()
{
	std::map<std::string, std::string> tempMap;

	tempMap["LaplacianOfGaussian"] = "<LaplacianOfGaussian sigma=\"10\" min=\"0\" />\n";
	tempMap["InvertIntensity"] = "<InvertIntensity />\n";
	tempMap["DownSample"] = "<DownSample />\n";
	tempMap["OtsuBinarize"] = "<OtsuBinarize num_thresholds=\"2\" num_in_foreground=\"1\" fgrnd_dark=\"0\" />\n";
	tempMap["ManualThreshold"] =  "<ManualThreshold threshold=\"22\" binary=\"0\" />\n";
	tempMap["RemoveConnectedComponents"] = "<RemoveConnectedComponents minObjSize=\"10\" />\n";
	tempMap["BinaryThinning"] = "<BinaryThinning />\n";
	tempMap["DanielssonDistanceMap"] = "<DanielssonDistanceMap />\n";
	tempMap["MedianFilter"] = "<MedianFilter radiusX=\"3\" radiusY=\"3\" radiusZ=\"0\" />\n";
	tempMap["MinErrorThresholding"] = "<MinErrorThresholding />\n";
	tempMap["GraphCutBinarize"] = "<GraphCutBinarize zyDivs=\"1\" />\n";
	tempMap["OpeningFilter"] = "<OpeningFilter radius=\"3\" />\n";
	tempMap["ClosingFilter"] = "<ClosingFilter radius=\"3\" />\n";
	tempMap["CannyEdgeDetection"] = "<CannyEdgeDetection variance=\"1.0\" upperThreshold=\"6\" lowerThreshold=\"3\" />\n";
	tempMap["DiscreteGaussian"] = "<DiscreteGaussian varX=\"1.0\" varY=\"1.0\" varZ=\"1.0\" maxError=\"0.1\" />\n";
	tempMap["SaveVTKPoints"] = "<SaveVTKPoints filename=\"points.vtk\" xyFactor=\"1\" min=\"255\" max=\"255\" />\n";

	return tempMap;
}

void Preprocess::RunPipe(std::string filename)
{
	TiXmlDocument doc;
	if ( !doc.LoadFile( filename.c_str() ) )
		return;

	TiXmlElement* rootElement = doc.FirstChildElement();
	const char* docname = rootElement->Value();
	if ( strcmp( docname, "Preprocess" ) != 0 )
		return;

	TiXmlElement* parentElement = rootElement->FirstChildElement();
	while (parentElement)
	{
		const char * parent = parentElement->Value();
		if ( strcmp( parent, "LaplacianOfGaussian" ) == 0 )
		{
			int sigma = 10;
			int min = 0;
			parentElement->QueryIntAttribute("sigma", &sigma);
			parentElement->QueryIntAttribute("min", &min);
			std::cout << "Starting LOG...";
			this->LaplacianOfGaussian(sigma, min);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "InvertIntensity" ) == 0 )
		{
			std::cout << "Starting InvertIntensity...";
			this->InvertIntensity();
			std::cout << "done\n";
		}
		else if( strcmp( parent, "DownSample" ) == 0 )
		{
			std::cout << "Starting DownSample...";
			this->DownSample();
			std::cout << "done\n";
		}
		else if( strcmp( parent, "OtsuBinarize" ) == 0 )
		{
			int numThresh = 2, numFore = 1;
			int fgrndDark = 0;
			parentElement->QueryIntAttribute("num_thresholds",&numThresh);
			parentElement->QueryIntAttribute("num_in_foreground",&numFore);
			parentElement->QueryIntAttribute("fgrnd_dark", &fgrndDark);
			std::cout << "Starting OtsuBinarize...";
			this->OtsuBinarize(numThresh,numFore, (bool)fgrndDark);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "ManualThreshold" ) == 0 )
		{
			int threshold = 128;
			int binary = 0;
			parentElement->QueryIntAttribute("threshold", &threshold);
			parentElement->QueryIntAttribute("binary", &binary);
			std::cout << "Starting Manual Threshold of " << threshold << "...";
			this->ManualThreshold(threshold, (bool)binary);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "RemoveConnectedComponents" ) == 0 )
		{
			int minObjSize = 1000;
			parentElement->QueryIntAttribute("minObjSize", &minObjSize);
			std::cout << "Starting RemoveCC...";
			this->RemoveConnectedComponents(minObjSize);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "BinaryThinning" ) == 0 )
		{
			std::cout << "Starting BinaryThinning...";
			this->BinaryThinning();
			std::cout << "done\n";
		}
		else if( strcmp( parent, "DanielssonDistanceMap" ) == 0 )
		{
			std::cout << "Starting DanielssonDistanceMap...";
			this->DanielssonDistanceMap();
			std::cout << "done\n";
		}
		else if( strcmp( parent, "MedianFilter" ) == 0 )
		{
			int radiusX=2, radiusY=3, radiusZ=0;
			parentElement->QueryIntAttribute("radiusX",&radiusX);
			parentElement->QueryIntAttribute("radiusY",&radiusY);
			parentElement->QueryIntAttribute("radiusZ",&radiusZ);
			std::cout << "Starting MedianFilter...";
			this->MedianFilter(radiusX,radiusY,radiusZ);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "MinErrorThresholding" ) == 0 )
		{
			float alpha_B, alpha_F, P_I;
			std::cout << "Starting MinErrorThresholding...";
			this->MinErrorThresholding(&alpha_B, &alpha_F, &P_I);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "GraphCutBinarize" ) == 0 )
		{
			int xyDivs=1;//, zDivs=1;
			parentElement->QueryIntAttribute("xyDivs",&xyDivs);
			std::cout << "Starting GraphCutBinarize...";
			this->GraphCutBinarize(false,xyDivs);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "OpeningFilter" ) == 0 )
		{
			int radius=3;
			parentElement->QueryIntAttribute("radius",&radius);
			std::cout << "Starting OpeningFilter...";
			this->OpeningFilter(radius);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "ClosingFilter" ) == 0 )
		{
			int radius=3;
			parentElement->QueryIntAttribute("radius",&radius);
			std::cout << "Starting ClosingFilter...";
			this->ClosingFilter(radius);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "CannyEdgeDetection" ) == 0 )
		{
			float variance=1.0;
			float upperThreshold = 6;
			float lowerThreshold = 3;
			parentElement->QueryFloatAttribute("variance",&variance);
			parentElement->QueryFloatAttribute("upperThreshold",&upperThreshold);
			parentElement->QueryFloatAttribute("lowerThreshold",&lowerThreshold);
			std::cout << "Starting CannyEdgeDetection...";
			this->CannyEdgeDetection(variance, upperThreshold, lowerThreshold);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "DiscreteGaussian") == 0 )
		{
			float varX=1.0, varY=1.0, varZ=1.0, maxError=0.1;
			parentElement->QueryFloatAttribute("varX",&varX);
			parentElement->QueryFloatAttribute("varY",&varY);
			parentElement->QueryFloatAttribute("varZ",&varZ);
			parentElement->QueryFloatAttribute("maxError",&maxError);
			std::cout << "Starting DiscreteGaussianFilter...";
			this->DiscreteGaussianFilter(varX, varY, varZ, maxError);
			std::cout << "done\n";
		}
		else if( strcmp( parent, "SaveVTKPoints") == 0 )
		{
			const char * filename = parentElement->Attribute("filename");
			float xyFactor = 1.0;
			int min=255, max=255;
			parentElement->QueryFloatAttribute("xyFactor", &xyFactor);
			parentElement->QueryIntAttribute("min",&min);
			parentElement->QueryIntAttribute("max",&max);
			std::cout << "Saving VTK Points...";
			if(!filename)
				this->SaveVTKPoints("points.vtk", xyFactor, min, max);
			else
				this->SaveVTKPoints(filename, xyFactor, min, max);
			std::cout << "done\n";
		}

		parentElement = parentElement->NextSiblingElement();
	} // end while(parentElement)
	//doc.close();
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

//Cuts the size of the image in half!!!
void Preprocess::DownSample()
{
	int size1 = (int)myImg->GetLargestPossibleRegion().GetSize()[0];
	int size2 = (int)myImg->GetLargestPossibleRegion().GetSize()[1];
	int size3 = (int)myImg->GetLargestPossibleRegion().GetSize()[2];

	int newSize1 = (int)(size1/2);
	int newSize2 = (int)(size2/2);
	int newSize3 = size3;

	ImageType3D::Pointer newImg = ImageType3D::New();
	ImageType3D::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	newImg->SetOrigin( origin );

	ImageType3D::IndexType start = {{ 0,0,0 }};
	ImageType3D::SizeType  size = {{ newSize1, newSize2, newSize3 }};
	ImageType3D::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	newImg->SetRegions( region ); 
	newImg->Allocate();

	for(int z=0; z<newSize3; ++z)
	{
		for(int y=0; y<newSize2; ++y)
		{
			for(int x=0; x<newSize1; ++x)
			{
				ImageType3D::IndexType newI = {{ x,y,z }};
				ImageType3D::IndexType oldI = {{ 2*x,2*y,z }};
				newImg->SetPixel( newI, myImg->GetPixel(oldI) );
			}
		}
	}

	myImg = newImg;
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

void Preprocess::DiscreteGaussianFilter(float varX, float varY, float varZ, float maxError)
{
	typedef itk::DiscreteGaussianImageFilter<ImageType3D, FloatImageType3D> FilterType;
	FilterType::Pointer filter = FilterType::New();

	FilterType::ArrayType maxErr;
    maxErr.Fill(maxError);
    filter->SetMaximumError( maxErr );

	FilterType::ArrayType variance;
	variance[0] = varX;
	variance[1] = varY;
	variance[2] = varZ;
	filter->SetVariance(variance);

	filter->SetInput( myImg );

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

void Preprocess::InvertIntensity(void)
{
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

void Preprocess::OpeningFilter( int radius )
{
	typedef itk::BinaryBallStructuringElement<PixelType,3> StructuringElementType;
	typedef itk::GrayscaleMorphologicalOpeningImageFilter<ImageType3D, ImageType3D, StructuringElementType > OpenFilterType;
	OpenFilterType::Pointer grayscaleOpen = OpenFilterType::New();
	StructuringElementType structuringElement;
	structuringElement.CreateStructuringElement();
	structuringElement.SetRadius(radius);
	grayscaleOpen->SetKernel( structuringElement );	 
	grayscaleOpen->SetInput(myImg);
	try
	{
		grayscaleOpen->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Exception caught: " << err << std::endl;
	}
	myImg = grayscaleOpen->GetOutput();
}

void Preprocess::ClosingFilter( int radius )
{
	typedef itk::BinaryBallStructuringElement<PixelType,3> StructuringElementType;
	typedef itk::GrayscaleMorphologicalClosingImageFilter<ImageType3D, ImageType3D, StructuringElementType > FilterType;
	FilterType::Pointer close = FilterType::New();
	StructuringElementType structuringElement;
	structuringElement.CreateStructuringElement();
	structuringElement.SetRadius(radius);
	close->SetKernel( structuringElement );	 
	close->SetInput(myImg);
	try
	{
		close->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Exception caught: " << err << std::endl;
	}
	myImg = close->GetOutput();
}

void Preprocess::CannyEdgeDetection(float variance, float upperThreshold, float lowerThreshold)
{
	typedef itk::CastImageFilter< ImageType3D, FloatImageType3D > CastFilterType;
	CastFilterType::Pointer cast = CastFilterType::New();
	cast->SetInput( myImg );

	typedef itk::CannyEdgeDetectionImageFilter< FloatImageType3D, FloatImageType3D > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( cast->GetOutput() );
	filter->SetUpperThreshold(upperThreshold);		//Threshold for detected edges = threshold
	filter->SetLowerThreshold(lowerThreshold);		//Threshold for detected edges = threshold/2
	//filter->SetThreshold(threshold);		//Lowest allowed value in the output image
	filter->SetVariance(variance);			//For Gaussian smoothing
	//filter->SetMaximumError(.01f);		//For Gaussian smoothing

	try
	{
		filter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "Exception caught: " << err << std::endl;
	}

	myImg = RescaleFloatToImageType( filter->GetOutput() );
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
	
	unsigned short initNumObjects = ccfilter->GetObjectCount();
	unsigned short finalNumObjects = relabel->GetNumberOfObjects();
	unsigned short removedNumObjects = initNumObjects - finalNumObjects;
	std::cout << "Removed " << removedNumObjects << " of " << initNumObjects << " objects...";

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
	typedef itk::BinaryThinningImageFilter3D< ImageType3D, ImageType3D > FilterType;
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

void Preprocess::SaveVTKPoints(std::string filename, float xyFactor, int min, int max)
{
	std::vector<ImageType3D::IndexType> nodes;

	typedef itk::ImageRegionIteratorWithIndex< ImageType3D > IteratorType;
	IteratorType it( myImg, myImg->GetLargestPossibleRegion() );
	for( it.GoToBegin(); !it.IsAtEnd(); ++it )
	{
		ImageType3D::PixelType pix = it.Get();
		if( pix >= min && pix <= max )
		{
			nodes.push_back( it.GetIndex() );
		}
	}

	FILE * fout = fopen(filename.c_str(), "w");
	if (fout == NULL)
		return;

	fprintf(fout, "# vtk DataFile Version 3.0\n");
	fprintf(fout,"Points between %d and %d, using xyFactor %f\n", min, max, xyFactor);
	fprintf(fout,"ASCII\n");
	fprintf(fout,"DATASET POLYDATA\n");

	int num_nodes = (int)nodes.size();
	fprintf(fout,"POINTS %d float\n",num_nodes);
	for(int i=0; i<num_nodes; ++i)
	{
		ImageType3D::IndexType nd = nodes.at(i);
		fprintf(fout,"%f %f %f\n", (float)nd[0]*xyFactor, (float)nd[1]*xyFactor, (float)nd[2]);
	}

	fprintf(fout,"VERTICES %d %d\n", num_nodes, num_nodes*2);
	for(int i=0; i<num_nodes; ++i)
	{
		fprintf(fout,"%d %d\n", 1, i);
	}

	fclose(fout);
}

//void Preprocess::MinErrorThresholding(float *alpha_B, float *alpha_A, float *P_I)
void Preprocess::MinErrorThresholding(float *alpha_B, float *alpha_A, float *P_I, bool overwrite)
{
	//Binarize
	typedef itk::MinErrorThresholdImageFilter< ImageType3D, ImageType3D >  FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( myImg );
	filter->SetNumberOfHistogramBins(256);
	try
	{
		filter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ITK FILTER ERROR: " << err << std::endl;
	}

	*alpha_B = (float)filter->GetAlphaLeft();
	*alpha_A = (float)filter->GetAlphaRight();
	*P_I = (float)filter->GetPriorLeft();
	
	if(overwrite)
		myImg = filter->GetOutput();
}

void Preprocess::GraphCutBinarize(bool shiftDown, int xyDivs)
{
	float alpha_B, alpha_F, P_I;
	MinErrorThresholding(&alpha_B, &alpha_F, &P_I, false);

	//Some times you need to shift the means down.
	if(shiftDown)
	{
		alpha_F = std::max((alpha_F)/2,((alpha_B)+(alpha_F))/2);
		alpha_B = (alpha_B)/1.5;
	}

	//Compute the poisson probs 
	double F_H[256], B_H[256];
	for(int i=0; i<256; i++)
    {
		if( i >= alpha_F )
			F_H[i] = (1-P_I)*ComputePoissonProb((int)alpha_F,alpha_F);
		else
			F_H[i] = (1-P_I)*ComputePoissonProb(i,alpha_F);

		if( i <= alpha_B )
			B_H[i] = P_I*ComputePoissonProb(int(alpha_B),alpha_B);
		else
			B_H[i] = P_I*ComputePoissonProb(i,alpha_B);
    }

	ImageType3D::SizeType size = myImg->GetLargestPossibleRegion().GetSize();

	for(int i=0; i<(int)size[0]; i+=(int)size[0]/xyDivs)
	{
		for(int j=0; j<(int)size[1]; j+=(int)size[1]/xyDivs)
		{
			ImageType3D::IndexType rIndexStart;
			rIndexStart[0] = i;
			rIndexStart[1] = j;
			rIndexStart[2] = 0;
			ImageType3D::SizeType rIndexEnd;		//Plus 1
			rIndexEnd[0] = (int)(i + size[0]/xyDivs + 1);
			rIndexEnd[1] = (int)(j + size[1]/xyDivs + 1);
			if(rIndexEnd[0] > size[0])
				rIndexEnd[0] = size[0];
			if(rIndexEnd[1] > size[1])
				rIndexEnd[1] = size[1];

			ImageType3D::SizeType rSize;
			rSize[0] = rIndexEnd[0] - rIndexStart[0];
			rSize[1] = rIndexEnd[1] - rIndexStart[1];
			rSize[2] = size[2];

			//long num_nodes = size[0]*size[1]*size[2];
			//long num_edges = 3*(size[0]-1)*(size[1]-1)*(size[2]-1);
			long num_nodes = rSize[0]*rSize[1]*rSize[2];
			long num_edges = 3*(rSize[0]-1)*(rSize[1]-1)*(rSize[2]-1);

			//Construct a graph:
			typedef Graph_B<short,short,short> GraphType;
			GraphType * g = new GraphType(/*estimated # of nodes*/ num_nodes, /*estimated # of edges*/ num_edges);

			ImageType3D::RegionType rRegion;
			rRegion.SetIndex(rIndexStart);
			rRegion.SetSize(rSize);

			typedef itk::ImageRegionIteratorWithIndex< ImageType3D > IteratorType;
			IteratorType it( myImg, rRegion );

			//ADD NODES:
			for( it.GoToBegin(); !it.IsAtEnd(); ++it )
			{
				int intst = (int)it.Get();
				ImageType3D::IndexType index = it.GetIndex();
				//int curr_node = (index[2]*size[1]*size[0])+(index[1]*size[0])+index[0];
				int curr_node = ((index[2]-rIndexStart[2])*rSize[1]*rSize[0])+((index[1]-rIndexStart[1])*rSize[0])+(index[0]-rIndexStart[0]);

				//First Add Nodes
				double Df = -log( F_H[intst] );  //it was multiplied by .5                              
				if(Df > 1000.0)
					Df = 1000;
				double Db = -log( B_H[intst] );                  
				if(Db > 1000.0)
					Db=1000;     			
					        				
				g->add_node();
				g->add_tweights( curr_node, /* capacities */ Df, Db );
			}

			//ADD EDGES:
			double sig = 30.0; //30.0;
			double w = 10.0; //10.0;
			for( it.GoToBegin(); !it.IsAtEnd(); ++it )
			{
				ImageType3D::IndexType index = it.GetIndex();
				//int curr_node = (index[2]*size[1]*size[0])+(index[1]*size[0])+index[0];
				int curr_node = ((index[2]-rIndexStart[2])*rSize[1]*rSize[0])+((index[1]-rIndexStart[1])*rSize[0])+(index[0]-rIndexStart[0]);
				int intst = (int)it.Get();

				for(int i=0; i<3; ++i)
				{
					ImageType3D::IndexType index2 = index;
					index2[i]++;
					if((int)index2[i] < (int)rSize[i])
					{
						//int nbr_node = (index2[2]*size[1]*size[0])+(index2[1]*size[0])+index2[0];
						int nbr_node = ((index2[2]-rIndexStart[2])*rSize[1]*rSize[0])+((index2[1]-rIndexStart[1])*rSize[0])+(index2[0]-rIndexStart[0]);
						int intst2 = myImg->GetPixel(index2);
						double Dn = w*exp(-pow((double)intst-intst2,2)/(2*pow(sig,2)));				
						g -> add_edge( curr_node, nbr_node,    /* capacities */  Dn, Dn );
					}
				}
			}

			//Compute the maximum flow:
			g->maxflow();		//Alex DO NOT REMOVE

			//Now Binarize
			for( it.GoToBegin(); !it.IsAtEnd(); ++it )
			{
				ImageType3D::IndexType index = it.GetIndex();
				//int curr_node = (index[2]*size[1]*size[0])+(index[1]*size[0])+index[0];
				int curr_node = ((index[2]-rIndexStart[2])*rSize[1]*rSize[0])+((index[1]-rIndexStart[1])*rSize[0])+(index[0]-rIndexStart[0]);

				if( g->what_segment(curr_node) == GraphType::SOURCE )
					it.Set(0);
				else
					it.Set(255);
			}

			delete g;

		}//end for j
	}//end for i
}

double Preprocess::ComputePoissonProb(double intensity, double alpha)
{
    /*this is the equation
      P = (alpha^intensity)*exp(-alpha)/factorial(intensity);
      however, since alpha and the intensity could be large, computing P in that
      way will result in infinity values from (alpha^intensity) and
      factorial(intensity) as a result of matlab's limitations of the data types*/

    //here is the solution
    double A, P;
    A = exp(-alpha);
    P = 1;
    for (int i=1; i<= intensity; i++)
        P = P * (alpha/i);
    
    P = P*A;

	if(P < std::numeric_limits<long double>::epsilon())
		P = std::numeric_limits<long double>::epsilon();
    
    return P;
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

	ImageType3D::IndexType extractStart = {{ 0,0,slice }};
	ImageType3D::SizeType extractSize = {{ size1, size2, 0 }};
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

	ImageType3D::IndexType index = {{ 0,0,0 }};
	ImageType3D::SizeType size = {{ size1, size2, 1 }};
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








