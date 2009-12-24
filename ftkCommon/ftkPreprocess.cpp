/*
 *  ftkPreprocess.cpp
 *  Farsight
 *
 *  Created by RAGHAV on 12/3/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "ftkPreprocess.h"




ftkPreprocess::ftkPreprocess()
{

std::cout<<"I was here"<<std::endl;

}






ftk::Image::Pointer ftkPreprocess::GADiffusion(void)
{
		
	     typedef itk::GradientAnisotropicDiffusionImageFilter<InpImageType,FloatImageType> FilterType;
		 typedef itk::CastImageFilter<FloatImageType,InpImageType> ReverseCastFilter;	
		 FilterType::Pointer filter = FilterType::New();
         filter->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));
		 filter->SetTimeStep(filterParams[0]);
		 filter->SetNumberOfIterations(filterParams[1]);
		 filter->SetConductanceParameter(filterParams[2]);
		
		std::cout << "Applying Anisotropic Diffusion Filter.........: " << std::endl;
		  
		try
		{
			filter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
		std::cerr << "Exception caught: " << err << std::endl;
		}
 
		 FloatImageType::Pointer imf = filter->GetOutput();
		 typedef itk::ImageRegionIterator<FloatImageType> IRI;
		 IRI iter(imf,imf->GetLargestPossibleRegion());
						
		 for(iter.GoToBegin();!iter.IsAtEnd(); ++iter)
			{
				float value = iter.Get();
				value = ((value < 0)?0:((value>255)?255:value));
				iter.Set(value);
			}

			ReverseCastFilter::Pointer rfilter = ReverseCastFilter::New();
			rfilter->SetInput(filter->GetOutput());
			rfilter->Update();
		
		   std::vector<unsigned char> color(3,255);
		   ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
		   filtImg->AppendChannelFromData3D(rfilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
		   myImg = filtImg;			
		   return myImg;								
																	
			/* typedef itk::ImageFileWriter<imagetype > WriterType;
			   WriterType::Pointer writer = WriterType::New();
			   std::string fName = this->FN+"_median"+".tif";
			   writer->SetFileName(fName.c_str());
			   writer->SetInput(mfilter->GetOutput());
			   writer->Update();
	 */	 
			//QApplication::restoreOverrideCursor();			 
	}


ftk::Image::Pointer ftkPreprocess::MedianFilter(void)
{
	 
	 typedef itk::MedianImageFilter<InpImageType,InpImageType> FilterType;
     
	 FilterType::Pointer mFilter = FilterType::New();
	 InpImageType::SizeType indexRadius; 
	 indexRadius[0] = this->filterParams[0]; // radius along x 
	 indexRadius[1] = filterParams[1]; // radius along y 
	 indexRadius[2] = filterParams[2]; // radius along y 
	

	mFilter->SetRadius( indexRadius );  
	mFilter->SetInput(myImg->GetItkPtr<InpPixelType>(0,0) );

    std::cout << "Applying Median Filter.........: " << std::endl;
	
	try
    {
		mFilter->Update();
    }
	catch( itk::ExceptionObject & err )
    {
		std::cerr << "Exception caught: " << err << std::endl;
    }

    
		
	/* typedef itk::ImageFileWriter<imagetype > WriterType;
	WriterType::Pointer writer = WriterType::New();
	std::string fName = this->FN+"_median"+".tif";
	writer->SetFileName(fName.c_str());
	writer->SetInput(mfilter->GetOutput());
	writer->Update();
	 */	 
	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(mFilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;
}



ftk::Image::Pointer ftkPreprocess::SigmoidFilter(void)
{
	
   // Declare the anisotropic diffusion vesselness filter
    typedef itk::SigmoidImageFilter< InpImageType,InpImageType>  SigmoidImFilter;

  // Create a vesselness Filter
    SigmoidImFilter::Pointer SigmoidFilter = SigmoidImFilter::New();
    SigmoidFilter->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));
	SigmoidFilter->SetAlpha(filterParams[0]); 
	SigmoidFilter->SetBeta( filterParams[1] ); 
    SigmoidFilter->SetOutputMinimum( filterParams[2] ); 
    SigmoidFilter->SetOutputMaximum( filterParams[3]); 

	std::cout << "Applying Sigmoid Filter.........: " << std::endl;

  try
    {
    SigmoidFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    }

	// Add the filter to the segmentation view
	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(SigmoidFilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;

    }



ftk::Image::Pointer ftkPreprocess::GrayscaleErode(void)
{
	
	typedef itk::BinaryBallStructuringElement<InpPixelType,3> StructuringElementType;
	typedef itk::GrayscaleErodeImageFilter<InpImageType,InpImageType,StructuringElementType > ErodeFilterType;
	ErodeFilterType::Pointer  grayscaleErode  = ErodeFilterType::New();
	StructuringElementType  structuringElement;
	structuringElement.CreateStructuringElement();
	structuringElement.SetRadius(filterParams[0]);
	grayscaleErode->SetKernel(  structuringElement );	 
	grayscaleErode->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));

  std::cout << "Applying Grayscale Erosion Filter.........: " << std::endl;

  try
    {
    grayscaleErode->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    }

	// Add the filter to the segmentation view
	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(grayscaleErode->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;
}




ftk::Image::Pointer ftkPreprocess::GrayscaleDilate(void)
{
	
	typedef itk::BinaryBallStructuringElement<InpPixelType,3> StructuringElementType;
	typedef itk::GrayscaleDilateImageFilter<InpImageType,InpImageType,StructuringElementType> DilateFilterType;	
	DilateFilterType::Pointer  grayscaleDilate  = DilateFilterType::New();
	StructuringElementType  structuringElement;
	structuringElement.CreateStructuringElement();
	structuringElement.SetRadius(filterParams[0]);
	grayscaleDilate->SetKernel(structuringElement);	 
	grayscaleDilate->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));


  std::cout << "Applying Grayscale Dilation Filter.........: " << std::endl;

  try
    {
    grayscaleDilate->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    }

	// Add the filter to the segmentation view
	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(grayscaleDilate->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	//Update the segmentation view
	myImg = filtImg;
	return myImg;
    }



ftk::Image::Pointer ftkPreprocess::GrayscaleOpen(void)
{
	typedef itk::BinaryBallStructuringElement<InpPixelType,3> StructuringElementType;
	typedef itk::GrayscaleMorphologicalOpeningImageFilter<InpImageType,InpImageType,StructuringElementType >  OpenFilterType;
	OpenFilterType::Pointer  grayscaleOpen  = OpenFilterType::New();
	StructuringElementType  structuringElement;
	structuringElement.CreateStructuringElement();
	structuringElement.SetRadius(filterParams[0]);
	grayscaleOpen->SetKernel(  structuringElement );	 
	grayscaleOpen->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));


	std::cout << "Applying Grayscale Opening Filter.........: " << std::endl;

	try
		{
			grayscaleOpen->Update();
		}
	catch( itk::ExceptionObject & err )
		{
			std::cerr << "Exception caught: " << err << std::endl;
		}

	// Add the filter to the segmentation view
	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(grayscaleOpen->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;

    }



ftk::Image::Pointer ftkPreprocess::GrayscaleClose(void)
{
	typedef itk::BinaryBallStructuringElement<InpPixelType,3> StructuringElementType;
	typedef itk::GrayscaleMorphologicalClosingImageFilter<InpImageType,InpImageType,StructuringElementType > CloseFilterType;

	CloseFilterType::Pointer  grayscaleClose  = CloseFilterType::New();
	StructuringElementType  structuringElement;
	structuringElement.CreateStructuringElement();
	structuringElement.SetRadius(filterParams[0]);
	grayscaleClose->SetKernel(structuringElement );	 
	grayscaleClose->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));


	std::cout << "Applying Grayscale Closing Filter.........: " << std::endl;

	try
		{
			grayscaleClose->Update();
		}
	catch( itk::ExceptionObject & err )
		{
			std::cerr << "Exception caught: " << err << std::endl;
		}

	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(grayscaleClose->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;

    }
 




 ftk::Image::Pointer ftkPreprocess::Resample(void)
{
		
		typedef   float     InternalPixelType;
		typedef itk::Image< InternalPixelType, 3 >   InternalImageType;	
		typedef itk::IntensityWindowingImageFilter< 
                                  InpImageType, 
                                  InternalImageType >  IntensityFilterType;
		
		  typedef itk::RecursiveGaussianImageFilter< 
                                InpImageType,
                                InpImageType > GaussianFilterType;
		
		 GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
		 GaussianFilterType::Pointer smootherY = GaussianFilterType::New();

		 smootherX->SetInput( myImg->GetItkPtr<InpPixelType>(0,0) );
         smootherY->SetInput( smootherX->GetOutput() );
		
		
		 InpImageType::Pointer inputImage = myImg->GetItkPtr<InpPixelType>(0,0);
         const InpImageType::SpacingType& inputSpacing = inputImage->GetSpacing();
		 const double isoSpacing = vcl_sqrt( inputSpacing[2] * inputSpacing[0] );
         smootherX->SetSigma( isoSpacing );
		 smootherY->SetSigma( isoSpacing );
		 
		 
		smootherX->SetDirection( 0 );
		smootherY->SetDirection( 1 );
		smootherX->SetNormalizeAcrossScale( true );
		smootherY->SetNormalizeAcrossScale( true );

		
		typedef   unsigned char   OutputPixelType;
		typedef itk::Image< OutputPixelType,3 >   OutputImageType;
        typedef itk::ResampleImageFilter<
                InpImageType, OutputImageType >  ResampleFilterType;
        ResampleFilterType::Pointer resampler = ResampleFilterType::New();	
				
		typedef itk::IdentityTransform< double,3>  TransformType;
        TransformType::Pointer transform = TransformType::New();
        transform->SetIdentity();
		resampler->SetTransform( transform );
						
		typedef itk::LinearInterpolateImageFunction< 
                          InpImageType, double >  InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		resampler->SetInterpolator( interpolator );
		resampler->SetDefaultPixelValue( 255 ); // highlight regions without source				
		resampler->SetOutputOrigin( inputImage->GetOrigin() );
		resampler->SetOutputDirection(inputImage->GetDirection() );
									
		  InpImageType::SizeType   inputSize = 
                    inputImage->GetLargestPossibleRegion().GetSize();
  
		typedef InpImageType::SizeType::SizeValueType SizeValueType;

		const double dx = inputSize[0] * inputSpacing[0] / isoSpacing;
		const double dy = inputSize[1] * inputSpacing[1] / isoSpacing;
		const double dz = (inputSize[2] - 1 ) * inputSpacing[2] / isoSpacing;
		InpImageType::SizeType   size;

		size[0] = static_cast<SizeValueType>( dx );
		size[1] = static_cast<SizeValueType>( dy );
		size[2] = static_cast<SizeValueType>( dz );	
																											
		resampler->SetSize( size );																																																																											
		resampler->SetInput( smootherY->GetOutput() );
		
		std::cout << "Resampling the image.........: " << std::endl;
		
		resampler->Update();
		
		std::vector<unsigned char> color(3,255);						
		ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
		filtImg->AppendChannelFromData3D(resampler->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	    myImg = filtImg;
		return myImg;		
	
	}


ftk::Image::Pointer ftkPreprocess::CurvAnisotropicDiffusion(void)
{
		 
		 typedef itk::CurvatureAnisotropicDiffusionImageFilter<InpImageType,FloatImageType> FilterType;
		 typedef itk::CastImageFilter<FloatImageType,InpImageType> ReverseCastFilter;	
		 FilterType::Pointer cfilter = FilterType::New();
         cfilter->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));
		 cfilter->SetTimeStep(filterParams[0]);
		 cfilter->SetNumberOfIterations(filterParams[1]);
		 cfilter->SetConductanceParameter(filterParams[2]);
		
		std::cout << "Applying Curvature Anisotropic Diffusion Filter.........: " << std::endl;
		  
		try
		{
			cfilter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
		std::cerr << "Exception caught: " << err << std::endl;
		}
 
		FloatImageType::Pointer imf = cfilter->GetOutput();
		typedef itk::ImageRegionIterator<FloatImageType> IRI;
		IRI iter(imf,imf->GetLargestPossibleRegion());
						
		for(iter.GoToBegin();!iter.IsAtEnd(); ++iter)
			{
				float value = iter.Get();
				value = ((value < 0)?0:((value>255)?255:value));
				iter.Set(value);
			}

		ReverseCastFilter::Pointer rfilter = ReverseCastFilter::New();
		rfilter->SetInput(cfilter->GetOutput());
		rfilter->Update();
		
		std::vector<unsigned char> color(3,255);
		ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
		filtImg->AppendChannelFromData3D(rfilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
		myImg = filtImg;
		return myImg;			
										
																		
			/* typedef itk::ImageFileWriter<imagetype > WriterType;
			   WriterType::Pointer writer = WriterType::New();
			   std::string fName = this->FN+"_median"+".tif";
			   writer->SetFileName(fName.c_str());
			   writer->SetInput(mfilter->GetOutput());
			   writer->Update();*/
		 
}





ftk::Image::Pointer ftkPreprocess::OpeningbyReconstruction(void)
{
	typedef itk::BinaryBallStructuringElement<InpPixelType,3> StructuringElementType;
	typedef itk::OpeningByReconstructionImageFilter<InpImageType,InpImageType,StructuringElementType > ROpenFilter;

	ROpenFilter::Pointer  RgrayscaleOpen  = ROpenFilter::New();
	StructuringElementType  structuringElement;
	structuringElement.CreateStructuringElement();
	structuringElement.SetRadius(filterParams[0]);
	RgrayscaleOpen->SetKernel(structuringElement );	 
	RgrayscaleOpen->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));
	if(filterParams[1]){
	RgrayscaleOpen->PreserveIntensitiesOn();
					   }
	else{
	RgrayscaleOpen->PreserveIntensitiesOff();	
					   
		}				   
					   
					   
	if(filterParams[2]){
		RgrayscaleOpen->FullyConnectedOn ();
					   }
	else{
		RgrayscaleOpen->FullyConnectedOff ();	
		}
					   
					   
	std::cout << "Applying Grayscale Opening by Reconstruction Filter.........: " << std::endl;

	try
		{
			RgrayscaleOpen->Update();
		}
	catch( itk::ExceptionObject & err )
		{
			std::cerr << "Exception caught: " << err << std::endl;
		}

	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(RgrayscaleOpen->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;

    }



ftk::Image::Pointer ftkPreprocess::ClosingbyReconstruction(void)
{
	typedef itk::BinaryBallStructuringElement<InpPixelType,3> StructuringElementType;
	typedef itk::ClosingByReconstructionImageFilter<InpImageType,InpImageType,StructuringElementType > RCloseFilter;

	RCloseFilter::Pointer  RgrayscaleClose  = RCloseFilter::New();
	StructuringElementType  structuringElement;
	structuringElement.CreateStructuringElement();
	structuringElement.SetRadius(filterParams[0]);
	RgrayscaleClose->SetKernel(structuringElement );	 
	RgrayscaleClose->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));
	if(filterParams[1]){
	RgrayscaleClose->PreserveIntensitiesOn();
					   }
	else{
	RgrayscaleClose->PreserveIntensitiesOff();				   
		}				   
					   
					   
	if(filterParams[2]){
		RgrayscaleClose->FullyConnectedOn ();
					   }
	else{
		RgrayscaleClose->FullyConnectedOff ();	
		}
					   
					   
	std::cout << "Applying Grayscale Closing by Reconstruction Filter.........: " << std::endl;

	try
		{
			RgrayscaleClose->Update();
		}
	catch( itk::ExceptionObject & err )
		{
			std::cerr << "Exception caught: " << err << std::endl;
		}

	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(RgrayscaleClose->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;

    }


ftk::Image::Pointer ftkPreprocess::MeanFilter(void)
{
	 
	 typedef itk::MeanImageFilter<InpImageType,InpImageType> FilterType;     
	 FilterType::Pointer mFilter = FilterType::New();
	 InpImageType::SizeType indexRadius; 
	 indexRadius[0] = filterParams[0]; // radius along x 
	 indexRadius[1] = filterParams[1]; // radius along y 
	 indexRadius[2] = filterParams[2]; // radius along y 

	mFilter->SetRadius( indexRadius ); 
	mFilter->SetInput(myImg->GetItkPtr<InpPixelType>(0,0) );
	
    std::cout << "Applying Mean Filter.........: " << std::endl;
	
	try
    {
		mFilter->Update();
    }
	catch( itk::ExceptionObject & err )
    {
		std::cerr << "Exception caught: " << err << std::endl;
    }

	/* typedef itk::ImageFileWriter<imagetype > WriterType;
	WriterType::Pointer writer = WriterType::New();
	std::string fName = this->FN+"_median"+".tif";
	writer->SetFileName(fName.c_str());
	writer->SetInput(mfilter->GetOutput());
	writer->Update();
	 */	 

	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(mFilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;
}




ftk::Image::Pointer ftkPreprocess::LaplacianFilter(void)
{
   

  typedef unsigned char    CharPixelType;  //IO
  typedef double          RealPixelType;  //Operations

  const    unsigned int    Dimension = 2;

  typedef itk::Image<CharPixelType, Dimension>    CharImageType;
  typedef itk::Image<RealPixelType, Dimension>    RealImageType;
  
  typedef itk::ImageFileReader< CharImageType >  ReaderType;
  typedef itk::ImageFileWriter< CharImageType >  WriterType;

  typedef itk::CastImageFilter<InpImageType, DoubleImageType> CastToDoubleFilterType;
  typedef itk::CastImageFilter<DoubleImageType, InpImageType> CastToCharFilterType;

  typedef itk::RescaleIntensityImageFilter<DoubleImageType, DoubleImageType> RescaleFilter;
  
  typedef itk::LaplacianImageFilter< 
                              DoubleImageType, 
                              DoubleImageType >    LaplacianFilter;

  typedef itk::ZeroCrossingImageFilter<
                              DoubleImageType, 
                              DoubleImageType>     ZeroCrossingFilter;

  CastToDoubleFilterType::Pointer toDouble = CastToDoubleFilterType::New();
  CastToCharFilterType::Pointer toChar = CastToCharFilterType::New();
  RescaleFilter::Pointer rescale = RescaleFilter::New();

  //Setting the ITK pipeline filter
  
  LaplacianFilter::Pointer lapFilter = LaplacianFilter::New();
  ZeroCrossingFilter::Pointer zeroFilter = ZeroCrossingFilter::New();  
  
  //The output of an edge filter is 0 or 1
  rescale->SetOutputMinimum(   0 );
  rescale->SetOutputMaximum( 255 );

  toDouble->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));
  toChar->SetInput( rescale->GetOutput() );
  //writer->SetInput( toChar->GetOutput() );

  //Edge Detection by Laplacian Image Filter:

  lapFilter->SetInput( toDouble->GetOutput() );
  zeroFilter->SetInput(lapFilter->GetOutput() );
  rescale->SetInput(zeroFilter->GetOutput() );

  try
    {
    rescale->Update();
    }
  catch( itk::ExceptionObject & err )
    { 
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
	} 
	
	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(rescale->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;
}



ftk::Image::Pointer ftkPreprocess::ThreeDSmoothingRGFilter(void)
{
	 
	 typedef itk::RecursiveGaussianImageFilter<FloatImageType,FloatImageType> RGFilterType;     
	 RGFilterType::Pointer RGFilter = RGFilterType::New();
	 
	 RGFilterType::Pointer filterX = RGFilterType::New(); 
	 RGFilterType::Pointer filterY = RGFilterType::New();
	 RGFilterType::Pointer filterZ = RGFilterType::New();
	 
	 filterX->SetDirection( 0 ); // 0 --> X direction 
	 filterY->SetDirection( 1 ); // 1 --> Y direction 
	 filterZ->SetDirection( 2 ); // 2 --> Z direction 
	 
	 
	 filterX->SetOrder( RGFilterType::ZeroOrder ); 
	 filterY->SetOrder( RGFilterType::ZeroOrder ); 
	 filterZ->SetOrder( RGFilterType::ZeroOrder ); 
	 
	 
	 filterX->SetInput(myImg->GetItkPtr<FloatPixelType>(0,0)); 
	 filterY->SetInput( filterX->GetOutput() ); 
	 filterZ->SetInput( filterY->GetOutput() );	
	 
	 
	 filterX->SetSigma( filterParams[0] ); 
	 filterY->SetSigma( filterParams[0] ); 
	 filterZ->SetSigma( filterParams[0] ); 
 
    if(filterParams[2]){ 
	 filterX->SetNormalizeAcrossScale( true ); 
	 filterY->SetNormalizeAcrossScale( true ); 
	 filterZ->SetNormalizeAcrossScale( true ); 
						}
	else{
	
	 filterX->SetNormalizeAcrossScale( false ); 
	 filterY->SetNormalizeAcrossScale( false ); 
	 filterZ->SetNormalizeAcrossScale( false );	
	
		}

	 filterX->SetNumberOfThreads (filterParams[1]); 
	 filterY->SetNumberOfThreads (filterParams[1] ); 
	 filterZ->SetNumberOfThreads (filterParams[1] );	
 	 
	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(filterZ->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;
}



ftk::Image::Pointer ftkPreprocess::NormalizeImage(void)
{
	 
	 typedef itk::NormalizeImageFilter<InpImageType,FloatImageType> FilterType;
     
	 FilterType::Pointer nFilter = FilterType::New();
	 	
    std::cout << "Applying Normalizing Filter.........: " << std::endl;
	
	try
    {
		nFilter->Update();
    }
	catch( itk::ExceptionObject & err )
    {
		std::cerr << "Exception caught: " << err << std::endl;
    }

    

	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(nFilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;
}



ftk::Image::Pointer ftkPreprocess::ShiftScale(void)
{	 
	 typedef itk::ShiftScaleImageFilter<InpImageType,InpImageType> FilterType;
	 FilterType::Pointer SSFilter = FilterType::New();
	 std::cout<<filterParams[0]<<std::endl;
	 SSFilter->SetShift(filterParams[0]);
	 SSFilter->SetScale(filterParams[1]);
				
    std::cout << "Applying the Shift and Scale Filter.........: " << std::endl;
	
	try
    {
		SSFilter->Update();
    }
	catch( itk::ExceptionObject & err )
    {
		std::cerr << "Exception caught: " << err << std::endl;
    }

	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(SSFilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;
}




ftk::Image::Pointer ftkPreprocess::SobelEdgeDetection(void)
{
	 
	 typedef itk::SobelEdgeDetectionImageFilter<InpImageType,InpImageType> FilterType;
     
	 FilterType::Pointer SEFilter = FilterType::New();
	 SEFilter->SetNumberOfThreads(filterParams[0]);

    std::cout << "Applying the Sobel Edge Detetction Filter.........: " << std::endl;
	
	try
    {
		SEFilter->Update();
    }
	catch( itk::ExceptionObject & err )
    {
		std::cerr << "Exception caught: " << err << std::endl;
    }

    

	std::vector<unsigned char> color(3,255);
    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
	filtImg->AppendChannelFromData3D(SEFilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
	myImg = filtImg;
	return myImg;
	
}





ftk::Image::Pointer ftkPreprocess::CurvatureFlow(void)
{
		
	     typedef itk::CurvatureFlowImageFilter<InpImageType,FloatImageType> FilterType;
		 typedef itk::CastImageFilter<FloatImageType,InpImageType> ReverseCastFilter;	
		 FilterType::Pointer filter = FilterType::New();
         filter->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));
		 filter->SetTimeStep(filterParams[0]);
		 filter->SetNumberOfIterations(filterParams[1]);
		 		
		std::cout << "Applying Curvature Flow Filter.........: " << std::endl;
		  
		try
		{
			filter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
		std::cerr << "Exception caught: " << err << std::endl;
		}
 
		 FloatImageType::Pointer imf = filter->GetOutput();
		 typedef itk::ImageRegionIterator<FloatImageType> IRI;
		 IRI iter(imf,imf->GetLargestPossibleRegion());
						
		 for(iter.GoToBegin();!iter.IsAtEnd(); ++iter)
			{
				float value = iter.Get();
				value = ((value < 0)?0:((value>255)?255:value));
				iter.Set(value);
			}

			ReverseCastFilter::Pointer rfilter = ReverseCastFilter::New();
			rfilter->SetInput(filter->GetOutput());
			rfilter->Update();
		
		   std::vector<unsigned char> color(3,255);
		   ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
		   filtImg->AppendChannelFromData3D(rfilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
		   myImg = filtImg;			
		   return myImg;								
																	
			/* typedef itk::ImageFileWriter<imagetype > WriterType;
			   WriterType::Pointer writer = WriterType::New();
			   std::string fName = this->FN+"_median"+".tif";
			   writer->SetFileName(fName.c_str());
			   writer->SetInput(mfilter->GetOutput());
			   writer->Update();
	 */	 
			//QApplication::restoreOverrideCursor();			 
	}


ftk::Image::Pointer ftkPreprocess::MinMaxCurvatureFlow(void)
{
		
	     typedef itk::MinMaxCurvatureFlowImageFilter<InpImageType,FloatImageType> FilterType;
		 typedef itk::CastImageFilter<FloatImageType,InpImageType> ReverseCastFilter;	
		 FilterType::Pointer filter = FilterType::New();
         filter->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));
		 filter->SetTimeStep(filterParams[0]);
		 filter->SetNumberOfIterations(filterParams[1]);
		 filter->SetStencilRadius(filterParams[2]);
		 		
		std::cout << "Applying Curvature Flow Filter.........: " << std::endl;
		  
		try
		{
			filter->Update();
		}
		catch( itk::ExceptionObject & err )
		{
		std::cerr << "Exception caught: " << err << std::endl;
		}
 
		 FloatImageType::Pointer imf = filter->GetOutput();
		 typedef itk::ImageRegionIterator<FloatImageType> IRI;
		 IRI iter(imf,imf->GetLargestPossibleRegion());
						
		 for(iter.GoToBegin();!iter.IsAtEnd(); ++iter)
			{
				float value = iter.Get();
				value = ((value < 0)?0:((value>255)?255:value));
				iter.Set(value);
			}

			ReverseCastFilter::Pointer rfilter = ReverseCastFilter::New();
			rfilter->SetInput(filter->GetOutput());
			rfilter->Update();
		
		   std::vector<unsigned char> color(3,255);
		   ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
		   filtImg->AppendChannelFromData3D(rfilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
		   myImg = filtImg;			
		   return myImg;								
																	
			/* typedef itk::ImageFileWriter<imagetype > WriterType;
			   WriterType::Pointer writer = WriterType::New();
			   std::string fName = this->FN+"_median"+".tif";
			   writer->SetFileName(fName.c_str());
			   writer->SetInput(mfilter->GetOutput());
			   writer->Update();
	 */	 
			//QApplication::restoreOverrideCursor();			 
	}
	
	
	
	
	
ftk::Image::Pointer ftkPreprocess::VesselFilter(void)
{
	    // typedef itk::AnisotropicDiffusionVesselEnhancementFilter<InpImageType,FloatImageType> FilterType;
		 /* typedef itk::CastImageFilter<FloatImageType,InpImageType> ReverseCastFilter;	
				  FilterType::Pointer filter = FilterType::New();
				  filter->SetInput(myImg->GetItkPtr<InpPixelType>(0,0));
				  filter->SetSigmaMin(filterParams[0]);
				  filter->SetSigmaMax(filterParams[1]);
				  filter->SetNumberOfSigmaSteps(filterParams[2]);
				  filter->SetNumberOfIterations(filterParams[3]);

						 
				 std::cout << "Applying Anisotropic Diffusion Vessel Enhancement Filter.........: " << std::endl;
				   
				 try
				 {
					 filter->Update();
				 }
				 catch( itk::ExceptionObject & err )
				 {
				 std::cerr << "Exception caught: " << err << std::endl;
				 }
		  
				  FloatImageType::Pointer imf = filter->GetOutput();
				  typedef itk::ImageRegionIterator<FloatImageType> IRI;
				  IRI iter(imf,imf->GetLargestPossibleRegion());
								 
				  for(iter.GoToBegin();!iter.IsAtEnd(); ++iter)
					 {
						 float value = iter.Get();
						 value = ((value < 0)?0:((value>255)?255:value));
						 iter.Set(value);
					 }

					 ReverseCastFilter::Pointer rfilter = ReverseCastFilter::New();
					 rfilter->SetInput(filter->GetOutput());
					 rfilter->Update();
				 
				    std::vector<unsigned char> color(3,255);
				    ftk::Image::Pointer filtImg = ftk::Image::New();	 	 
				    filtImg->AppendChannelFromData3D(rfilter->GetOutput()->GetBufferPointer(), itk::ImageIOBase::UCHAR, sizeof(unsigned char), myImg->GetImageInfo()->numColumns, myImg->GetImageInfo()->numRows, myImg->GetImageInfo()->numZSlices, "gray", color, true);
				    myImg = filtImg;			
		  */
return myImg;
}
	
















