#include<iostream>
#include<string>
#include<sstream>


#include "itkImageFileWriter.h"
#include "itkStatisticsImageFilter.h"
#include "itkMinimumProjectionImageFilter.h"
#include "itkMaximumProjectionImageFilter.h"
#include "ImageProc.h"
#include "vnl/vnl_sample.h"

typedef itk::StatisticsImageFilter<SpineImageType>  StatisticsType;
typedef itk::ImageRegionIterator<SpineImageType> IteratorType;

typedef itk::ImageFileWriter< RGBImageType >  RGBWriterType;


typedef itk::ImageFileWriter< SpineImage2DType > SpineImage2DWriterType;
typedef itk::ImageFileWriter< BinaryImage2DType > BinaryImage2DWriterType;

//template <typename T> 
//void ImStats(typename T::Pointer im)
//{
//	typedef typename itk::StatisticsImageFilter<T>  StatisticsType;
//	StatisticsType::Pointer statistics = StatisticsType::New();
//	statistics->SetInput(im );
//	statistics->Update();
//	float imin = statistics->GetMinimum();
//	float imax = statistics->GetMaximum();
//	float imean = statistics->GetMean();
//	float istd = statistics->GetSigma();
//
//	std::cout << "Input Image Statistics: min:"<< imin << " max:" << imax << " mean:"<< imean << " std:" << istd << std::endl;
//}


void ImStatsSpeedImageType(SpeedImageType::Pointer im, std::string s)
{
	typedef itk::StatisticsImageFilter<SpeedImageType>  StatisticsType;
	StatisticsType::Pointer statistics = StatisticsType::New();
	statistics->SetInput(im );
	statistics->Update();
	float imin = statistics->GetMinimum();
	float imax = statistics->GetMaximum();
	float imean = statistics->GetMean();
	float istd = statistics->GetSigma();

	std::cout << s << " Image Statistics: min:"<< imin << " max:" << imax << " mean:"<< imean << " std:" << istd << std::endl;
}

void ImageDebugger::ImageStatistics(SpineImageType::Pointer & im3D)	{

	StatisticsType::Pointer statistics = StatisticsType::New();
	statistics->SetInput(im3D );
	statistics->Update();
	float imin = statistics->GetMinimum();
	float imax = statistics->GetMaximum();
	float imean = statistics->GetMean();
	float istd = statistics->GetSigma();

	std::cout << "Input Image Statistics: min:"<< imin << " max:" << imax << " mean:"<< imean << " std:" << istd << std::endl;
	if ((imean - imin) < (imax - imean))	
	{
		std::cout << "Inverting image (becomes dark foreground over bright background) & rescaling to [0,255]." << std::endl;
		IteratorType it(im3D, im3D->GetRequestedRegion());
		float imax2 = imean+2*istd;
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)	
		{
			float d = it.Get();
			d = (d > (imax2)) ? (imax2) : d;
			d = 1 - (d - imin)/(imax2 - imin);
			it.Set(255.0*d);
		}
	}
	else 
	{
		std::cout << "No inverting necessary. Rescaling to [0,255]." << std::endl;
		IteratorType it(im3D, im3D->GetRequestedRegion());
		float imin2 = imean-2*istd;
		for (it.GoToBegin(); !it.IsAtEnd(); ++it)	
		{
			float d = it.Get();
			d = (d < (imin2)) ? (imin2) : d;
			d = (d - imin2)/(imax - imin2);
			it.Set(255.0*d);
		}
	}
}

////void ImageDebugger::MIPGenerator()	{
////	if (spr_DEBUGGING < im_DBG && !NeedMIP)
////		return;
////	std::string MIPfilename = basefilename + std::string("Zproj.tif");
////	typedef itk::MinimumProjectionImageFilter<SpineImageType, SpineImage2DType> MIPfilterType;
////	MIPfilterType::Pointer MIPfilter = MIPfilterType::New();
////	MIPfilter->SetInput(im3D);
////	MIPfilter->SetProjectionDimension( 2 );
////	std::cout <<"MIP Generated!!" <<std::endl;
////
////	SpineImage2DWriterType::Pointer writer2D = SpineImage2DWriterType::New();
////	MIP = MIPfilter->GetOutput();
////	MIPfilter->Update();
////	writeImageToFile<SpineImage2DType>(MIP, MIPfilename);
////	////writer2D->SetInput(MIPfilter->GetOutput());
////	//////writer2D->SetInput(MIP);
////	////writer2D->SetFileName( MIPfilename );
////	////	
////	////try
////	////	{
////	////	writer2D->Update();
////	////	} 
////	////catch ( itk::ExceptionObject & excp )
////	////	{
////	////	std::cerr << excp << std::endl;
////	////	return;  // EXIT_FAILURE;
////	////	}
////}

void ImageDebugger::MaxIPGenerator(BinaryImageType::Pointer bin3D)	{
	// max proj for the binary CC image 
	if (spr_DEBUGGING < im_DBG)
		return;
	std::string MIPfilename = basefilename + std::string("Zccproj.tif");
	typedef itk::MaximumProjectionImageFilter<BinaryImageType, BinaryImage2DType> MIPfilterType;
	MIPfilterType::Pointer MIPfilter = MIPfilterType::New();
	MIPfilter->SetInput(bin3D);
	MIPfilter->SetProjectionDimension( 2 );
	std::cout <<"MIP Generated!!" <<std::endl;
	
	BinaryImage2DWriterType::Pointer writer2D = BinaryImage2DWriterType::New();
	bin2D = MIPfilter->GetOutput();
	MIPfilter->Update();
	writeImageToFile(bin2D, MIPfilename);
	////writer2D->SetInput(MIPfilter->GetOutput());
	////writer2D->SetFileName( MIPfilename );
	////	
	////try
	////	{
	////	writer2D->Update();
	////	} 
	////catch ( itk::ExceptionObject & excp )
	////	{
	////	std::cerr << excp << std::endl;
	////	return;  // EXIT_FAILURE;
	////	}
}

RGBImageType::Pointer ImageDebugger::CreateRGBImFrom2D(SpineImage2DType::Pointer im2d)
{
	RGBImageType::Pointer RGBImage = RGBImageType::New();
	RGBImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;

  //if (!MIP)
  //{
	 // NeedMIP = true;
	 // MIP = MIPGenerator(im3D);
  //}

  RGBImageType::SizeType size = im2d->GetRequestedRegion().GetSize();;

  RGBImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  RGBImage->SetRegions( region );
  RGBImage->Allocate();
  itk::ImageRegionIterator<RGBImageType> iterRGB(RGBImage,RGBImage->GetBufferedRegion());
  itk::ImageRegionConstIterator<SpineImage2DType> iter(im2d, im2d->GetBufferedRegion());

  RGBPixelType newPixel;
  for (iter.GoToBegin(), iterRGB.GoToBegin(); !iter.IsAtEnd(); ++iter, ++iterRGB)	{

		newPixel.SetRed(static_cast<unsigned char>(iter.Get()));
		newPixel.SetGreen(static_cast<unsigned char>(iter.Get()));
	    newPixel.SetBlue(static_cast<unsigned char>(iter.Get()));
	  	iterRGB.Set(newPixel);

  }

  return RGBImage;
}

////void ImageDebugger::Write2DRGBImageToFile(RGBImageType::Pointer rgbim2D, std::string fn)
////{
////	  RGBWriterType::Pointer writer = RGBWriterType::New();
////	  std::stringstream tmpstr;
////	  //tmpstr << basefilename << appendtofname ;
////	  writer->SetFileName( fn); //tmpstr.str() );
////	  writer->SetInput( rgbim2D );
////	   std::cout <<"writing rgbim file " << fn ;
////	  try
////	  {
////		  writer->Update();
////		  std::cout <<" ... done " << std::endl;
////	  }
////	  catch (itk::ExceptionObject & e)
////	  {
////		  std::cerr << "Exception in Writer: " << e << std::endl;
////		  exit(0);
////	  }
////}

void ImageDebugger::WriteTempObj2D(int trnum, int ringnum, IndexVecType *inring_f_idx)	{
	if (spr_DEBUGGING < ring_DBG)
		return;

  RGBImageType::Pointer RGBImage = CreateRGBImFrom2D(MIP);
  
  RGBPixelType pixelValue;
  RGBImageType::IndexType pixelIndex;



  std::cout << "		Writing ring "<<  ringnum << " pixels of trace " << trnum << " to MIP";
  IndexVecType::iterator                xiter;
  for (xiter  = inring_f_idx->begin(); xiter != inring_f_idx->end(); xiter++) {		
		pixelIndex[0] = static_cast<long int> ((*xiter)[0]);
		pixelIndex[1] = static_cast<long int> ((*xiter)[1]);

	    pixelValue[2] = 0;
		pixelValue[1] = 255;
		pixelValue[0] = 0;

		RGBImage->SetPixel(pixelIndex, pixelValue);
	}

  //RGBWriterType::Pointer writer = RGBWriterType::New();

  std::stringstream tmpstr;
  tmpstr << basefilename << "_ring_" << trnum << "_" << ringnum << ".png" ;
  writeImageToFile(RGBImage, tmpstr.str() );
  //writer->SetFileName( tmpstr.str() );
  //writer->SetInput( RGBImage );
  //try
  //  {
  //    writer->Update();
  //    std::cout <<"...done." << std::endl;
  //  }
  //catch (itk::ExceptionObject & e)
  //  {
  //    std::cerr << "Exception in Writer: " << e << std::endl;
  //    exit(0);
  //  }
}

RGBPixelType	GlobalDetector::GetColorMap(int loc)
{
	if (colormap.empty())
	{
		RGBPixelType px;
		
		// get max number of spine cands per dendrite
		GetSpCandStats();
		// Getting a colormap that is enough for all dendrites with their spine cands
		colormap.resize( maxmapsz + 1 );
		vnl_sample_reseed( 1031571 );
		for (unsigned short i=0; i < colormap.size(); ++i)
		{
			px.SetRed(
				static_cast<unsigned char>(255*vnl_sample_uniform( 0.3333, 1.0 ) ));
			px.SetGreen(
				static_cast<unsigned char>(255*vnl_sample_uniform( 0.3333, 1.0 ) ));
			px.SetBlue(
				static_cast<unsigned char>(255*vnl_sample_uniform( 0.3333, 1.0 ) ));
			colormap[i] = px;
		}
	}
	return colormap[loc];
}

void ImageDebugger::WriteSpCandidates2RGBImage(GlobalDetector *gd) 
{
  RGBImageType::Pointer RGBIm = CreateRGBImFrom2D(MIP);
	int validspc = gd->WriteDendSpCandMap2RGB(RGBIm);
  
	std::stringstream tmpstr;
	tmpstr << basefilename << "_totspcands_" << gd->GetSpCandTotal()<< "_valid_" << validspc << ".png" ;
    writeImageToFile(RGBIm,  tmpstr.str() );
//	Write2DRGBImageToFile(RGBIm, tmpstr.str() );
  //try
  //  {
		//RGBWriterType::Pointer writer = RGBWriterType::New();
		//writer->SetInput (RGBIm);
		//std::stringstream tmpstr;
		//tmpstr << basefilename << "_totspcands_" << spcandtotal<< "_valid_" << validspc << ".png" ;
		//writer->SetFileName( tmpstr.str() );
		//writer->Update();
  //  }
  //catch( itk::ExceptionObject & excep )
  //  {
		//std::cerr << "Exception caught !" << std::endl;
		//std::cerr << excep << std::endl;
  //  }
}



void ImageDebugger::PlotPtAsBoxOnRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color, RGBPointType pt, int plotboxx, int plotboxy)
{
	RGBIndexType rgbndx;

	for (rgbndx[1] = pt[1] - plotboxy; rgbndx[1] <= pt[1] + plotboxy; rgbndx[1]++)
		for (rgbndx[0] = pt[0] - plotboxx; rgbndx[0] <= pt[0] + plotboxx; rgbndx[0]++)
			if (rgbim2D->GetBufferedRegion().IsInside(rgbndx))
				rgbim2D->SetPixel(rgbndx, color);

}

void ImageDebugger::PlotPtAsBoxOnRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color, PointType pt, int plotboxx, int plotboxy)
{
	RGBPointType rgbpt;
	rgbpt[0] = pt[0];
	rgbpt[1] = pt[1];
	PlotPtAsBoxOnRGBImage(rgbim2D, color, rgbpt, plotboxx, plotboxy);
}

void ImageDebugger::PlotContAsBoxOnRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color, PointSetContainerType::Pointer pcont, int plotboxx, int plotboxy)
{
	PointSetType::PointsContainerIterator	pciter = pcont->Begin();
	PointType    pt;

	while (pciter != pcont->End()) 
	{
		pt = pciter->Value();
		PlotPtAsBoxOnRGBImage(rgbim2D, color, pt, plotboxx, plotboxy);
		pciter++;
	}
}


int ImageDebugger::PlotContOnRGBImage(RGBImageType::Pointer rgbim2D, RGBPixelType color, PointSetContainerType::Pointer pcont) 
{
	//if (!valid)
	//	return 0;

		PointSetType::PointsContainerIterator	pciter = pcont->Begin();
		PointSetType::PointType spp;
        RGBImageType::IndexType rgbndx;
		while (pciter != pcont->End()) 
		{
			spp = pciter->Value();
			rgbndx[0] = spp[0]; rgbndx[1] = spp[1];
			rgbim2D->SetPixel(rgbndx, color);
			pciter++;
		}
		
	return 1;
}

SpineImage2DType::Pointer MIPGenerator(SpineImageType::Pointer im)	
{	
	//std::string MIPfilename
	typedef itk::MinimumProjectionImageFilter<SpineImageType, SpineImage2DType> MIPfilterType;
	MIPfilterType::Pointer MIPfilter = MIPfilterType::New();
	MIPfilter->SetInput(im);
	MIPfilter->SetProjectionDimension( 2 );
	MIPfilter->Update();
	std::cout <<"MIP Generated!!" <<std::endl;
	SpineImage2DType::Pointer mipim = MIPfilter->GetOutput();
	mipim->DisconnectPipeline();
	//writeImageToFile<T0>(mipim, MIPfilename);
	return mipim;
}
 SpeedImage2DType::Pointer MaxIPGenerator(SpeedImageType::Pointer im)	
{	
	//std::string MIPfilename
	typedef itk::MaximumProjectionImageFilter<SpeedImageType, SpeedImage2DType> MIPfilterType;
	MIPfilterType::Pointer MIPfilter = MIPfilterType::New();
	MIPfilter->SetInput(im);
	MIPfilter->SetProjectionDimension( 2 );
	MIPfilter->Update();
	std::cout <<"MIP Generated!!" <<std::endl;
	SpeedImage2DType::Pointer mipim = MIPfilter->GetOutput();
	mipim->DisconnectPipeline();
	//writeImageToFile<T0>(mipim, MIPfilename);
	return mipim;
}
SpineImage2DType::Pointer MaxIPGenerator(SpineImageType::Pointer im)	
{	
	//std::string MIPfilename
	typedef itk::MaximumProjectionImageFilter<SpineImageType, SpineImage2DType> MIPfilterType;
	MIPfilterType::Pointer MIPfilter = MIPfilterType::New();
	MIPfilter->SetInput(im);
	MIPfilter->SetProjectionDimension( 2 );
	MIPfilter->Update();
	std::cout <<"MIP Generated!!" <<std::endl;
	SpineImage2DType::Pointer mipim = MIPfilter->GetOutput();
	//writeImageToFile<T0>(mipim, MIPfilename);
	mipim->DisconnectPipeline();
	return mipim;
}
int writeImageToFile(SpeedImage2DType::Pointer im, std::string fn)
{
	std::cout << "Writing  " << fn << std::endl;
        
	typedef itk::ImageFileWriter<SpeedImage2DType> WriterType;
    WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fn);
    writer->SetInput(im);
    try
    {
        writer->Update();
    }
    catch(itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught!" <<std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int writeImageToFile(SpineImage2DType::Pointer im, std::string fn)
{
	std::cout << "Writing  " << fn << std::endl;
        
	typedef itk::ImageFileWriter<SpineImage2DType> WriterType;
    WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fn);
    writer->SetInput(im);
    try
    {
        writer->Update();
    }
    catch(itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught!" <<std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int writeImageToFile(BinaryImage2DType::Pointer im, std::string fn)
{
	std::cout << "Writing  " << fn << std::endl;
        
	typedef itk::ImageFileWriter<BinaryImage2DType> WriterType;
    WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fn);
    writer->SetInput(im);
    try
    {
        writer->Update();
    }
    catch(itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught!" <<std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int writeImageToFile(RGBImageType::Pointer im, std::string fn)
{
	std::cout << "Writing  " << fn << std::endl;
        
	typedef itk::ImageFileWriter<RGBImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fn);
    writer->SetInput(im);
    try
    {
        writer->Update();
    }
    catch(itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught!" <<std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int writeImageToFile(SpineImageType::Pointer im, std::string fn)
{
	std::cout << "Writing  " << fn << std::endl;
        
	typedef itk::ImageFileWriter<SpineImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fn);
    writer->SetInput(im);
    try
    {
        writer->Update();
    }
    catch(itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught!" <<std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int writeImageToFile(SpeedImageType::Pointer im, std::string fn)
{
	std::cout << "Writing  " << fn << std::endl;
        
	typedef itk::ImageFileWriter<SpeedImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(fn);
    writer->SetInput(im);
    try
    {
        writer->Update();
    }
    catch(itk::ExceptionObject &err)
    {
        std::cerr << "ExceptionObject caught!" <<std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

////template <typename T, typename T0>
////typename T0::Pointer MIPGenerator(typename T::Pointer im)	{
////	
////	//std::string MIPfilename
////	typedef itk::MinimumProjectionImageFilter<T, T0> MIPfilterType;
////	typename MIPfilterType::Pointer MIPfilter = MIPfilterType::New();
////	MIPfilter->SetInput(im);
////	MIPfilter->SetProjectionDimension( 2 );
////	MIPfilter->Update();
////	std::cout <<"MIP Generated!!" <<std::endl;
////	T0::Pointer mipim = MIPfilter->GetOutput();
////	//writeImageToFile<T0>(mipim, MIPfilename);
////	return mipim;
////}
////
////template <typename T>
////int writeImageToFile(typename T::Pointer im, std::string fn)
////{
////	std::cout << "Writing  " << fn << std::endl;
////        
////	typedef typename itk::ImageFileWriter<T> WriterType;
////    typename WriterType::Pointer writer = WriterType::New();
////	writer->SetFileName(fn);
////    writer->SetInput(im);
////    try
////    {
////        writer->Update();
////    }
////    catch(itk::ExceptionObject &err)
////    {
////        std::cerr << "ExceptionObject caught!" <<std::endl;
////        std::cerr << err << std::endl;
////        return EXIT_FAILURE;
////    }
////    return EXIT_SUCCESS;
////}

#if 0
// manual MIPGenerator
// this code is slow
// resorting to itk filter above
	SpineImageType::SizeType sz3 = im3D->GetRequestedRegion().GetSize();
	SpineImage2DType::SizeType sz2 = {{sz3[0], sz3[1] }};

	MIP = SpineImage2DType::New(); 
	MIP->SetRegions(sz2 );
	MIP->Allocate();


	for (unsigned long x=0; x<sz3[0]; x++)	{
		for (unsigned long y=0; y<sz3[1]; y++)	{

			//ImageType3D::PixelType maxVal = 0.0;
			SpineImageType::PixelType minVal = 255.0;

			for (unsigned long z=0; z<sz3[2]; z++)	{
				SpineImageType::IndexType nd3 = {{x,y,z}};
				SpineImageType::PixelType val = im3D->GetPixel(nd3);
			//	maxVal = (maxVal > val) ? maxVal : val;
				minVal = (minVal < val) ? minVal : val;
			}

			SpineImage2DType::IndexType nd2 = {{x, y}};
			//im2D->SetPixel(nd2, maxVal);
			MIP->SetPixel(nd2, minVal);
		}
	}
#endif 

