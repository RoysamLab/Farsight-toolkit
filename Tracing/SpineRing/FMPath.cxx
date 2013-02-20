#if defined(_MSC_VER)
//Warning about: identifier was truncated to '255' characters in the debug information (MVC6.0 Debug)
#pragma warning( disable : 4786 )
#endif

// General includes
#include <string>
#include <iostream>

// ITK includes
#include "itkNumericTraits.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPolyLineParametricPath.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkArrivalFunctionToPathFilter.h"
//#include "itkSpeedFunctionToPathFilter.h"
#include "itkPathIterator.h"
#include "itkGradientDescentOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkIterateNeighborhoodOptimizer.h"
#include "itkStatisticsImageFilter.h"
//#include "itkFastMarchingUpwindGradientImageFilter.h"
#include "itkFastMarchingImageFilter.h"

#include "SpineRing.h"
#include "ImageProc.h"


PathType::Pointer SpCandidate::SpeedToFMPath(SpeedImageType::Pointer speed, 
									  PointType start, 
									  PointSetContainerType::Pointer endptscont, 
									  std::string outfn, 
									  RGBImageType::Pointer rgbim2d)
{
	//typedef float	PixelType;//speed
	typedef	unsigned char OutputPixelType;
	//typedef itk::Image< PixelType, spr_SPIMGDIMENSION	> ImageType;//speed
	typedef	itk::Image<	OutputPixelType, spr_SPIMGDIMENSION	> OutputImageType;
	//typedef itk::ImageFileReader<	ImageType >	ReaderType;
	typedef	itk::ImageFileWriter< OutputImageType >	WriterType;
	//typedef	itk::SpeedFunctionToPathFilter<	SpineImageType,	PathType > PathFilterType;
	//typedef	itk::SpeedFunctionToPathFilter<	SpeedImageType,	PathType > PathFilterType;
	typedef	itk::ArrivalFunctionToPathFilter<	SpeedImageType,	PathType > PathFilterType;
	typedef	PathFilterType::CostFunctionType::CoordRepType	CoordRepType;
	typedef	itk::PathIterator< OutputImageType,	PathType > PathIteratorType;

	//char* OutputFilename = "TestPath.tif";

	speed->DisconnectPipeline();

	///////////////////////////////////////////////////////////
	////////////////      F M Step        /////////////////////
	///////////////////////////////////////////////////////////
	typedef  itk::FastMarchingImageFilter< SpeedImageType, 
		SpeedImageType >    FastMarchingFilterType;
	FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
	fastMarching->SetInput( speed );
	typedef FastMarchingFilterType::NodeContainer           NodeContainer;
	typedef FastMarchingFilterType::NodeType                NodeType;
	NodeContainer::Pointer seeds = NodeContainer::New();
	//  Software Guide : EndCodeSnippet 
	SpeedImageType::IndexType  seedPosition;
	seedPosition[0] = start[0];
	seedPosition[1] = start[1];
	seedPosition[2] = start[2];
	NodeType node;
	const double seedValue = 0.0;
	node.SetValue( seedValue );
	node.SetIndex( seedPosition );
	seeds->Initialize();
	seeds->InsertElement( 0, node );
	fastMarching->SetTrialPoints(  seeds  );
	fastMarching->SetOutputSize( speed->GetBufferedRegion().GetSize() );
	/////////////////////////////////////////////
	// the stopping time might not be really needed
	//const double stoppingTime=200.0;
	//fastMarching->SetStoppingValue(  stoppingTime  );

	fastMarching->Update();

	//////////////////////////
	////////////////      F M End         /////////////////////
	///////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////
	////////////////      Get Best EndPoint    ////////////////
	///////////////////////////////////////////////////////////
	SpeedImageType::IndexType arrIndex;
	SpeedImageType::Pointer Arrival = fastMarching->GetOutput();
	float bestArrival = 10000.0;//some large value
	SpeedPixelType pixelVal;

	PointSetType::PointsContainerIterator	pciter = endptscont->Begin();
	PointType								pt, bestendpt;
	bestendpt = pciter->Value(); //any initial value
	while (pciter != endptscont->End()) 
	{
		pt = pciter->Value();
		bool isInside = Arrival->TransformPhysicalPointToIndex( pt, arrIndex );
		if ( isInside )
		{
			pixelVal = Arrival->GetPixel( arrIndex );
			if (pixelVal < bestArrival)
			{
				bestArrival = pixelVal;
				bestendpt = pt;
			}
		}
		pciter++;
	}
	////////////////    Best EndPoint Done     ////////////////
	///////////////////////////////////////////////////////////

	// Create interpolator
	typedef	itk::LinearInterpolateImageFunction<SpeedImageType, CoordRepType>
		InterpolatorType;
	InterpolatorType::Pointer interp =	InterpolatorType::New();

	// Create cost function
	PathFilterType::CostFunctionType::Pointer cost	=
		PathFilterType::CostFunctionType::New();
	cost->SetInterpolator( interp );

	// Create optimizer
	///////////////////////////////////////////////////////////////
	/*typedef	itk::GradientDescentOptimizer OptimizerType;
	OptimizerType::Pointer	optimizer =	OptimizerType::New();
	optimizer->SetNumberOfIterations( 3000 );
	optimizer->SetLearningRate(15.0);*/
	//////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////
	typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
	OptimizerType::Pointer optimizer	= OptimizerType::New();
	optimizer->SetNumberOfIterations(	3000 );//50000
	optimizer->SetMaximumStepLength( 1.0 );
	optimizer->SetMinimumStepLength( 0.5 );
	optimizer->SetRelaxationFactor( 0.5 ); //.999

	///////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////
	// Create IterateNeighborhoodOptimizer
	// typedef itk::IterateNeighborhoodOptimizer OptimizerType;
	//OptimizerType::Pointer optimizer = OptimizerType::New();
	//optimizer->MinimizeOn( );
	//optimizer->FullyConnectedOn( );
	//OptimizerType::NeighborhoodSizeType nbrsize( 3 );
	//float StepLengthFactor = .5;
	//for (unsigned int i=0; i<3; i++)
	//nbrsize[i] = speed->GetSpacing()[i]* StepLengthFactor;
	//optimizer->SetNeighborhoodSize( nbrsize );

	// Create path filter
	PathFilterType::Pointer pathFilter	= PathFilterType::New();
	pathFilter->SetInput( Arrival );
	pathFilter->SetCostFunction( cost );
	pathFilter->SetOptimizer( optimizer	);
	pathFilter->SetTerminationValue( 1 );
	pathFilter->SetPathEndPoint ( bestendpt);

	// Compute the path
	pathFilter->Update(	);
	PathType::Pointer path;
	// Rasterize path
	//for	(unsigned int i=0; i<pathFilter->GetNumberOfOutputs(); i++)
	// Get the path
	path	= pathFilter->GetOutput( 0 );
	// Check path is valid
	if ( path->GetVertexList()->Size() == 0	)
	{
		std::cout << "WARNING: Path	" <<  "	contains no	points!" <<	std::endl;
	}

	//////////////////////////////////////////////////////////
	/////////////// DEBUG STUFF
	//SpeedImageType::Pointer Arrival = pathFilter->GetArrivalFunc();
	if (spr_DEBUGGING >= path_DBG)
	{
		typedef	itk::ImageFileWriter< SpeedImageType >	WriterType2;
		WriterType2::Pointer writer2	= WriterType2::New();
		writer2->SetFileName( "Arrival.vtk"	);
		writer2->SetInput( Arrival );
		writer2->Update();
		typedef itk::StatisticsImageFilter<SpeedImageType>  StatisticsType;
		StatisticsType::Pointer statistics = StatisticsType::New();
		statistics->SetInput(Arrival );
		statistics->Update();
		float imin = statistics->GetMinimum();
		float imax = statistics->GetMaximum();
		float imean = statistics->GetMean();
		float istd = statistics->GetSigma();

		std::cout << "Arrival Image Statistics: min:"<< imin << " max:" << imax << " mean:"<< imean << " std:" << istd << std::endl;
		//SpeedImageType::IndexType arrIndex;
		for (int ii=0; ii<3; ii++)
			arrIndex[ii] = bestendpt[ii];
		std::cout << "End Arrival Val=" << Arrival->GetPixel(arrIndex) << std::endl;
		for (int ii=0; ii<3; ii++)
			arrIndex[ii] = start[ii];
		std::cout << "start Arrival Val=" << Arrival->GetPixel(arrIndex) << std::endl;
		for (int ii=0; ii<3; ii++)
			arrIndex[ii] = start[ii]+1;
		std::cout << "start+1 Arrival Val=" << Arrival->GetPixel(arrIndex) << std::endl;

		for (int ii=0; ii<3; ii++)
			arrIndex[ii] = bestendpt[ii];
		std::cout << "End speed Val=" << speed->GetPixel(arrIndex) << std::endl;
		for (int ii=0; ii<3; ii++)
			arrIndex[ii] = start[ii];
		std::cout << "start speed Val=" << speed->GetPixel(arrIndex) << std::endl;
		for (int ii=0; ii<3; ii++)
			arrIndex[ii] = start[ii]+1;
		std::cout << "start+1 speed Val=" << speed->GetPixel(arrIndex) << std::endl;
		//itk::ImageRegionIterator<SpeedImageType> arrit(Arrival, Arrival->GetRequestedRegion());
		//SpeedPixelType spdpxl;
		//for ( arrit.GoToBegin(); !arrit.IsAtEnd(); ++arrit)
		//{
		//	std::cout << arrit.Get(); //spregion = (spregion-m)/(M-m);  %0 to 1
		//}
		// Allocate	output image
	}
	/////////// END DEBUG STUFF    ///////////////////////////////
	//////////////////////////////////////////////////////////////
	VERBOSE2("Path Size = ", path->GetVertexList()->Size());
	if (path->GetVertexList()->Size() >= spr_MINPATHSIZE)
	{
		OutputImageType::Pointer output = OutputImageType::New();
		output->SetRegions(	speed->GetLargestPossibleRegion() );
		output->SetSpacing(	speed->GetSpacing()	);
		output->SetOrigin( speed->GetOrigin() );
		output->Allocate( );
		output->FillBuffer(	itk::NumericTraits<OutputPixelType>::Zero );
		//WriterType::Pointer writer0	= WriterType::New();
		//writer0->SetFileName( "TestPathIm.tif"	);
		//writer0->SetInput( output );
		//try
		//{
		//	writer0->Update();
		//}
		//catch(itk::ExceptionObject &excep)
		//{
		//    std::cerr << "ExceptionObject caught!" <<std::endl;
		//	std::cerr << excep << std::endl;
		//}
		RGBPixelType color;
		color.SetRed(200);
		color.SetGreen(70);
		color.SetBlue(150);

		//int pathcount=0;

		// Iterate path	and	convert	to image
		try
		{
			PathIteratorType it( output, path );
			//RGBIndexType rgbndx;
			RGBPointType rgbpt, rgbpt0;
			SpineImageType::PointType pt0, ppt;
			for	(it.GoToBegin(); !it.IsAtEnd();	++it)
			{
				it.Set(	itk::NumericTraits<OutputPixelType>::max() );
				OutputImageType::IndexType idx = it.GetIndex();
				rgbpt[0]=idx[0]; rgbpt[1]=idx[1];
				ppt[0]  = idx[0]; ppt[1] = idx[1]; ppt[2] = idx[2];
				pt0  = TranslatePts(ppt, 1);
				rgbpt0[0]=pt0[0]; rgbpt0[1]=pt0[1]; 
				DEBUGSTMT(ParentDet->imdbg->PlotPtAsBoxOnRGBImage(rgbim2d, color, rgbpt, 0, 0));
				ParentDet->imdbg->PlotPtAsBoxOnRGBImage(ParentDet->GetRGBPathsImage(), ParentDet->GetColorMap(candID) , rgbpt0, 0, 0);
				//rgbndx[0]=idx[0]; rgbndx[1]=idx[1];
				//if (rgbim2d->GetBufferedRegion().IsInside(rgbndx))
				//	rgbim2d->SetPixel(rgbndx, color);
				DEBUGSTMT(VERBOSE2("Path pt ", it.GetIndex()));
				//pathcount++;
			}
			validPath = true;
			//VERBOSE2("Path count = ", pathcount);
		}
		catch(itk::ExceptionObject &excep)
		{
			std::cerr << "ExceptionObject caught!" <<std::endl;
			std::cerr << excep << "Path Size= "<<path->GetVertexList()->Size()<<std::endl;
		}



		// Write output
		if (spr_DEBUGGING >= path_DBG)
		{
			DEBUGSTMT(writeImageToFile(rgbim2d, outfn+ ".png"));
			WriterType::Pointer writer	= WriterType::New();
			writer->SetFileName( outfn + ".tif");
			writer->SetInput( output );
			try
			{
				writer->Update();
			}
			catch(itk::ExceptionObject &excep)
			{
				std::cerr << "ExceptionObject caught!" <<std::endl;
				std::cerr << excep << std::endl;
			}
		}
	}
	return path;
}
