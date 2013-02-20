#include "itkRescaleIntensityImageFilter.h"
#include "SpineRing.h"
#include "ImageProc.h"

SpeedImageType::Pointer SpCandidate::GetSpeedFunc(SpeedImageType::Pointer Spregion )
{

	SpeedImageType::Pointer SpeedImage;
	
	//  The RescaleIntensityImageFilter type is declared below. This filter will
	//  renormalize image to 0.1 -> 1
	//
	typedef itk::RescaleIntensityImageFilter< 
		SpeedImageType, 
		SpeedImageType >   CastFilterType;

	CastFilterType::Pointer caster = CastFilterType::New();
	caster->SetOutputMinimum(   0.1 );
	caster->SetOutputMaximum( 1 );
	caster->SetInput( Spregion );
	SpeedImage = caster->GetOutput();
	caster->Update();
	return SpeedImage;
}

void SpCandidate::RegionToPath()
{

	//   stpt{1}  = brtpt; ok
	//endpt{1} = bmu;  never used
	//stpt{2}  = X;   never used
	//endpt{2} = mucand;//all cand mu pt coords  ok

	//??chkp{2}=0;



	//[x1 y1 z1 spregion startpt endpts options] =   regionmap(pxls, img, stpt{1}, endpt{2}, maxiter);
	//uses pxls, NbrMuPtsCont
	PointSetContainerType::Pointer  EndPtsCont  = TranslatePts(NbrMusPtCont,  -1);
	PointType                       Startpt     = TranslatePts(brtpt,  -1);

	SpeedImageType::Pointer         SpeedReg = GetSpeedSubIm();
	DEBUGSTMT(ImStatsSpeedImageType(SpeedReg, "speed Reg"));

	// I have implemeneted a few vairiants for computing the speed
	// however, the best one has been the simplest grayscale values
	// Reason: size of the spine makes any modification to 
	//  the intensities potentially corrupting.
	// SpeedImageType::Pointer         SpeedFunc   = GetSpeedFunc(SpeedReg, false, false, false, 2, -.3, 2);
	SpeedImageType::Pointer         SpeedFunc   = GetSpeedFunc(SpeedReg);
	std::stringstream tmpstr0, tmpstr;
	tmpstr0 << ParentDet->imdbg->basefilename << candID << "_" << traceID ;

	RGBImageType::Pointer rgbreg2D;//leave declaration even if not debugging

	if (spr_DEBUGGING >= path_DBG) 
	{
		ImStatsSpeedImageType(SpeedFunc, "my speed func");

		SpineImageType::IndexType spIndex;
		for (int ii=0; ii<3; ii++)
			spIndex[ii] = brtpt[ii];
		std::cout << "orig img brtpt Val=" << ParentDet->image->GetPixel(spIndex) << std::endl;

		for (int ii=0; ii<3; ii++)
			spIndex[ii] = Startpt[ii];
		std::cout << "start reg Val=" << SpeedReg->GetPixel(spIndex) << std::endl;

		SpeedImageType::IndexType spdIndex;
		for (int ii=0; ii<3; ii++)
			spdIndex[ii] = Startpt[ii];
		std::cout << "start speed Val=" << SpeedFunc->GetPixel(spdIndex) << std::endl;

		//SpineImageType::Pointer         SpeedFun    = GetSpeedFunction(spregion);

		


		typedef itk::RescaleIntensityImageFilter< 
			SpeedImageType, 
			SpineImageType >   CastFilterType;
		CastFilterType::Pointer caster = CastFilterType::New();
		// only casting to proper type. pixel values are the same
		caster->SetOutputMinimum(   0 );
		caster->SetOutputMaximum( 255 );
		caster->SetInput( SpeedReg );
		SpineImageType::Pointer spinereg = caster->GetOutput();
		caster->Update();

		tmpstr << tmpstr0.str() << "_CandReg3D.tif" ;

		writeImageToFile(spinereg, tmpstr.str());
		//SpeedImage2DType::Pointer Sp2DReg;
		SpineImage2DType::Pointer Sp2DReg;
		//Sp2DReg = MIPGenerator<SpeedImageType, SpeedImage2DType>(SpeedFun);
		Sp2DReg = MaxIPGenerator( spinereg);//	SpeedReg);
		rgbreg2D = ParentDet->imdbg->CreateRGBImFrom2D(Sp2DReg);
		RGBPixelType color;
		color.SetRed(255);
		color.SetGreen(0);
		color.SetBlue(0);
		ParentDet->imdbg->PlotPtAsBoxOnRGBImage(rgbreg2D, color, Startpt, 1, 1);
		color.SetRed(100);
		color.SetGreen(240);
		color.SetBlue(60);
		ParentDet->imdbg->PlotContAsBoxOnRGBImage(rgbreg2D, color, EndPtsCont, 1, 1);
		//writeImageToFile<SpeedImage2DType>(Sp2DReg, tmpstr.str());
		//writeImageToFile(Sp2DReg, tmpstr.str());
		tmpstr.str(std::string());
		tmpstr << tmpstr0.str() << "_CandRegRGB.png" ;
		writeImageToFile(rgbreg2D, tmpstr.str());

		tmpstr.str(std::string());
		tmpstr << tmpstr0.str() << "_SpeedFunc.mhd" ;
		writeImageToFile(SpeedFunc, tmpstr.str());
	}
#if 0		// READING TEST speedfun:
	typedef itk::ImageFileWriter< SpeedImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( "myspeedfun.mhd" );
	writer->SetInput( SpeedFun);
	writer->Update();

	typedef itk::ImageFileReader< SpeedImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( "Synthetic-04-Speed.mhd" );
	try
	{
		reader->Update();
		SpeedFunc = reader->GetOutput();
		SpeedFunc->DisconnectPipeline();
		Startpt[0]=114;Startpt[1]=66; Startpt[2]=11;
		ImStatsSpeedImageType(SpeedFunc, "SynthSpeed");

	}
	catch(itk::ExceptionObject &err)
	{
		std::cerr << "ExceptionObject caught!" <<std::endl;
		std::cerr << err << std::endl;
	}
#endif	
	tmpstr.str(std::string());
	tmpstr << tmpstr0.str() << "_Path3D.tif"  ;
	Path = SpeedToFMPath(SpeedFunc, Startpt, EndPtsCont, tmpstr.str(), rgbreg2D);
	///////////////////////////////////////////////////////////
	/// testing SpSubIm with caster only, no grad:
	//////////////////////////////////////////////////////////
	//typedef itk::RescaleIntensityImageFilter< 
	//	SpeedImageType, 
	//	SpeedImageType >   CastFilterType2;
	//CastFilterType2::Pointer caster2 = CastFilterType2::New();
	//// only casting to proper type. pixel values are the same
	//caster2->SetOutputMinimum(   0.001 );
	//caster2->SetOutputMaximum( 1 );
	//caster2->SetInput( SpeedReg );
	//SpeedImageType::Pointer speedfunc2 = caster2->GetOutput();
	//caster2->Update();
	//tmpstr.str(std::string());
	//Path = FMPath(speedfunc2, Startpt, EndPtsCont, tmpstr.str(), rgbim2D);
}
