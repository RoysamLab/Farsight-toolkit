#include "Voting_Seg.h"

Voting_Seg::Voting_Seg()
{
}

Voting_Seg::~Voting_Seg()
{
}

void Voting_Seg::Initialize()
{
	this->DefultPara();
}

void Voting_Seg::SetInputImage(DImageType_2D::Pointer im)
{
	this->InputImage = im;
	this->imagesize = InputImage->GetBufferedRegion().GetSize();
}

void Voting_Seg::SetParameters(std::vector<Parameter> para)
{
	this->parameters.voting_rmin = (int)para[0].value;
	this->parameters.voting_rmax = (int)para[1].value;
	this->parameters.voting_radius = (int)para[2].value;
	this->parameters.voting_gradthresh = (double)para[3].value;
	this->parameters.voting_scale = (double)para[4].value;
	this->parameters.smooth_Iteration = (unsigned int)para[5].value;
	this->parameters.smooth_conductance = (double)para[6].value;
	this->parameters.gradient_sigma = (double)para[7].value;
	this->parameters.sigmoid_alpha = (double)para[8].value;
	this->parameters.sigmoid_beta = (double)para[9].value;
	this->parameters.fastmarching_timethreshold = (int)para[10].value;
	this->parameters.ac_curvature = (double)para[11].value;
	this->parameters.ac_advect = (double)para[12].value;
	this->parameters.ac_iteration = (int)para[13].value;
	this->parameters.ac_propagation = (double)para[14].value;
	this->parameters.ac_maximumRMSE = (double)para[15].value;
}

void Voting_Seg::Update()
{
	this->Binarize_2D();
	this->SeedDetection_2D();
	this->LevelSet_2D();
}

void Voting_Seg::GetOutputImage()
{}

void Voting_Seg::DefultPara()
{}

void Voting_Seg::Binarize_2D()
{
}

void Voting_Seg::SeedDetection_2D()
{
	int hmin = 1;
	int hmax = 25;//35;//77; //15 was ok, 25 just to test
	int radius = 0.7*hmax;//25;//10;//77;
	double min_grad = 0.1;
	//threshold (for picking seed points, not necesary for now)
	double scale = 1.5; // Scale for computing the gradient using Gradient
	ftkVoting voteMain;
	voteMain.setParams(parameters.voting_rmin,parameters.voting_rmax,parameters.voting_radius,parameters.voting_gradthresh,parameters.voting_scale);
	voteMain.setPrefix("temp/");
	voteMain.compute(InputImage);
	this->VotingImage = voteMain.GetOutput();
	//this->WriteImage( this->VotingImage, "votingimage.tif" );
	this->SeedExtraction();
}

void Voting_Seg::Watershed_2D()
{
	typedef itk::CurvatureFlowImageFilter< DImageType_2D, DImageType_2D >CurvatureFlowImageFilterType;
	CurvatureFlowImageFilterType::Pointer smoothing = CurvatureFlowImageFilterType::New();
	smoothing->SetInput( this->InputImage );
	smoothing->SetNumberOfIterations( 5 );
	smoothing->SetTimeStep( 0.125 );
	this->WriteImage( smoothing->GetOutput(), "smoothingimage.tif" );

	typedef itk::ConfidenceConnectedImageFilter<DImageType_2D, DImageType_2D>ConnectedFilterType;
	ConnectedFilterType::Pointer confidenceConnected = ConnectedFilterType::New();
	confidenceConnected->SetInput( smoothing->GetOutput() );
	confidenceConnected->SetMultiplier( /*2.5*/ 1.5);
	confidenceConnected->SetNumberOfIterations( 10 );
	confidenceConnected->SetReplaceValue( 255 );
	confidenceConnected->SetInitialNeighborhoodRadius( 2 );
	confidenceConnected->ClearSeeds();
	for(int i = 0; i < this->seedVector.size(); i++)
	{
		confidenceConnected->AddSeed( this->seedVector[i] );
		//std::cout<<"seed "<<i<<" inddex[0] "<<this->seedVector[i][0]<<"and index[1] "<<this->seedVector[i][1]<<std::endl;
	}

	try
	{
		confidenceConnected->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	this->SegmentedImage = confidenceConnected->GetOutput();
	//this->WriteImage( this->SegmentedImage, "segmentedimage.tif" );
}

void Voting_Seg::WriteImage(DImageType_2D::Pointer image, const char* filename)
{
	typedef itk::RescaleIntensityImageFilter< DImageType_2D, USImageType_2D > RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput( image );
	rescaleFilter->SetOutputMaximum( 255 );
	rescaleFilter->SetOutputMinimum( 0 ); 
	try{
			rescaleFilter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
			std::cerr << "ExceptionObject caught in rescaling "<<filename<< " !" << std::endl;
			std::cerr << err << std::endl;
	}

	USWriterType::Pointer uswriter = USWriterType::New();
	uswriter->SetInput( rescaleFilter->GetOutput() );
	uswriter->SetFileName( filename );
	try{
			uswriter->Update();
	}
	catch( itk::ExceptionObject & err )
	{
			std::cerr << "ExceptionObject caught in writting "<<filename<< " !" << std::endl;
			std::cerr << err << std::endl;
	}
}

void Voting_Seg::SeedExtraction()
{	
	UCImageType_2D::IndexType start;
	start[0] = 0;
	start[1] = 0;
	UCImageType_2D::RegionType region;
	region.SetSize( this->imagesize );
	region.SetIndex( start );
	const UCImageType_2D::PixelType zeros = 0;
	

	UCImageType_2D::Pointer votingBinImage;
	votingBinImage= UCImageType_2D::New();
	votingBinImage->SetRegions( region );
	votingBinImage->Allocate();
	votingBinImage->FillBuffer( zeros );
	votingBinImage->Update();

	this->SeedImage = UCImageType_2D::New();
	this->SeedImage->SetRegions( region );	
	this->SeedImage->Allocate();		
	this->SeedImage->FillBuffer( zeros );	
	this->SeedImage->Update();

	this->SeedLabelImage = UCImageType_2D::New();
	this->SeedLabelImage->SetRegions( region );	
	this->SeedLabelImage->Allocate();		
	this->SeedLabelImage->FillBuffer( zeros );	
	this->SeedLabelImage->Update();

	DIteratorType_2D votingIterator(this->VotingImage, this->VotingImage->GetRequestedRegion());
	UCIteratorType_2D votingBinIterator(votingBinImage, votingBinImage->GetRequestedRegion());
	DIteratorType_2D::PixelType pixelvalue;
	const double threshold = 15000.0;
	for(long int i = 0; i < imagesize[0] * imagesize[1]; i++)
	{
		pixelvalue = votingIterator.Get();
		if (pixelvalue > threshold)
		{
			votingBinIterator.Set(255);
		}
		votingIterator++;
		votingBinIterator++;
	}

	const char* filename = "votingBinimage.tif";
	//int temp = writeImage< UCImageType_2D, UCImageType_2D >( votingBinImage, filename);

	typedef itk::ConnectedComponentImageFilter <UCImageType_2D, UCImageType_2D >ConnectedComponentImageFilterType;
	ConnectedComponentImageFilterType::Pointer connectFilter  = ConnectedComponentImageFilterType::New ();
	connectFilter->SetInput(votingBinImage);
	connectFilter->Update();
	filename = "connectedimage.tif";
	//temp = writeImage< UCImageType_2D, UCImageType_2D >( connectFilter->GetOutput(),filename  );

	typedef itk::LabelGeometryImageFilter< UCImageType_2D > LabelGeometryImageFilterType;
	LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	labelGeometryImageFilter->SetInput( connectFilter->GetOutput() );
	labelGeometryImageFilter->Update();
	LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
	LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
	
	int k = 1;
	this->seedNum = 0;
	this->seedVector.clear();
	for( allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++ )
	{
		LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
		DImageType_2D::PointType point;
		point = labelGeometryImageFilter->GetCentroid(labelValue);
		DImageType_2D::IndexType pixelindex;
		pixelindex[0]= point[0];
		pixelindex[1]= point[1];
		this->SeedImage->SetPixel(pixelindex,255);
		this->SeedLabelImage->SetPixel(pixelindex,k++);
		this->seedNum++;

		itk::Index<Dimension_2D> index;
		index[0] = point[0];
		index[1] = point[1];
		this->seedVector.push_back(index);
	}
	filename =  "seedimage.tif";
	int temp = writeImage< UCImageType_2D, UCImageType_2D >( this->SeedImage, filename );
}

void Voting_Seg::LevelSet_2D() 																	  
{
	const char* filename;
	int temp;

	std::cout<< "smoothing ....."<<std::endl;
	unsigned int smoothIteration = 10;
	double conductance = 9.0;
	typedef itk::CurvatureAnisotropicDiffusionImageFilter<DImageType_2D, DImageType_2D> SmoothingFilterType;
	SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
	smoothing->SetInput(this->InputImage);
	smoothing->SetTimeStep( 0.0625);
	smoothing->SetNumberOfIterations( this->parameters.smooth_Iteration);
	smoothing->SetConductanceParameter( this->parameters.smooth_conductance);
	smoothing->Update();
	filename = "smoothimage.nrrd";
	temp = writeImage< DImageType_2D >( smoothing->GetOutput(), filename);
	//temp = writeImage< DImageType_2D, UCImageType_2D >( smoothing->GetOutput(), filename);

	std::cout<< "gradient ....."<<std::endl;
	double sigma = 3.0;
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<DImageType_2D, DImageType_2D>  GradientFilterType;
	GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
	gradientMagnitude->SetInput( smoothing->GetOutput());
	gradientMagnitude->SetSigma( this->parameters.gradient_sigma);
	gradientMagnitude->Update();
	filename = "gradientimage.nrrd";
	temp = writeImage< DImageType_2D >( gradientMagnitude->GetOutput(), filename);
	//temp = writeImage< DImageType_2D, UCImageType_2D >( gradientMagnitude->GetOutput(), filename);

	std::cout<< "sigmoid....."<<std::endl;
	double alfa = -2;
	double beta = 8;
	typedef itk::SigmoidImageFilter<DImageType_2D, DImageType_2D> SigmoidFilterType;
	SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	sigmoid->SetAlpha( this->parameters.sigmoid_alpha);
	sigmoid->SetBeta(  this->parameters.sigmoid_beta);
	sigmoid->SetOutputMinimum( 0.0);
	sigmoid->SetOutputMaximum( 1.0);
	sigmoid->SetInput( gradientMagnitude->GetOutput());
	sigmoid->Update();
	filename = "sigmoidimage.nrrd";
	temp = writeImage< DImageType_2D >( sigmoid->GetOutput(), filename);
	//temp = writeImage< DImageType_2D, UCImageType_2D >( sigmoid->GetOutput(), filename);

	typedef itk::FastMarchingImageFilter< DImageType_2D, DImageType_2D > FastMarchingFilterType;
	typedef FastMarchingFilterType::NodeContainer    NodeContainer;
    typedef FastMarchingFilterType::NodeType      NodeType;
	FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
	NodeContainer::Pointer seeds = NodeContainer::New();
	seeds->Initialize();
	for( int i = 0; i< this->seedNum; i++ )
	{
		double seedValue = -5;
		DImageType_2D::IndexType  seedPosition;
		seedPosition[0] = this->seedVector[i][0];
		seedPosition[1] = this->seedVector[i][1];
		NodeType node;
		node.SetValue( seedValue);
		node.SetIndex( seedPosition);
		seeds->InsertElement( i, node);
	}

	std::cout<< "fastmatching....."<<std::endl;
	int timethreshold = 15; 
	//fastMarching->SetInput( sigmoid->GetOutput() );
	fastMarching->SetTrialPoints(  seeds);
	fastMarching->SetOutputSize(this->imagesize );
	fastMarching->SetStoppingValue(  this->parameters.fastmarching_timethreshold);
	fastMarching->Update();
	fastMarching->SetSpeedConstant( 1.0 );
	filename = "fastmarching.nrrd";
	temp = writeImage< DImageType_2D >( fastMarching->GetOutput(), filename);
	//temp = writeImage< DImageType_2D, UCImageType_2D >( fastMarching->GetOutput(), filename);

	std::cout<< "Geodesic Active contour....."<<std::endl;
	double curvatureScaling = 1.0; 
	double advectScaling = 1.0;
	double PropagationScaling = 2.0;
	double rmsThres;
	int maxIterations = 10;
	typedef itk::GeodesicActiveContourLevelSetImageFilter< DImageType_2D, DImageType_2D, DPixelType > GeodesicActiveContourFilterType;
	GeodesicActiveContourFilterType::Pointer geodesicActiveContour = GeodesicActiveContourFilterType::New();
	geodesicActiveContour->SetInput(fastMarching->GetOutput());
	geodesicActiveContour->SetFeatureImage(sigmoid->GetOutput());
	geodesicActiveContour->SetPropagationScaling( this->parameters.ac_propagation);
	geodesicActiveContour->SetCurvatureScaling( this->parameters.ac_curvature);
	geodesicActiveContour->SetAdvectionScaling( this->parameters.ac_advect);
	geodesicActiveContour->SetMaximumRMSError( this->parameters.ac_maximumRMSE);
	geodesicActiveContour->SetNumberOfIterations( this->parameters.ac_iteration);
	geodesicActiveContour->Update();

	std::cout << "No. elpased iterations: " << geodesicActiveContour->GetElapsedIterations() << std::endl;
	std::cout << "RMS Error: " << geodesicActiveContour->GetRMSChange() << std::endl;

	std::cout<< "Thresholding..."<<endl;
	filename = "activecontourimage.nrrd";
	temp = writeImage< DImageType_2D >( geodesicActiveContour->GetOutput(), filename);
	//temp = writeImage< DImageType_2D, UCImageType_2D >( geodesicActiveContour->GetOutput(), filename);

	typedef itk::BinaryThresholdImageFilter< DImageType_2D, UCImageType_2D> BinaryThresholdingFilterType;
	BinaryThresholdingFilterType::Pointer thresholder = BinaryThresholdingFilterType::New();
	thresholder->SetLowerThreshold( -10000);
	thresholder->SetUpperThreshold(0.0);                                                           
	thresholder->SetOutsideValue( 0);
	thresholder->SetInsideValue(255);
	thresholder->SetInput( geodesicActiveContour->GetOutput());
	filename = "thresholderimage.tif";
	temp = writeImage< UCImageType_2D >( thresholder->GetOutput(), filename);
	//temp = writeImage< UIImageType_2D, UCImageType_2D >( thresholder->GetOutput(), filename);

	/// Label image, recaculate centroids

	//LabelFilterType::Pointer label = LabelFilterType::New();
	//label->SetInput(thresholder->GetOutput());

	//RelabelFilterType::Pointer relabel = RelabelFilterType::New();
	//relabel->SetInput( label->GetOutput());
	//relabel->SetMinimumObjectSize( minObjSize);  

	//relabel->Update();
	//SegmentedImageType::Pointer somas = relabel->GetOutput();

	//LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
	//labelGeometryImageFilter->SetInput( somas);
	//labelGeometryImageFilter->CalculatePixelIndicesOff();
	//labelGeometryImageFilter->CalculateOrientedBoundingBoxOff();
	//labelGeometryImageFilter->CalculateOrientedLabelRegionsOff();
	//labelGeometryImageFilter->CalculateOrientedIntensityRegionsOff();
	//labelGeometryImageFilter->Update();
	//LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
	//LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt =  allLabels.begin();

	//somaCentroids.clear();
	//for( allLabelsIt++; allLabelsIt != allLabels.end(); allLabelsIt++ )
	//{
	//	LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
	//	vtkSmartPointer<vtkVariantArray> row = vtkSmartPointer<vtkVariantArray>::New();
	//	LabelGeometryImageFilterType::LabelPointType point = labelGeometryImageFilter->GetCentroid(labelValue);
	//	itk::Index<3> index;
	//	index[0] = point[0];
	//	index[1] = point[1];
	//	index[2] = point[2];
	//	somaCentroids.push_back(index);
	//}
	//return somas;
}
