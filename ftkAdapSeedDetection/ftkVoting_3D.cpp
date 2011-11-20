// ############################################################################################################################################################################
#include "ftkVoting_3D.h"


int ftkWPoint3D::nx = 0;
int ftkWPoint3D::ny = 0;
int ftkWPoint3D::nz = 0;

int ftkCone3D::nx = 0;
int ftkCone3D::ny = 0;
int ftkCone3D::nz = 0;
int ftkCone3D::contador_1 = 0;

int ftkBins3D::nx = 0;
int ftkBins3D::ny = 0;
int ftkBins3D::nz = 0;

// Static constants
const double ftkVoting_3D::epsilon((double)0.001);
const double ftkVoting_3D::pi = (double)3.1415927;

// ############################################################################################################################################################################
ftkVoting_3D::ftkVoting_3D(){

	// Initial Parameters
	//reg_type (have not underestand what is this for)
	_hmin = 1;
	_hmax = 20;
	_radius = 4;
	_min_grad = 0.01;
	//_threshold (for picking seed points, not necesary for now)
	_scale = 1.5; // Scale for computing the gradient using DoG
	//_zc_only TRUE - voting points should also be the zero-crossing points of the image; FALSE - otherwise

	// Default Parameters for quantizing the direction of voting
	ntheta = 256;
	delta_theta = 2*pi/ntheta;
	_NN_dir = 1000;
	_sigmaG = 1.5;

}

//// ############################################################################################################################################################################
void ftkVoting_3D::setParams(	int hmin, int hmax,	int radius,	double min_grad, double scale ){
	_hmin = hmin;
	_hmax = hmax;
	_radius = radius;
	_min_grad = min_grad;
	_scale = scale; // Scale for computing the gradient using DoG
}

// ############################################################################################################################################################################
void ftkVoting_3D::compute(nftkVotingGlobal::InputImageType_3D::Pointer inputImage)
{
	nx = inputImage->GetLargestPossibleRegion().GetSize()[0];
	ny = inputImage->GetLargestPossibleRegion().GetSize()[1];
	nz = inputImage->GetLargestPossibleRegion().GetSize()[2];

	_voting_points = std::vector<VPoint3D>(); // These are the points that actually vote
	_voting_points.reserve(npix/2);		// To speed up, reserve some space (on average what will be used)



	std::cout << std::endl << "--> COMPUTING DERIVATIVES <--";

	// Derivatives using Gaussian
	typedef itk::RecursiveGaussianImageFilter<nftkVotingGlobal::InputImageType_3D,nftkVotingGlobal::InputImageType_3D > ShortFilterType;
	typedef itk::ImageDuplicator< nftkVotingGlobal::InputImageType_3D > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();

	//************** 1st order derivatives ***************************
	ShortFilterType::Pointer dx = ShortFilterType::New();
	dx->SetDirection( 0 );//dx works along x
	dx->SetSigma( _sigmaG );
	// Specify the derivative order in each direction
	dx->SetFirstOrder();
	// build the processing pipe
	dx->SetInput( inputImage );

	// Calculate the derivative in X
	std::cout << std::endl << "Derivative X. ";
	_t1_begin=clock();
	try{
		dx->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	duplicator->SetInputImage( dx->GetOutput() );
	try{
		duplicator->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	nftkVotingGlobal::InputImageType_3D::Pointer votingDirX_3D = duplicator->GetOutput(); //first order derivative along x
	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(_t1_end,_t1_begin)) << " s";


	// Calculate the derivative in Y
	std::cout << std::endl << "Derivative Y. ";
	_t1_begin=clock();
	dx->SetDirection( 1 );//dx works along y
	try{
		dx->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	try{
		duplicator->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	nftkVotingGlobal::InputImageType_3D::Pointer votingDirY_3D = duplicator->GetOutput();
	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(_t1_end,_t1_begin)) << " s";


	// Calculate the derivative in Z
	std::cout << std::endl << "Derivative Z. ";
	_t1_begin=clock();
	dx->SetDirection( 2 );//dx works along z
	try{
		dx->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	try{
		duplicator->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	nftkVotingGlobal::InputImageType_3D::Pointer votingDirZ_3D = duplicator->GetOutput();
	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(_t1_end,_t1_begin)) << " s";


	nftkVotingGlobal::InputImageType_3D::PixelType * votingDirX_3DArray = votingDirX_3D->GetBufferPointer();
	nftkVotingGlobal::InputImageType_3D::PixelType * votingDirY_3DArray = votingDirY_3D->GetBufferPointer();
	nftkVotingGlobal::InputImageType_3D::PixelType * votingDirZ_3DArray = votingDirZ_3D->GetBufferPointer();

	// Gradient magnitud
	// Duplicate image X (to create a copy)
	typedef itk::ImageDuplicator< nftkVotingGlobal::InputImageType_3D > DuplicatorType_2;
	DuplicatorType_2::Pointer duplicator_2 = DuplicatorType_2::New();
	duplicator_2->SetInputImage( votingDirX_3D );
	try{
		duplicator_2->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	nftkVotingGlobal::InputImageType_3D::Pointer votingDirXYZ_3D = duplicator_2->GetOutput(); //first order derivative along z

	nftkVotingGlobal::InputImageType_3D::PixelType * votingDirXYZ_3DArray = votingDirXYZ_3D->GetBufferPointer();

	// Gradient magnitude and normalize gradient X, Y and Z
#pragma omp parallel for
	for( int tt=0; tt<npix; tt++ ){
		votingDirXYZ_3DArray[tt] = sqrt(votingDirX_3DArray[tt]*votingDirX_3DArray[tt]+votingDirY_3DArray[tt]*votingDirY_3DArray[tt]+votingDirZ_3DArray[tt]*votingDirZ_3DArray[tt]);
		if(votingDirXYZ_3DArray[tt]>epsilon){
			votingDirX_3DArray[tt] = votingDirX_3DArray[tt]/votingDirXYZ_3DArray[tt];
			votingDirY_3DArray[tt] = votingDirY_3DArray[tt]/votingDirXYZ_3DArray[tt];
			votingDirZ_3DArray[tt] = votingDirZ_3DArray[tt]/votingDirXYZ_3DArray[tt];
		}
	}

	// Calculate the Maximum value of the gradient
	typedef itk::MinimumMaximumImageCalculator < nftkVotingGlobal::InputImageType_3D > ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
	imageCalculatorFilter->SetImage( votingDirXYZ_3D);
	imageCalculatorFilter->Compute();

	nftkVotingGlobal::InputImageType_3D::PixelType maxVotVal = imageCalculatorFilter->GetMaximum();

	// Scale the values of the gradient by the maximum value of the gradiente, to ensure that the maximum value is 1

	std::cout << std::endl << "Scalling the values of the gradiet (0 to 1) -> ";
	_t1_begin=clock();
#pragma omp parallel for
	for(int i=0; i<npix; i++) {
		if (votingDirXYZ_3DArray[i]<0.01)
		{
			votingDirXYZ_3DArray[i] = 0;
		}
		else
		{
			votingDirXYZ_3DArray[i] /= maxVotVal;
		}
	}
	std::cout << "Done.";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(_t1_end,_t1_begin)) << " s";

	//string filenameDerivativeXYZ_Thre = "output\\out_ftkDerivativeXYZ_Thre.tif";
	//if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(votingDirXYZ_3D, filenameDerivativeXYZ_Thre.c_str() )){
	//	std::std::cout<<std::endl<<"\tProblema escribiendo";
	//}




	//std::cout<<std::endl<<"Write the derivative images";
	//
	//	// Save the derivatives
	//	string filenameDerivativeX = "output\\out_ftkDerivativeX.tif";
	//	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(votingDirX_3D, filenameDerivativeX.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//	string filenameDerivativeY = "output\\out_ftkDerivativeY.tif";
	//	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(votingDirY_3D, filenameDerivativeY.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//	string filenameDerivativeZ = "output\\out_ftkDerivativeZ.tif";
	//	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(votingDirZ_3D, filenameDerivativeZ.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//	string filenameDerivativeXYZ = "output\\out_ftkDerivativeXYZ.tif";
	//	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(votingDirXYZ_3D, filenameDerivativeXYZ.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}




	//std::cout<<std::endl<<"Threshold the gradients and store the results"; // The thresholded gradient magnitude will become the voting pixels
	//
	//	for( unsigned int tt=0; tt<npix; tt++ ){
	//		votingDirX_3DArray[tt] = abs(votingDirX_3DArray[tt]);
	//		if( votingDirX_3DArray[tt] < 0.01 ){
	//			votingDirX_3DArray[tt] = 0;
	//		}
	//		//std::cout<<std::endl<<votingDirZ_3DArray[tt];
	//	}
	//	string filenameDerivativeX_Thre = "output\\out_ftkDerivativeX_Thre.tif";
	//	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(votingDirX_3D, filenameDerivativeX_Thre.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//
	//	for( unsigned int tt=0; tt<npix; tt++ ){
	//		votingDirY_3DArray[tt] = abs(votingDirY_3DArray[tt]);
	//		if( votingDirY_3DArray[tt] < 0.01 ){
	//			votingDirY_3DArray[tt] = 0;
	//		}
	//		//std::cout<<std::endl<<votingDirZ_3DArray[tt];
	//	}
	//	string filenameDerivativeY_Thre = "output\\out_ftkDerivativeY_Thre.tif";
	//	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(votingDirY_3D, filenameDerivativeY_Thre.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//
	//	for( unsigned int tt=0; tt<npix; tt++ ){
	//		votingDirZ_3DArray[tt] = abs(votingDirZ_3DArray[tt]);
	//		if( votingDirZ_3DArray[tt] < 0.01 ){
	//			votingDirZ_3DArray[tt] = 0;
	//		}
	//		//std::cout<<std::endl<<votingDirZ_3DArray[tt];
	//	}
	//	string filenameDerivativeZ_Thre = "output\\out_ftkDerivativeZ_Thre.tif";
	//	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(votingDirZ_3D, filenameDerivativeZ_Thre.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//
	//	for( unsigned int tt=0; tt<npix; tt++ ){
	//		if( votingDirXYZ_3DArray[tt] < 0.01 ){
	//			votingDirXYZ_3DArray[tt] = 0;
	//		}
	//		//std::cout<<std::endl<<votingDirZ_3DArray[tt];
	//	}
	//	string filenameDerivativeXYZ_Thre = "output\\out_ftkDerivativeXYZ_Thre.tif";
	//	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(votingDirXYZ_3D, filenameDerivativeXYZ_Thre.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}


	std::cout << std::endl << "Now we are going to compute the cones";
	computeCones(_hmin, _hmax, _radius);
	std::cout << " -> Done.";
	//computeCones_prob(_hmin, _hmax, _radius); // Just for one moment


	// Padding
	bw = sqrt((double)(_radius*_radius + _hmax*_hmax)) + 3; // Que es esto ??
	const int bw2 = 2*bw;

	// Calculate the voting poitns
	std::cout << std::endl << "Computing the voting points -> ";
	_t1_begin=clock();
	int max1=0;
	int max2=0;
// This can be parallel
	std::pair< int, int > indx;
	VPoint3D vp;
	for(int z = 0; z < nz; z++)
	{
		for(int y = 0; y < ny; y++)
		{
			for(int x = 0; x < nx; x++)
			{
				int k = x+nx*y+nx*ny*z;
				if (votingDirXYZ_3DArray[k] > epsilon) // This is the original one, only vote if the mag of gradient is greater than epsilon (and gra_min) previously thresholded
				{
					indx = computeAngleIndex_3D(votingDirX_3DArray[k], votingDirY_3DArray[k], votingDirZ_3DArray[k]); // NO ESTAN INCLUYENDO LOS CEROS PORQUE ??
					vp.x = x + bw;
					vp.y = y + bw;
					vp.z = z + bw;
					vp.mag = votingDirXYZ_3DArray[k]; 
					vp.angIndex = indx;
					_voting_points.push_back(vp);

					if(indx.first>max1)
						max1=indx.first;
					if(indx.second>max2)
						max2=indx.second;
				}
			}
		}
	}
	std::cout << "Done.";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(_t1_end,_t1_begin)) << " s";

	std::cout << std::endl << "The number of voting points: "<<_voting_points.size();
	std::cout << std::endl << "MAX ANGLE: "<<max1<<" "<<max2;


	// Calculate the voting poitns
	std::cout << std::endl << "Calculate Voting Points 3D -> ";
	_t1_begin=clock();
// This can be parallel
	int direction_vote;
	VPoint3D vp_3D;
	for(int z = 0; z < nz; z++)
	{
		for(int y = 0; y < ny; y++)
		{
			for(int x = 0; x < nx; x++)
			{
				int k = x+nx*y+nx*ny*z;
				if (votingDirXYZ_3DArray[k] > epsilon) // This is the original one, only vote if the mag of gradient is greater than epsilon (and gra_min) previously thresholded
				{
					vp_3D.x = x + bw;
					vp_3D.y = y + bw;
					vp_3D.z = z + bw;
					vp_3D.mag = votingDirXYZ_3DArray[k]; 
					vp_3D.posOri = k; // This is the original position in the original image
					_voting_points_3D.push_back(vp_3D);
				}
			}
		}
	}
	std::cout << "Done.";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(_t1_end,_t1_begin)) << " s";

	// Calculate the direction of the voting poitns
	std::cout << std::endl << "Quantize Voting Points Directions in 3D -> ";
	_t1_begin=clock();
// This can be parallel
	int limit = _voting_points_3D.size();
#pragma omp parallel for
	for( int i = 0; i<limit; ++i )\
	{
		int k = _voting_points_3D[i].posOri;
		_voting_points_3D[i].direc_vote = computeDirecIndex_3D(votingDirX_3DArray[k], votingDirY_3DArray[k], votingDirZ_3DArray[k]); 
	}
	std::cout << "Done.";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(_t1_end,_t1_begin)) << " s";





	nx += bw2;
	ny += bw2;
	nz += bw2;
	npix = nx*ny*nz;

	ftkWPoint3D::setImageSize(nx,ny,nz);
	ftkBins3D::setImageSize(nx,ny,nz);
	ftkCone3D::setImageSize(nx,ny,nz);

	vote();

	// Switch back to the original size
	nx -= bw2;
	ny -= bw2;
	nz -= bw2;
	npix = nx*ny*nz;

	//_votingVotes = VotingDirType_3D::New();
	//_votingVotes->SetRegions( inputImage->GetRequestedRegion() );
	//_votingVotes->Allocate();

	// Copy sum votes to 














	//	// Derivative XXXXXXXXXX
	//	VotingDirType_3D::Pointer votingDirX = VotingDirType_3D::New();
	//	votingDirX->SetRegions( inputImage->GetRequestedRegion() );
	//	votingDirX->Allocate();
	//
	//	itk::SobelOperator< VotingDirPixelType, 2 > sobelOperator;
	//	sobelOperator.SetDirection( 0 );
	//	sobelOperator.CreateDirectional();
	//	typedef itk::ImageRegionIterator< VotingDirType_3D > IteratorType3;
	//	IteratorType3 out2( votingDirX, votingDirX->GetRequestedRegion() );
	//
	//	typedef itk::ConstNeighborhoodIterator< nftkVotingGlobal::InputImageType > NeighborhoodIteratorType2;
	//	typedef itk::ImageRegionIterator< nftkVotingGlobal::InputImageType > IteratorType2;
	//	NeighborhoodIteratorType2::RadiusType radius2 = sobelOperator.GetRadius();
	//	NeighborhoodIteratorType2 it2( radius2, inputImage, inputImage->GetRequestedRegion() );
	//	itk::NeighborhoodInnerProduct< nftkVotingGlobal::InputImageType > innerProduct;
	//
	//	for (it2.GoToBegin(), out2.GoToBegin(); !it2.IsAtEnd(); ++it2, ++out2)
	//	{
	//		out2.Set( innerProduct( it2, sobelOperator ) );
	//	}
	//
	//	// Derivative YYYYYYYYYY
	//	VotingDirType_3D::Pointer votingDirY = VotingDirType_3D::New();
	//	votingDirY->SetRegions( inputImage->GetRequestedRegion() );
	//	votingDirY->Allocate();
	//
	//	itk::SobelOperator<VotingDirPixelType, 2> sobelOperatorY;
	//	sobelOperatorY.SetDirection( 1 );
	//	sobelOperatorY.CreateDirectional();
	//	typedef itk::ImageRegionIterator< VotingDirType_3D > IteratorType4;
	//	IteratorType4 out3( votingDirY, votingDirY->GetRequestedRegion() );
	//
	//	typedef itk::ConstNeighborhoodIterator< nftkVotingGlobal::InputImageType > NeighborhoodIteratorType3;
	//	typedef itk::ImageRegionIterator< nftkVotingGlobal::InputImageType > IteratorType3;
	//	NeighborhoodIteratorType3::RadiusType radius3 = sobelOperatorY.GetRadius();
	//	NeighborhoodIteratorType3 it3( radius3, inputImage, inputImage->GetRequestedRegion() );
	//	itk::NeighborhoodInnerProduct< nftkVotingGlobal::InputImageType > innerProduct2;
	//
	//	for (it3.GoToBegin(), out3.GoToBegin(); !it3.IsAtEnd(); ++it3, ++out3)
	//	{
	//		out3.Set( innerProduct2( it3, sobelOperatorY ) );
	//	}
	//
	//	// Magnitude Image
	//	VotingDirType_3D::Pointer votingMagImage = VotingDirType_3D::New(); //VotingDirType_3D = double
	//	votingMagImage->SetRegions( inputImage->GetRequestedRegion()); // IMPORTANTE PARA CREAR UNA IMAGEN NUEVA EN BASE A UNA QUE YA EXISTE EN VEZ DE PONERME A LEER LOS TAMANOS Y LAS REGIONES
	//	votingMagImage->Allocate();
	//
	//	typedef itk::ImageRegionIteratorWithIndex< VotingDirType_3D > ITVotingMag;
	//	ITVotingMag iVotingMag(votingMagImage, votingMagImage->GetLargestPossibleRegion() );
	//
	//	typedef itk::ImageRegionIteratorWithIndex< VotingDirType_3D > ITVotingDir;
	//	ITVotingDir iVotingDirX(votingDirX, votingDirX->GetLargestPossibleRegion() );
	//	ITVotingDir iVotingDirY(votingDirY, votingDirY->GetLargestPossibleRegion() );
	//
	//	for ( iVotingMag.GoToBegin(),iVotingDirX.GoToBegin(),iVotingDirY.GoToBegin(); !iVotingMag.IsAtEnd(); ++iVotingMag, ++iVotingDirX, ++iVotingDirY ){
	//		iVotingMag.Set(sqrt( iVotingDirX.Get()*iVotingDirX.Get() + iVotingDirY.Get()*iVotingDirY.Get() ));
	//		if(iVotingMag.Get()>epsilon){
	//			iVotingDirX.Set(iVotingDirX.Get()/iVotingMag.Get());
	//			iVotingDirY.Set(iVotingDirY.Get()/iVotingMag.Get());
	//		}
	//	}
	//
	//	// Magnitude Image Binary (just to store in file the voting pixels)
	//	VotingDirType_3D::Pointer votingMagImage_bin = VotingDirType_3D::New();
	//	votingMagImage_bin->SetRegions( inputImage->GetRequestedRegion() );
	//	votingMagImage_bin->Allocate();
	//
	//	// Canny edge detection
	//	typedef itk::CannyEdgeDetectionImageFilter <VotingDirType_3D, VotingDirType_3D> CannyEdgeDetectionImageFilterType;
	//	CannyEdgeDetectionImageFilterType::Pointer cannyFilter = CannyEdgeDetectionImageFilterType::New();
	//	cannyFilter->SetInput(inputImage);
	//	cannyFilter->SetVariance( 2.5 );
	//	cannyFilter->SetUpperThreshold( 0.03/*0.0238*/ );
	//	cannyFilter->SetLowerThreshold( 0.005/*0.0175*/ );
	//	cannyFilter->Update();
	//
	//	// Save the canny edge result
	//	string filenameCanny = "output\\out_ImageOfCanny.jpg";
	//	if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(cannyFilter->GetOutput(), filenameCanny.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//
	//
	//
	//	// Arrays of data
	//	nftkVotingGlobal::InputImageType::PixelType * votingImaArray = inputImage->GetBufferPointer();
	//	VotingDirType_3D::PixelType * votingDirXArray = votingDirX->GetBufferPointer();
	//	VotingDirType_3D::PixelType * votingDirYArray = votingDirY->GetBufferPointer();
	//	VotingDirType_3D::PixelType * votingMagArray = votingMagImage->GetBufferPointer();
	//	VotingDirType_3D::PixelType * votingMagArray_bin = votingMagImage_bin->GetBufferPointer();
	//	VotingDirType_3D::PixelType * votingCannyArray = cannyFilter->GetOutput()->GetBufferPointer();
	//	
	//
	//	//double maxgrad = votingMagImag
	//
	//	typedef itk::MinimumMaximumImageCalculator < VotingDirType_3D > ImageCalculatorFilterType;
	//	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
	//	imageCalculatorFilter->SetImage( votingMagImage );
	//	imageCalculatorFilter->Compute();
	//
	//	VotingDirType_3D::PixelType maxVotVal = imageCalculatorFilter->GetMaximum();
	//
	//	// Scale the values of the gradient by the maximum value of the gradiente, to ensure that the maximum value is 1
	//	for(int i=0; i<npix; i++) {
	//			if (votingMagArray[i]<_min_grad)
	//			{
	//				votingMagArray[i] = 0;
	//			}
	//			else
	//			{
	//				votingMagArray[i] /= maxVotVal;
	//			}
	//	}
	//
	//
	//	// Testing to store the resulting voting image base on the gradient
	//	// Save the canny edge result
	//	string filenameGradVot = "output\\out_ImageOfGradVot.jpg";
	//	if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(votingMagImage, filenameGradVot.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//	// Put 1 if mag array is different to zero
	//	for(int i=0; i<npix; i++) {
	//			votingMagArray_bin[i] = 0;
	//			if (votingMagArray[i]!=0)
	//			{
	//				votingMagArray_bin[i] = 1;
	//			}
	//	}
	//	string filenameGradVot_bin = "output\\out_ImageOfGradVot_bin.jpg";
	//	if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(votingMagImage_bin, filenameGradVot_bin.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//	// Put 1 if mag array intersect with canny edge is di
	//	for(int i=0; i<npix; i++) {
	//			votingMagArray_bin[i] = 0;
	//			if (votingMagArray[i]!=0 && votingCannyArray[i]!=0)
	//			{
	//				votingMagArray_bin[i] = 1;
	//			}
	//	}
	//	string filenameGradVotInterCanny = "output\\out_ImageOfGradVotInterCanny.jpg";
	//	if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(votingMagImage_bin, filenameGradVotInterCanny.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//
	//
	//
	//
	//
	//
	//	// Padding
	//	bw = sqrt((double)(_radius*_radius + _hmax*_hmax)) + 3; // Que es esto ??
	//	const int bw2 = 2*bw;
	//
	//
	//	// Calculate the voting poitns
	////int stop1;
	//	int indx;
	//	VPoint2D vp;
	//	int x, y, k;
	//	for(y = 0; y < ny; y++)
	//	{
	//		for(x = 0; x < nx; x++)
	//		{
	//			k = x+nx*y;
	//			if (votingMagArray[k] > epsilon) // This is the original one, only vote if the mag of gradient is greater than epsilon (and gra_min) previously thresholded
	//			//if(votingCannyArray[k]!=0) // This is new, use the canny edge detection to vote form this values (IT DOES NOT WORK TO USE THE CANNY EDGE DETECTION, INCREIBLE)
	//			{
	//				if ((indx=computeAngleIndex(votingDirXArray[k], votingDirYArray[k])) >= 0 ) // NO ESTAN INCLUYENDO LOS CEROS PORQUE ??
	//				{
	//					vp.x = x + bw;
	//					vp.y = y + bw;
	//					vp.mag = votingMagArray[k]; 
	//					vp.angIndex = indx;
	//					_voting_points.push_back(vp);
	//				}	
	//			}
	//		}
	//	}
	////std::cin>>stop1;
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//	nx += bw2;
	//	ny += bw2;
	//	npix = nx*ny;
	//
	//	ftkWPoint2D::setImageSize(nx,ny);
	//	ftkBins2D::setImageSize(nx,ny);
	//	ftkCone2D::setImageSize(nx,ny);
	//
	//	vote();
	//
	//	// Switch back to the original size
	//	nx -= bw2;
	//	ny -= bw2;
	//	npix = nx*ny;
	//
	//	//_votingVotes = VotingDirType_3D::New();
	//	//_votingVotes->SetRegions( inputImage->GetRequestedRegion() );
	//	//_votingVotes->Allocate();
	//
	//	// Copy sum votes to 


}


//void ftkVoting_3D::compute(const Image& I)
//{
//	//std::cout<<"Executing local version of circular voting "<<" Thresh :"<<_threshold<<std::endl;
//	nx = I.nx();
//	ny = I.ny();
//	npix = nx*ny;
//
//	_voting_points = std::std::vector<VPoint2D>();
//	_voting_points.reserve(npix/2);
//
//	const double reg_sign = _reg_type==DarkReg ? -1 : 1;
//
//	doubleImage Dx = dX(I, _scale);
//	doubleImage Dy = dY(I, _scale);
//	doubleImage G(nx, ny);
//	for(int i=0; i<npix; i++) {
//		Dx(i) *= reg_sign;
//		Dy(i) *= reg_sign;
//		double g = sqrt(Dx(i)*Dx(i)+Dy(i)*Dy(i));
//		G(i) = g;
//		if (g>epsilon) {
//			Dx(i) /= g;
//			Dy(i) /= g;
//		}
//	}
//
//	double maxgrad = G.max();
//	if (maxgrad>epsilon) 
//		for(int i=0; i<npix; i++) {
//			if (G(i)<_min_grad)
//				G(i) = 0;
//			else if( !_abs_thresh )
//				G(i) /= maxgrad;
//		}
//
//		if (_zc_only) {
//			ByteImage ZC = edge(features(I, _scale, BrightZC));    
//			for(int i=0; i<npix; i++) 
//				if (!ZC(i))
//					G(i) = 0;
//		}
//
//		// Weighting   
//		if (_weighting_method==PowMethod && fabs(_weighting_param-1)>epsilon) 
//			for(int i=0; i<npix; i++) 
//				if (G(i)>epsilon)
//					G(i) = pow(G(i), _weighting_param);
//				else if (_weighting_method==ExpMethod) 
//					for(int i=0; i<npix; i++) 
//						if (G(i)>epsilon)
//							G(i) = exp( G(i) *  _weighting_param );
//
//////		const int nx_1 = nx-1; // Para que esto ?????
//////		const int ny_1 = ny-1;
//
//		bw = sqrt((double)(_radius*_radius + _hmax*_hmax)) + 3; // Que es esto ??
//		const int bw2 = 2*bw;
//
//		int indx;
//		VPoint2D vp;
//		int x, y, k;
//		for(y = 0; y < ny; y++)
//			for(x = 0; x < nx; x++) {      
//				k = I.pos(x, y);
//				if (G(k) > epsilon) {
//					if ((indx=computeAngleIndex(Dx(k), Dy(k))) >= 0 ) {
//						vp.x = x + bw;
//						vp.y = y + bw;
//						vp.mag = G(k); 
//						vp.angIndex = indx;
//						_voting_points.push_back(vp);
//					}	
//				}
//			}
//
//			//std::cout << "#Voting points = " << _voting_points.size() << std::endl;
//
//			nx += bw2;
//			ny += bw2;
//			npix = nx*ny;
//
//			WPoint2D::setImageSize(nx, ny);
//			ConePlane2D::setImageSize(nx, ny);
//			Cone2D::setImageSize(nx, ny);
//
//			vote();
//
//			// Switch back to the original size
//			nx -= bw2;
//			ny -= bw2;
//			npix = nx*ny;
//			_vote = doubleImage(nx, ny);
//
//			double val, maxval=0;
//
//			for(y = 0; y < ny; y++)
//				for(x = 0; x < nx; x++) {
//					val = _sum(x+bw, y+bw);
//					_vote(x, y) = val;
//					maxval = max<double>(maxval, val);
//				}
//				//  return;
//				int nc;
//				//std::cout<<"Voting Threshold "<<_threshold<<" max Value "<<maxval<<" thresh * maxval:" <<_threshold * maxval<<std::endl;
//				// ByteImage bIm = threshold(_vote,_threshold*maxval) ;
//				// ByteImage bIm = threshold(_vote,_threshold) ;
//				ByteImage bIm(_vote.nx(),_vote.ny(),_vote.nz());
//				for (int count=0;count<bIm.npix();count++){
//					if( _abs_thresh )
//						bIm.set(count,_vote.pixel(count) > _threshold?255:0 );
//					else
//						bIm.set(count,_vote.pixel(count) > _threshold*maxval?255:0 );
//				}
//
//				ShortImage label = regions(bIm, nc) ;
//				std::vector<RegionFeature> rf = regionFeatures(label, nc, _vote);
//				_vote.reset(0);
//				_center_weights.clear();
//				for(std::vector<RegionFeature>::iterator it=rf.begin(); it!=rf.end(); it++) {
//					if (_vote.inImage(it->xc, it->yc)) {
//						_vote(it->xc, it->yc) = it->mean;
//						_centers.push_back(Point( it->xc, it->yc ));
//						_center_weights.push_back(it->mean);
//					}   
//				}
//
//				for(std::vector<VPoint2D>::iterator it=_voting_points.begin(); it!=_voting_points.end(); it++) {
//					it->x -= bw;
//					it->y -= bw;
//					it->xc -= bw;
//					it->yc -= bw;
//				}
//}
//
// ############################################################################################################################################################################
ftkVoting_3D::~ftkVoting_3D()
{
}

// ############################################################################################################################################################################
void ftkVoting_3D::setPrefix(const std::string& p)
{
	//_prefix = p;
}

//// ############################################################################################################################################################################
//int ftkVoting_3D::computeAngleIndex(double dx, double dy) const 
//{
//	// NEEDS REVIEW
//	double a = acos(dx);
//	if (dy<0) {
//		a = 2*pi-a;
//	}
//
//	//int indx_theta = nftkVotingGlobal::round_double(a/delta_theta); // De por si esta divisino no veo como pueda ser mas de 255
//	//return max(0, min(ntheta-1, indx_theta)); // Se asegura que el maximo devuelto sea 255 (256 nunca va a salir de aca)
//
//	int indx_theta = floor(a/delta_theta); // De por si esta divisino no veo como pueda ser mas de 255
//	if(indx_theta<0)
//	{int rr;
//	std::cout<<std::endl<<"Error en compute angle: "<<indx_theta<<" "<<dx<<" "<<dy;
//	std::cin>>rr;
//	}
//	return max(0, min(ntheta-1, indx_theta));
//}

// ############################################################################################################################################################################
std::pair<int,int> ftkVoting_3D::computeAngleIndex_3D(double dx, double dy, double dz) const 
{

	//dx = i/dmag;
	//dy = j/dmag;
	//dz = kk/dmag;

	int phi_qu;
	int tetaa_qu;
	double phi;
	double tetaa;

	//           if( (dx==0) && (dy==0) )
	//		{
	//               betaa_qu = 0;
	//               tetaa_qu = 0;
	//               if(dz<0)
	//			{
	//                   betaa_qu = 0;
	//                   tetaa_qu = _conesPru_3D;
	//			}
	//               //ban = 1;
	//		}
	//           else if( (dx==0) && (dz==0) )
	//		{
	//               betaa_qu = 127;
	//               tetaa_qu = 127;
	//               if(dy<0)
	//			{
	//                   betaa_qu = 192;
	//                   tetaa_qu = 127;
	//			}
	//               //ban = 1;
	//		}
	//           else if( (dy==0) && (dz==0) )
	//		{
	//               betaa_qu = 0;
	//               tetaa_qu = 127;
	//               if(dx<0)
	//			{
	//                   betaa_qu = _conesPru_3D;
	//                   tetaa_qu = 127; 
	//			}
	//               //ban = 1;
	//		}
	//           else
	//		{
	//               betaa = atan(dy/dx);
	//               if(dx==0)
	//			{
	//                   betaa = pi/2;
	//                   if(dy<0)
	//				{
	//                       betaa = 3*pi/4;
	//				}
	//			}
	//               else if(dy==0)
	//			{
	//                   betaa = 0;
	//                   if(dx<0)
	//				{
	//                       betaa = pi;
	//				}
	//			}
	//               else if(dy>0 && dx<0 )
	//			{
	//                   betaa = pi + betaa;
	//			}
	//               else if(dy<0 && dx>0 )
	//			{
	//                   betaa = 2*pi + betaa;
	//			}
	//               else if(dy<0 && dx<0 )
	//			{
	//                   betaa = pi + betaa;
	//			}
	//               
	//               tetaa = acos(dz);
	//               
	//               betaa_qu = floor(betaa/delta_theta);
	//               tetaa_qu = floor(tetaa/delta_theta);
	//               //ban = 1;
	//		}

	tetaa = acos(dz);

	phi = atan(dy/dx);
	if( (dx==0) && (dy==0) )
	{
		phi = 0; // By default
	}
	else if(dx==0)
	{
		phi = pi/2;
		if(dy<0)
		{
			phi = 3*pi/2;
		}
	}
	else if(dy==0)
	{
		phi = 0;
		if(dx<0)
		{
			phi = pi;
		}
	}
	else if((dy>0) && (dx<0) )
	{
		phi = pi + phi;
	}
	else if((dy<0) && (dx>0) )
	{
		phi = 2*pi + phi;
	}
	else if((dy<0) && (dx<0) )
	{
		phi = pi + phi;
	}




	tetaa_qu = floor(tetaa/delta_theta);
	phi_qu = floor(phi/delta_theta);

	std::pair<int,int> indx_theta;
	indx_theta.first = tetaa_qu; // De por si esta divisino no veo como pueda ser mas de 255
	indx_theta.second = phi_qu; // De por si esta divisino no veo como pueda ser mas de 255

	//if( indx_theta.second==24)
	//{
	//	std::cout<<"48: "<<dx<<" "<<dy<<" "<<dz<<" "<<phi<<" "<<phi_qu<<" "<<tetaa<<" "<<tetaa_qu;
	//}


	if(indx_theta.first<0 || indx_theta.first>127 || indx_theta.second<0 || indx_theta.second>255 )
	{int rr;
	std::cout<<std::endl<<"Error en compute angle: "<<indx_theta.first<<" "<<indx_theta.second<<" "<<dx<<" "<<dy<<" "<<dz;
	std::cin>>rr;
	}


	//// NEEDS REVIEW
	//double a = acos(dx);
	//if (dy<0) {
	//	a = 2*pi-a;
	//}

	//// NEEDS REVIEW
	//double b = acos(dz);
	//if (dz<0) {
	//	b = pi-b;
	//}

	////int indx_theta = nftkVotingGlobal::round_double(a/delta_theta); // De por si esta divisino no veo como pueda ser mas de 255
	////return max(0, min(ntheta-1, indx_theta)); // Se asegura que el maximo devuelto sea 255 (256 nunca va a salir de aca)

	//pair<int,int> indx_theta;
	//indx_theta.first = floor(a/delta_theta); // De por si esta divisino no veo como pueda ser mas de 255
	//indx_theta.second = floor(b/(pi/ntheta)); // De por si esta divisino no veo como pueda ser mas de 255


	//if(indx_theta.first<0 || indx_theta.first>255 || indx_theta.second<0 || indx_theta.second>_conesPru_3D )
	//{int rr;
	//std::cout<<std::endl<<"Error en compute angle: "<<indx_theta.first<<" "<<indx_theta.second<<" "<<dx<<" "<<dy<<" "<<dz;
	//std::cin>>rr;
	//}


	//return max(0, min(ntheta-1, indx_theta));
	return indx_theta;
}

// ############################################################################################################################################################################
void inline ftkCone3D::vote(VotingDirType_3D::PixelType * p, const VPoint3D& vp)
{
	//	//int stop2;
	//	iterator from = begin(); //Pointer to Bin
	//	iterator to = end();     
	//	for(iterator it = from; it != to; it++) 
	//	{
	//		for(ftkBins2D::iterator it1=it->begin(); it1!=it->end(); ++it1) {
	//			p[it1->off] += vp.mag * it1->w;
	//			//p[it1->off] += 1; //Numer of votes, not using magnitude
	//			contador_1++;
	//		}      
	//	}
}

//// ############################################################################################################################################################################
//void inline ftkCone3D::vote_dir(VotingDirType_3D::PixelType * p, const VPoint3D& vp) //Vota y guarda direcciones de los votos
//{
//	//int stop2;
//	iterator from = begin(); //Pointer to Bin
//	iterator to = end();     
//	for(iterator it = from; it != to; it++) 
//	{
//		for(ftkBins2D::iterator it1=it->begin(); it1!=it->end(); ++it1) {
//			p[it1->off] += vp.mag * it1->w;
//			//p[it1->off] += 1; //Numer of votes, not using magnitude
//			contador_1++;
//		}      
//	}
//}

// ############################################################################################################################################################################
void inline ftkCone3D::vote(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int& dist)
{
	//iterator binRequired = begin();
	//ftkBins2D::iterator from = binRequired[dist].begin();
	//ftkBins2D::iterator to = binRequired[dist].end();

	//for(ftkBins2D::iterator it = from; it != to; it++) 
	//{
	//	p[it->off] += vp.mag * it->w; 
	//	//p[it->off] += 1; // Number of votes, not using magnitude
	//	contador_1++;
	//}
}

// ############################################################################################################################################################################
void inline ftkCone3D::vote(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int& dist, int &mag) // Incluye la probabilidad
{

	//iterator binRequired = begin();
	//ftkBins2D::iterator from = binRequired[dist].begin();
	//ftkBins2D::iterator to = binRequired[dist].end();

	//for(ftkBins2D::iterator it = from; it != to; it++) 
	//{
	//	//p[it->off] += vp.mag;// * (1/(double)mag)/*it->w*/; 
	//	p[it->off] += vp.mag * (1/(double)mag); 
	//	//p[it->off] += (1/(double)mag); 
	//	//p[it->off] += (1/(double)mag)/*it->w*/; 
	//	//std::cout<<" "<<mag;
	//	//p[it->off] += 1; // Number of votes, not using magnitude
	//	contador_1++;
	//}
}

// ############################################################################################################################################################################
void inline ftkCone3D::vote_dir(std::vector<std::vector<int> >& p_dir, VotingDirType_3D::PixelType * p, const VPoint3D& vp, int& dist, int& offset_1)
{

	//iterator binRequired = begin();
	//ftkBins2D::iterator from = binRequired[dist].begin();
	//ftkBins2D::iterator to = binRequired[dist].end();

	//for(ftkBins2D::iterator it = from; it != to; it++) 
	//{
	//	p[it->off] += vp.mag * it->w;
	//	p_dir.at(it->off+offset_1).at(vp.angIndex) += 1; // angIndex 0 255
	//	//p_dir[it->off][vp.angIndex] += 1;
	//	//p[it->off] += 1; // Number of votes, not using magnitude
	//	contador_1++;
	//}
}

//// ############################################################################################################################################################################
void inline ftkVoting_3D::updateDirection(VPoint3D& vp)
{
	//ftkBins2D::iterator max_point;
	//double maxvote = 0;
	//for( int tt=_voteDirec[vp.angIndex].first; tt!=_voteDirec[vp.angIndex].second; tt=(tt+1)%ntheta )
	//{
	//	ftkCone2D::iterator from = _conesPru[tt].begin();
	//	ftkCone2D::iterator to = _conesPru[tt].end();
	//	//_conesPru[tt].vote(p,vp);

	//	VotingDirType_3D::PixelType * votingSumArray = _votingSumVotes->GetBufferPointer()+vp.pos;
	//	for(ftkCone2D::iterator it = from; it != to; it++) 
	//	{
	//		for(ftkBins2D::iterator it1=it->begin(); it1!=it->end(); it1++)
	//		{
	//			if ( votingSumArray[it1->off]>maxvote) 
	//			{
	//				maxvote = max(maxvote, votingSumArray[it1->off]);
	//				max_point = it1;
	//			}    
	//		}
	//	} 
	//}

	//if (maxvote<epsilon)
	//	return;

	//double dx = max_point->x;
	//double dy = max_point->y;

	//double r = sqrt(dx*dx+dy*dy);
	//if (r>epsilon) {
	//	vp.xc = vp.x + max_point->x;
	//	vp.yc = vp.y + max_point->y;
	//	vp.angIndex = computeAngleIndex(dx/r, dy/r);
	//	//std::cout<<"\n"<<vp.angIndex;
	//}
}


//// ############################################################################################################################################################################
void inline ftkVoting_3D::updateDirection_prob(VPoint3D& vp)
{
	//ftkBins2D::iterator max_point;
	//double maxvote = 0;
	//for( int tt=_voteDirec[vp.angIndex].first; tt!=_voteDirec[vp.angIndex].second; tt=(tt+1)%ntheta ) // Update direction for each voting point
	//{
	//	ftkCone2D::iterator from = _conesPru[tt].begin();
	//	ftkCone2D::iterator to = _conesPru[tt].end();
	//	//_conesPru[tt].vote(p,vp);

	//	VotingDirType_3D::PixelType * votingSumArray = _votingSumVotes->GetBufferPointer()+vp.pos;
	//	for(ftkCone2D::iterator it = from; it != to; it++)  // Iterate to the cones
	//	{
	//		for(ftkBins2D::iterator it1=it->begin(); it1!=it->end(); it1++) // Iterate to the bins of each cone
	//		{
	//			if ( votingSumArray[it1->off]>maxvote) 
	//			{
	//				maxvote = max(maxvote, votingSumArray[it1->off]);
	//				max_point = it1;
	//			}    
	//		}
	//	} 
	//}

	//if (maxvote<epsilon)
	//	return;

	//double dx = max_point->x;
	//double dy = max_point->y;

	//double r = sqrt(dx*dx+dy*dy);
	//if (r>epsilon) {
	//	vp.xc = vp.x + max_point->x;
	//	vp.yc = vp.y + max_point->y;
	//	vp.angIndex = computeAngleIndex(dx/r, dy/r);
	//	//std::cout<<"\n"<<vp.angIndex;
	//}
}


//// ############################################################################################################################################################################
void ftkVoting_3D::computeCones(int hmin, int hmax, int radius)
{
	std::cout << std::endl << "Computing the cones in the classic way";
	clock_t begin_1=clock();
	_conesPru_3D_new = std::vector < ftkCone3D >(_NN_dir); // Todos los conos posibles

	double dlong = pi*(3-sqrt(5.0));
	double dz = 2.0/_NN_dir;
	double long_ = 0;
	double z = 1-dz/2;

	for( int uu=0; uu<_NN_dir; uu++ ) // Phi
	{
		double r = sqrt(1-z*z);
		_conesPru_3D_new[uu].dxx = cos(long_)*r;
		_conesPru_3D_new[uu].dyy = sin(long_)*r;
		_conesPru_3D_new[uu].dzz = z;
		z    = z - dz;
		long_ = long_ + dlong;

		for( int uuu=0; uuu<hmax-hmin+1; uuu++ ) 
		{
			ftkBins3D bin;
			_conesPru_3D_new[uu].push_back(bin);
			// Direction of evenly distributed directions
		}
	}

	int z1_new;
	int R_new;
	double R_dou_new;
	int R_quan_new;
	double x_nor_new, y_nor_new, z_nor_new;
	//int maxmax_pos_ = 0;

	ftkWPoint3D wp_new;
	for( int xx=-hmax; xx<=hmax; ++xx )
	{
		for( int yy=-hmax; yy<=hmax; ++yy )
		{
			z1_new = (int)floor((double)sqrt((double)hmax*hmax-(double)xx*xx-(double)yy*yy));
			for( int zz=-z1_new; zz<=z1_new; ++zz )
			{

				R_dou_new = (double)sqrt((double)zz*zz+(double)yy*yy+(double)xx*xx);
				R_new = (int)ceil(R_dou_new);

				if( R_new>=hmin ) // Can be done more efficiently but for now is ok the speed
				{
					R_quan_new = R_new-hmin;
					//std::cout<<std::endl<<R_quan_new;
					x_nor_new = ((double)xx)/R_dou_new;
					y_nor_new = ((double)yy)/R_dou_new;
					z_nor_new = ((double)zz)/R_dou_new;

					wp_new.x = xx;
					wp_new.y = yy;
					wp_new.z = zz;
					wp_new.w = 1; // In case of deciding to put some weight

					double maxmax_dot = 0;
					int maxmax_pos = 0;

					for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
					{
						double dotprod = _conesPru_3D_new[uu].dxx*x_nor_new+_conesPru_3D_new[uu].dyy*y_nor_new+_conesPru_3D_new[uu].dzz*z_nor_new;
						if ( maxmax_dot<dotprod)
						{
							maxmax_dot = dotprod;
							maxmax_pos = uu;
						}
					}
					_conesPru_3D_new[maxmax_pos][R_quan_new].push_back(wp_new);
				}
			}
		}
	}

	_voteDirec_3D_new = std::vector< std::vector< std::vector< int > > >(_NN_dir,std::vector< std::vector <int> >(10));
	double max_dotpro_delta = atan((double)_radius/(double)_hmax)/10; // Max span

	//INEFFICIENT, 
	//#pragma omp parallel for 
	for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	{
		for( int uu2=0; uu2<_NN_dir; ++uu2 ) // All possible directions
		{
			double dotproduct = acos(_conesPru_3D_new[uu].dxx*_conesPru_3D_new[uu2].dxx+_conesPru_3D_new[uu].dyy*_conesPru_3D_new[uu2].dyy+_conesPru_3D_new[uu].dzz*_conesPru_3D_new[uu2].dzz);

			int angle = ceil(dotproduct/max_dotpro_delta);
			if( (0<=angle) && (angle<10) )
			{
				_voteDirec_3D_new[uu][angle].push_back(uu2);
			}
		}
	}
	std::cout << " -> Done";
	clock_t end_1=clock();
	std::cout << "Time elapsed: " << double(nftkVotingGlobal::diffclock(end_1,begin_1)) << " s";




	// Parellel to calculate the (pixels) on a given direction

	std::cout << std::endl << "Computing the cones in the parallel way";
	clock_t begin_2=clock();

	_conesPru_3D_new_para = std::vector < ftkCone3D >(_NN_dir); // Todos los conos posibles

	double dlong_para = pi*(3-sqrt(5.0));
	double dz_para = 2.0/_NN_dir;
	double long__para = 0;
	double z_para = 1-dz_para/2;

	for( int uu=0; uu<_NN_dir; uu++ ) // Phi
	{
		double r_para = sqrt(1-z_para*z_para);
		_conesPru_3D_new_para[uu].dxx = cos(long__para)*r_para;
		_conesPru_3D_new_para[uu].dyy = sin(long__para)*r_para;
		_conesPru_3D_new_para[uu].dzz = z_para;
		z_para    = z_para - dz_para;
		long__para = long__para + dlong_para;

		for( int uuu=0; uuu<hmax-hmin+1; uuu++ ) 
		{
			ftkBins3D bin_para;
			_conesPru_3D_new_para[uu].push_back(bin_para);
			// Direction of evenly distributed directions
		}
	}

#pragma omp parallel for private(z1_new,R_new,R_dou_new,R_quan_new,x_nor_new, y_nor_new, z_nor_new)
	for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	{
		double maxmax_dot = 0;
		int maxmax_pos = 0;
		double dotprod = 0;

		for( int xx=-hmax; xx<=hmax; ++xx )
		{
			for( int yy=-hmax; yy<=hmax; ++yy )
			{
				z1_new = (int)floor((double)sqrt((double)hmax*hmax-(double)xx*xx-(double)yy*yy));
				for( int zz=-z1_new; zz<=z1_new; ++zz )
				{
					R_dou_new = (double)sqrt((double)zz*zz+(double)yy*yy+(double)xx*xx);
					R_new = (int)ceil(R_dou_new);

					if( R_new>=hmin ) // Can be done more efficiently but for now is ok the speed
					{
						//std::cout<<std::endl<<R_quan_new;
						x_nor_new = ((double)xx)/R_dou_new;
						y_nor_new = ((double)yy)/R_dou_new;
						z_nor_new = ((double)zz)/R_dou_new;

						//for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
						dotprod = _conesPru_3D_new_para[uu].dxx*x_nor_new+_conesPru_3D_new_para[uu].dyy*y_nor_new+_conesPru_3D_new_para[uu].dzz*z_nor_new;
						if ( maxmax_dot<dotprod)
						{
							maxmax_dot = dotprod;
							maxmax_pos = uu;
							R_quan_new = R_new-hmin;

							ftkWPoint3D wp_new;
							wp_new.x = xx;
							wp_new.y = yy;
							wp_new.z = zz;
							wp_new.w = 1; // In case of deciding to put some weight
						}
					}
				}
			}
		}
		_conesPru_3D_new_para[maxmax_pos][R_quan_new].push_back(wp_new);
	}

	_voteDirec_3D_new_para = std::vector< std::vector< std::vector< int > > >(_NN_dir,std::vector< std::vector <int> >(10));

	double max_dotpro_delta_para = atan((double)_radius/(double)_hmax)/10; // Max span


	//INEFFICIENT, 
	//#pragma omp parallel for 
	for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	{
		for( int uu2=0; uu2<_NN_dir; ++uu2 ) // All possible directions
		{
			double dotproduct = acos(_conesPru_3D_new_para[uu].dxx*_conesPru_3D_new_para[uu2].dxx+_conesPru_3D_new_para[uu].dyy*_conesPru_3D_new_para[uu2].dyy+_conesPru_3D_new_para[uu].dzz*_conesPru_3D_new_para[uu2].dzz);

			int angle = ceil(dotproduct/max_dotpro_delta_para);
			if( (0<=angle) && (angle<10) )
			{
				_voteDirec_3D_new_para[uu][angle].push_back(uu2);
			}
		}
	}
	std::cout << " -> Done. ";
	clock_t end_2=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(end_2,begin_2)) << " s";





	// Parellel to calculate the (pixels) on a given direction

	std::cout << std::endl << "Computing the cones in the parallel way without qua";
	clock_t begin_3=clock();

	_conesPru_3D_new_para_notqu = std::vector < ftkCone3D_new >(_NN_dir); // Todos los conos posibles

	double dlong_para_notqu = pi*(3-sqrt(5.0));
	double dz_para_notqu = 2.0/_NN_dir;
	double long__para_notqu = 0;
	double z_para_notqu = 1-dz_para_notqu/2;

	for( int uu=0; uu<_NN_dir; uu++ ) // Phi
	{
		double r_para_notqu = sqrt(1-z_para_notqu*z_para_notqu);
		_conesPru_3D_new_para_notqu[uu].dxx = cos(long__para_notqu)*r_para_notqu;
		_conesPru_3D_new_para_notqu[uu].dyy = sin(long__para_notqu)*r_para_notqu;
		_conesPru_3D_new_para_notqu[uu].dzz = z_para_notqu;
		z_para_notqu    = z_para_notqu - dz_para_notqu;
		long__para_notqu = long__para_notqu + dlong_para_notqu;

		for( int uuu=0; uuu<hmax-hmin+1; uuu++ ) 
		{
			ftkWPoint3D bin_para_notqu;
			_conesPru_3D_new_para_notqu[uu].push_back(bin_para_notqu);
			// Direction of evenly distributed directions
		}
	}

#pragma omp parallel for private(z1_new,R_new,R_dou_new,R_quan_new,x_nor_new, y_nor_new, z_nor_new)
	for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	{
		double maxmax_dot = 0;
		int maxmax_pos = 0;
		double dotprod = 0;

		for( int xx=-hmax; xx<=hmax; ++xx )
		{
			for( int yy=-hmax; yy<=hmax; ++yy )
			{
				z1_new = (int)floor((double)sqrt((double)hmax*hmax-(double)xx*xx-(double)yy*yy));
				for( int zz=-z1_new; zz<=z1_new; ++zz )
				{
					R_dou_new = (double)sqrt((double)zz*zz+(double)yy*yy+(double)xx*xx);
					R_new = (int)ceil(R_dou_new);

					if( R_new>=hmin ) // Can be done more efficiently but for now is ok the speed
					{
						//std::cout<<std::endl<<R_quan_new;
						x_nor_new = ((double)xx)/R_dou_new;
						y_nor_new = ((double)yy)/R_dou_new;
						z_nor_new = ((double)zz)/R_dou_new;

						//for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
						dotprod = _conesPru_3D_new_para[uu].dxx*x_nor_new+_conesPru_3D_new_para[uu].dyy*y_nor_new+_conesPru_3D_new_para[uu].dzz*z_nor_new;
						if ( maxmax_dot<dotprod)
						{
							maxmax_dot = dotprod;
							maxmax_pos = uu;
							R_quan_new = R_new-hmin;

							ftkWPoint3D wp_new;
							wp_new.x = xx;
							wp_new.y = yy;
							wp_new.z = zz;
							wp_new.w = 1; // In case of deciding to put some weight
						}
					}
				}
			}
		}
		_conesPru_3D_new_para_notqu[maxmax_pos].push_back(wp_new);
	}

	_voteDirec_3D_new_para_notqu = std::vector< std::vector< std::vector< int > > >(_NN_dir,std::vector< std::vector <int> >(10));
	double max_dotpro_delta_para_notqu = atan((double)_radius/(double)_hmax)/10; // Max span

	//INEFFICIENT, 
	//#pragma omp parallel for 
	for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	{
		for( int uu2=0; uu2<_NN_dir; ++uu2 ) // All possible directions
		{
			double dotproduct = acos(_conesPru_3D_new_para_notqu[uu].dxx*_conesPru_3D_new_para_notqu[uu2].dxx+_conesPru_3D_new_para_notqu[uu].dyy*_conesPru_3D_new_para_notqu[uu2].dyy+_conesPru_3D_new_para_notqu[uu].dzz*_conesPru_3D_new_para_notqu[uu2].dzz);

			int angle = ceil(dotproduct/max_dotpro_delta_para_notqu);
			if( (0<=angle) && (angle<10) )
			{
				_voteDirec_3D_new_para_notqu[uu][angle].push_back(uu2);
			}
		}
	}
	std::cout << " -> Done. ";
	clock_t end_3=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(end_3,begin_3)) << " s";








	//TO TEST THE SIZE OF THE VECTORS
	//for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	//{
	//	for( int uu_angle=0; uu_angle<10; ++uu_angle)
	//	{
	//		std::cout<<std::endl<<uu<<" "<<uu_angle<<" "<<_voteDirec_3D_new[uu][uu_angle].size();


	//	}
	//}



	//for( int uu=0; uu<ntheta/2; uu++ ) // Phi // REFERENCE DIRECTION
	//{
	//	for( int uu2=0; uu2<ntheta; uu2++ ) // Theta
	//	{
	//		int temp_pri =uu2;

	//		for( int vv=0; vv<ntheta/2; vv++ ) // Phi
	//		{
	//			for( int vv2=0; vv2<ntheta; vv2++ ) // Theta
	//			{			
	//				double dotproduct = _conesPru_3D[vv][vv2].dxx*_conesPru_3D[uu][temp_pri].dxx+_conesPru_3D[vv][vv2].dyy*_conesPru_3D[uu][temp_pri].dyy+_conesPru_3D[vv][vv2].dzz*_conesPru_3D[uu][temp_pri].dzz;

	//				double angle = floor(dotproduct/max_dotpro_delta);
	//				if( (0<=angle) && (angle<10) )
	//				{
	//					pair<int,int> wpp;
	//					wpp.first = vv;
	//					wpp.second = vv2;
	//					_voteDirec_3D_theta[uu][temp_pri][angle].push_back(wpp);
	//				}
	//			}
	//		}
	//	}
	//}






















	////_conesPru = std::vector<std::vector<ftkCone3D>(ntheta)>(ntheta/2); // Todos los conos posibles
	//_conesPru_3D = std::vector < std::vector < ftkCone3D > >(ntheta/2, std::vector<ftkCone3D>(ntheta)); // Todos los conos posibles
	//for( int uu=0; uu<ntheta/2; uu++ ) // Phi
	//{
	//	for( int uu2=0; uu2<ntheta; uu2++ ) // Theta
	//	{
	//		for( int uuu=0; uuu<hmax-hmin+1; uuu++ ) 
	//		{
	//			ftkBins3D bin;
	//			_conesPru_3D[uu][uu2].push_back(bin);
	//		}
	//	}
	//}


	//int z1;
	//int R;
	//double R_dou;
	//int R_quan;
	//pair<int,int> ang_quan;
	//double x_nor, y_nor, z_nor;
	//double x_nor2, y_nor2;
	//int countt=0;

	//ftkWPoint3D wp;
	//for( int xx=-hmax; xx<=hmax; ++xx )
	//{
	//	for( int yy=-hmax; yy<=hmax; ++yy )
	//	{
	//		z1 = (int)floor((double)sqrt((double)hmax*hmax-(double)xx*xx-(double)yy*yy));
	//		for( int zz=0; zz<=z1; ++zz )
	//		{
	//			R_dou = (double)sqrt((double)zz*zz+(double)yy*yy+(double)xx*xx);
	//			R = (int)ceil(R_dou);

	//			if( R>=hmin ) // Can be done more efficiently but for now is ok the speed
	//			{

	//				R_quan = R-hmin;
	//				//std::cout<<std::endl<<R_quan;
	//				x_nor = ((double)xx)/R_dou;
	//				y_nor = ((double)yy)/R_dou;
	//				z_nor = ((double)zz)/R_dou;
	//				ang_quan = computeAngleIndex_3D(x_nor, y_nor, z_nor);
	//				//ftkWPoint2D wp;
	//				wp.x = xx;
	//				wp.y = yy;
	//				wp.z = zz;
	//				wp.w = 1; // In case of deciding to put some weight
	//				_conesPru_3D[ang_quan.first][ang_quan.second][R_quan].push_back(wp);

	//			}
	//		}
	//	}
	//}

	//// Put the direction of the center pixel dxx, dyy, dzz
	//for( int uu=0; uu<ntheta/2; uu++ ) // Phi
	//{
	//	for( int uu2=0; uu2<ntheta; uu2++ ) // Theta
	//	{
	//		_conesPru_3D[uu][uu2].dxx = sin(uu*delta_theta)*cos(uu2*delta_theta);
	//		_conesPru_3D[uu][uu2].dyy = sin(uu*delta_theta)*sin(uu2*delta_theta);
	//		_conesPru_3D[uu][uu2].dzz = cos(uu*delta_theta);
	//	}
	//}








	// Now calculate the span
	//double span = sqrt((double)hmax*hmax+(double)radius*radius);
	//double dx_dou = hmax/span;
	//double dy_dou = radius/span;
	//_intSpan = computeAngleIndex(dx_dou, dy_dou);
	//_voteDirec = std::vector< pair< int,int > > (ntheta);


	////std::vector < std::vector < std::vector<int> > > _voteDirec_3D_phi = std::vector < std::vector < std::vector<int> > >(ntheta/2, std::vector< std::vector<int> >(ntheta) );
	//std::vector < std::vector < std::vector< std::vector< pair<int,int> > > > > _voteDirec_3D_theta = std::vector < std::vector < std::vector< std::vector< pair<int,int> > > > >(ntheta/2, std::vector< std::vector< std::vector< pair<int,int> > > >(ntheta, std::vector< std::vector< pair<int,int> > >(10) ) );

	//radius = 10;
	//double span = sqrt((double)hmax*hmax+(double)radius*radius);
	//double max_dotpro = hmax/sqrt((double)hmax*hmax+(double)radius*radius);
	//double max_dotpro_delta = max_dotpro/(10);

	//	std::cout<<std::endl<<"START";
	//
	////#pragma omp parallel for
	//	for( int uu=0; uu<ntheta/2; uu++ ) // Phi // REFERENCE DIRECTION
	//	{
	//		for( int uu2=0; uu2<ntheta; uu2++ ) // Theta
	//		{
	//			int temp_pri =uu2;
	//
	//			for( int vv=0; vv<ntheta/2; vv++ ) // Phi
	//			{
	//				for( int vv2=0; vv2<ntheta; vv2++ ) // Theta
	//				{			
	//					double dotproduct = _conesPru_3D[vv][vv2].dxx*_conesPru_3D[uu][temp_pri].dxx+_conesPru_3D[vv][vv2].dyy*_conesPru_3D[uu][temp_pri].dyy+_conesPru_3D[vv][vv2].dzz*_conesPru_3D[uu][temp_pri].dzz;
	//
	//					double angle = floor(dotproduct/max_dotpro_delta);
	//					if( (0<=angle) && (angle<10) )
	//					{
	//						pair<int,int> wpp;
	//						wpp.first = vv;
	//						wpp.second = vv2;
	//						_voteDirec_3D_theta[uu][temp_pri][angle].push_back(wpp);
	//					}
	//				}
	//			}
	//		}
	//	}
	//
	//	std::cout<<std::endl<<"END";
	//	std::cin>>rert;

	//for( int uu=0; uu<ntheta/2; uu++ ) // Phi // REFERENCE DIRECTION
	//{
	//	for( int uu2=0; uu2<ntheta; uu2++ ) // Theta
	//	{
	//		for( int uu3=0; uu3<10; ++uu3 )
	//		{
	//			std::cout<<std::endl<<uu<<" "<<uu2<<" "<<uu3<<": "<<_voteDirec_3D_theta[uu][uu2][uu3].size();
	//		}
	//	}
	//}

















	//// CAREFULL WITH THIS 30
	//int numberofcercles = 30; 
	//_conesPru = std::vector<ftkCone3D>(numberofcercles); // Todos los conos posibles
	//for( int uu=0; uu<numberofcercles; uu++ )
	//{
	//	for( int uuu=0; uuu<hmax-hmin+1; uuu++ )
	//	{
	//		ftkBins3D bin;
	//		_conesPru[uu].push_back(bin);
	//	}
	//}

	//int z1;
	//int R;
	//double R_dou;
	//int R_quan;
	//pair<int,int> ang_quan;
	//double x_nor, y_nor, z_nor;
	//double x_nor2, y_nor2;
	//int countt=0;

	//ftkWPoint3D wp;
	//for( int xx=0; xx<=hmax; xx++ )
	//{
	//	for( int yy=0; yy<=hmax; yy++ )
	//	{
	//		z1 = (int)floor((double)sqrt((double)hmax*hmax-(double)xx*xx-(double)yy*yy));
	//		for( int zz=1; zz<=z1; zz++ )
	//		{
	//			R_dou = (double)sqrt((double)zz*zz+(double)yy*yy+(double)xx*xx);
	//			R = (int)ceil(R_dou);

	//			if( R>= hmin) // Can be done more efficiently but for now is ok the speed
	//			{
	//				
	//				R_quan = R-hmin;
	//				//std::cout<<std::endl<<R_quan;
	//				x_nor = ((double)xx)/R_dou;
	//				y_nor = ((double)yy)/R_dou;
	//				z_nor = ((double)zz)/R_dou;
	//				ang_quan = computeAngleIndex_3D(x_nor, y_nor, z_nor);
	//				//ftkWPoint2D wp;
	//				wp.x = xx;
	//				wp.y = yy;
	//				wp.z = zz;
	//				wp.w = 1; // In case of deciding to put some weight
	//				if( ang_quan.second < numberofcercles )
	//				{
	//					_conesPru[ang_quan.second][R_quan].push_back(wp);
	//				}

	//			x_nor2 = -x_nor;
	//			y_nor2 = y_nor;
	//			ang_quan = computeAngleIndex_3D(x_nor2, y_nor2, z_nor); // It is like ang_quan + 128, + 256 etc...
	//			wp.x = -xx;
	//			wp.y = yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//				if( ang_quan.second < numberofcercles )
	//				{
	//					_conesPru[ang_quan.second][R_quan].push_back(wp);
	//				}

	//			x_nor2 = x_nor;
	//			y_nor2 = -y_nor;
	//			ang_quan = computeAngleIndex_3D(x_nor2, y_nor2, z_nor); // It is like ang_quan + 128, + 256 etc...
	//			wp.x = xx;
	//			wp.y = -yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//				if( ang_quan.second < numberofcercles )
	//				{
	//					_conesPru[ang_quan.second][R_quan].push_back(wp);
	//				}

	//			x_nor2 = -x_nor;
	//			y_nor2 = -y_nor;
	//			ang_quan = computeAngleIndex_3D(x_nor2, y_nor2, z_nor); // It is like ang_quan + 128, + 256 etc...
	//			wp.x = -xx;
	//			wp.y = -yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//				if( ang_quan.second < numberofcercles )
	//				{
	//					_conesPru[ang_quan.second][R_quan].push_back(wp);
	//				}



	//			}
	//		}
	//	}
	//}

	//for( int xx=hmin; xx<=hmax; xx++ )
	//{

	//		R_dou = xx;
	//		R_quan = R_dou-hmin;
	//		x_nor = ((double)xx)/R_dou;
	//		y_nor = 0;
	//		ang_quan = computeAngleIndex_3D(x_nor, y_nor, z_nor);
	//		//ftkWPoint2D wp;
	//		wp.x = 0;
	//		wp.y = 0;
	//		wp.z = xx;
	//		wp.w = 1; // In case of deciding to put some weight
	//		if( ang_quan.second < numberofcercles )
	//		{
	//			_conesPru[ang_quan.second][R_quan].push_back(wp);
	//		}
	//}

	////countt = 0;
	//for( int tt=0; tt<_conesPru.size(); tt++ )
	//{
	//	_conesPru[tt].setOffset();
	//}






	//// Now calculate the span
	//double span = sqrt((double)hmax*hmax+(double)radius*radius);
	//double dx_dou = hmax/span;
	//double dy_dou = radius/span;
	//_intSpan = computeAngleIndex(dx_dou, dy_dou);
	//_voteDirec = std::vector< pair< int,int > > (ntheta);
	//for( int tt=0; tt<ntheta; tt++ )
	//{
	//	pair<int,int> spann;
	//	if(tt-_intSpan<0)
	//	{
	//		spann.first = (tt-_intSpan+ntheta)%ntheta;
	//	}
	//	else
	//	{
	//		spann.first = (tt-_intSpan)%ntheta;
	//	}
	//	//std::cout<<std::endl<<"\t"<<spann.first;
	//	spann.second = (tt+_intSpan)%ntheta;
	//	_voteDirec[tt] = spann;
	//}


	//std::cout<<"Cones computed";
	////std::cin>>stop11;






	//int stop11;

	//_conesPru = std::vector<ftkCone2D>(ntheta); // Todos los conos posibles
	//for( int uu=0; uu<ntheta; uu++ )
	//{
	//	for( int uuu=0; uuu<hmax-hmin+1; uuu++ )
	//	{
	//		ftkBins2D bin;
	//		_conesPru[uu].push_back(bin);
	//	}
	//}

	//int y1;
	//int R;
	//double R_dou;
	//int R_quan;
	//int ang_quan;
	//double x_nor, y_nor;
	//double x_nor2, y_nor2;
	//int countt=0;

	//ftkWPoint2D wp;
	//#pragma omp parallel for num_threads(8) private(y1)
	//for( int xx=1; xx<=hmax; xx++ )
	//{
	//	y1 = (int)floor((double)sqrt((double)hmax*hmax-(double)xx*xx));
	//	#pragma omp parallel for num_threads(8) private(wp,R_dou,R,R_quan,x_nor,y_nor,ang_quan)
	//	for( int yy=1; yy<=y1; yy++ )
	//	{
	//		R_dou = (double)sqrt((double)yy*yy+(double)xx*xx);
	//		R = (int)ceil(R_dou);
	//		if( R>= hmin) // Can be done more efficiently but for now is ok the speed
	//		{
	//			
	//			R_quan = R-hmin;
	//			//std::cout<<std::endl<<R_quan;
	//			x_nor = ((double)xx)/R_dou;
	//			y_nor = ((double)yy)/R_dou;
	//			ang_quan = computeAngleIndex(x_nor, y_nor);
	//			//ftkWPoint2D wp;
	//			wp.x = xx;
	//			wp.y = yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//			//std::cout<<std::endl<<ang_quan<<" "<<R_quan;
	//			_conesPru[ang_quan][R_quan].push_back(wp);

	//			x_nor2 = -x_nor;
	//			y_nor2 = y_nor;
	//			ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 128, + 256 etc...
	//			wp.x = -xx;
	//			wp.y = yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//			_conesPru[ang_quan][R_quan].push_back(wp);

	//			x_nor2 = x_nor;
	//			y_nor2 = -y_nor;
	//			ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 128, + 256 etc...
	//			wp.x = xx;
	//			wp.y = -yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//			_conesPru[ang_quan][R_quan].push_back(wp);

	//			x_nor2 = -x_nor;
	//			y_nor2 = -y_nor;
	//			ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 128, + 256 etc...
	//			wp.x = -xx;
	//			wp.y = -yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//			_conesPru[ang_quan][R_quan].push_back(wp);

	//			//std::cout<<std::endl<<"coutn: "<<countt++;

	//			//std::cout<<std::endl<<" "<<R_quan<<" "<<ang_quan<<" "<<R_dou<<"-> "<<xx<<","<<yy<<"|"<<-xx<<" "<<yy<<"|"<<xx<<" "<<-yy<<"|"<<-xx<<" "<<-yy;
	//		}
	//	}
	//}
	////ftkWPoint2D wp;
	//#pragma omp parallel for num_threads(16) private(wp,R_dou,R_quan,x_nor,y_nor,ang_quan)
	//for( int xx=hmin; xx<=hmax; xx++ )
	//{

	//		R_dou = xx;
	//		R_quan = R_dou-hmin;
	//		x_nor = ((double)xx)/R_dou;
	//		y_nor = 0;
	//		ang_quan = computeAngleIndex(x_nor, y_nor);
	//		//ftkWPoint2D wp;
	//		wp.x = xx;
	//		wp.y = 0;
	//		wp.w = 1; // In case of deciding to put some weight
	//		_conesPru[ang_quan][R_quan].push_back(wp);

	//		x_nor = -((double)xx)/R_dou;
	//		y_nor = 0;
	//		ang_quan = computeAngleIndex(x_nor, y_nor);
	//		wp.x = -xx;
	//		wp.y = 0;
	//		wp.w = 1; // In case of deciding to put some weight
	//		_conesPru[ang_quan][R_quan].push_back(wp);

	//		x_nor = 0;
	//		y_nor = ((double)xx)/R_dou;
	//		ang_quan = computeAngleIndex(x_nor, y_nor);
	//		wp.x = 0;
	//		wp.y = xx;
	//		wp.w = 1; // In case of deciding to put some weight
	//		_conesPru[ang_quan][R_quan].push_back(wp);

	//		x_nor = 0;
	//		y_nor = -((double)xx)/R_dou;
	//		ang_quan = computeAngleIndex(x_nor, y_nor);
	//		wp.x = 0;
	//		wp.y = -xx;
	//		wp.w = 1; // In case of deciding to put some weight
	//		_conesPru[ang_quan][R_quan].push_back(wp);


	//		//std::cout<<std::endl<<"coutn: "<<countt++;

	//		//std::cout<<std::endl<<" "<<R_quan<<" "<<ang_quan<<" "<<R_dou<<"-> "<<xx;

	//}


	////std::cout<<"Terminio: ";
	////int ght;
	////std::cin>>ght;


	////countt = 0;
	//for( int tt=0; tt<_conesPru.size(); tt++ )
	//{
	//	_conesPru[tt].setOffset();
	//}

	//// Now calculate the span
	//double span = sqrt((double)hmax*hmax+(double)radius*radius);
	//double dx_dou = hmax/span;
	//double dy_dou = radius/span;
	//_intSpan = computeAngleIndex(dx_dou, dy_dou);
	//_voteDirec = std::vector< pair< int,int > > (ntheta);
	//for( int tt=0; tt<ntheta; tt++ )
	//{
	//	pair<int,int> spann;
	//	if(tt-_intSpan<0)
	//	{
	//		spann.first = (tt-_intSpan+ntheta)%ntheta;
	//	}
	//	else
	//	{
	//		spann.first = (tt-_intSpan)%ntheta;
	//	}
	//	//std::cout<<std::endl<<"\t"<<spann.first;
	//	spann.second = (tt+_intSpan)%ntheta;
	//	_voteDirec[tt] = spann;
	//}


	//std::cout<<"Cones computed";
	////std::cin>>stop11;

}


//// ############################################################################################################################################################################
void ftkVoting_3D::computeCones_prob(int hmin, int hmax, int radius)
{

	//int stop11;

	//_conesPru_prob = std::vector<ftkCone2D>(ntheta); // Todos los conos posibles
	//for( int uu=0; uu<ntheta; uu++ )
	//{
	//	for( int uuu=0; uuu<hmax-hmin+1; uuu++ )
	//	{
	//		ftkBins2D bin;
	//		_conesPru_prob[uu].push_back(bin);
	//	}
	//}

	//int y1;
	//int R;
	//double R_dou;
	//int R_quan;
	//int ang_quan;
	//double x_nor, y_nor;
	//double x_nor2, y_nor2;
	//int countt=0;

	//ftkWPoint2D wp;
	////#pragma omp parallel for num_threads(8) private(y1)
	//for( int xx=1; xx<=hmax; xx++ )
	//{
	//	y1 = (int)floor((double)sqrt((double)hmax*hmax-(double)xx*xx));
	//	//#pragma omp parallel for num_threads(8) private(wp,R_dou,R,R_quan,x_nor,y_nor,ang_quan)
	//	for( int yy=1; yy<=y1; yy++ )
	//	{
	//		R_dou = (double)sqrt((double)yy*yy+(double)xx*xx);
	//		R = (int)ceil(R_dou);
	//		if( R>= hmin) // Can be done more efficiently but for now is ok the speed
	//		{
	//			
	//			R_quan = R-hmin;
	//			//std::cout<<std::endl<<R_quan;
	//			x_nor = ((double)xx)/R_dou;
	//			y_nor = ((double)yy)/R_dou;
	//			ang_quan = computeAngleIndex(x_nor, y_nor);
	//			//ftkWPoint2D wp;
	//			wp.x = xx;
	//			wp.y = yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//			//std::cout<<std::endl<<ang_quan<<" "<<R_quan;
	//			_conesPru_prob[ang_quan][R_quan].push_back(wp);

	//			x_nor2 = -x_nor;
	//			y_nor2 = y_nor;
	//			ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 128, + 256 etc...
	//			wp.x = -xx;
	//			wp.y = yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//			_conesPru_prob[ang_quan][R_quan].push_back(wp);

	//			x_nor2 = x_nor;
	//			y_nor2 = -y_nor;
	//			ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 128, + 256 etc...
	//			wp.x = xx;
	//			wp.y = -yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//			_conesPru_prob[ang_quan][R_quan].push_back(wp);

	//			x_nor2 = -x_nor;
	//			y_nor2 = -y_nor;
	//			ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 128, + 256 etc...
	//			wp.x = -xx;
	//			wp.y = -yy;
	//			wp.w = 1; // In case of deciding to put some weight
	//			_conesPru_prob[ang_quan][R_quan].push_back(wp);

	//			//std::cout<<std::endl<<"coutn: "<<countt++<<" BBB";

	//			//std::cout<<std::endl<<" "<<R<<" "<<R_quan<<" "<<ang_quan<<" "<<R_dou<<"-> "<<xx<<","<<yy<<"|"<<-xx<<" "<<yy<<"|"<<xx<<" "<<-yy<<"|"<<-xx<<" "<<-yy<<" ACA";
	//		}
	//	}
	//}
	//std::cout<<std::endl<<"gggg";
	////ftkWPoint2D wp;
	////#pragma omp parallel for num_threads(16) private(wp,R_dou,R_quan,x_nor,y_nor,ang_quan)
	//for( int xx=hmin; xx<=hmax; xx++ )
	//{

	//		R_dou = xx;
	//		R_quan = R_dou-hmin;
	//		x_nor = ((double)xx)/R_dou;
	//		y_nor = 0;
	//		ang_quan = computeAngleIndex(x_nor, y_nor);
	//		//ftkWPoint2D wp;
	//		wp.x = xx;
	//		wp.y = 0;
	//		wp.w = 1; // In case of deciding to put some weight
	//		_conesPru_prob[ang_quan][R_quan].push_back(wp);

	//		x_nor = -((double)xx)/R_dou;
	//		y_nor = 0;
	//		ang_quan = computeAngleIndex(x_nor, y_nor);
	//		wp.x = -xx;
	//		wp.y = 0;
	//		wp.w = 1; // In case of deciding to put some weight
	//		_conesPru_prob[ang_quan][R_quan].push_back(wp);

	//		x_nor = 0;
	//		y_nor = ((double)xx)/R_dou;
	//		ang_quan = computeAngleIndex(x_nor, y_nor);
	//		wp.x = 0;
	//		wp.y = xx;
	//		wp.w = 1; // In case of deciding to put some weight
	//		_conesPru_prob[ang_quan][R_quan].push_back(wp);

	//		x_nor = 0;
	//		y_nor = -((double)xx)/R_dou;
	//		ang_quan = computeAngleIndex(x_nor, y_nor);
	//		wp.x = 0;
	//		wp.y = -xx;
	//		wp.w = 1; // In case of deciding to put some weight
	//		_conesPru_prob[ang_quan][R_quan].push_back(wp);

	//		//std::cout<<std::endl<<"coutn: "<<countt++;

	//		//std::cout<<std::endl<<" "<<R_quan<<" "<<ang_quan<<" "<<R_dou<<"-> "<<xx;

	//}


	//std::cout<<std::endl<<"Computing probability";

	////for( int tt=0; tt<_conesPru_prob.size(); tt++ )
	////{
	////	std::cout<<std::endl<<"tamanos: "<<tt<<" "<<_conesPru_prob[tt].size()<<" -> ";
	////	for( int uu=0; uu<_conesPru_prob[tt].size(); ++uu )
	////	{
	////		std::cout<<" "<<_conesPru_prob[tt][uu].size();
	////	}
	////	
	////}

	////////////////////////////////////////////////
	////std::cout<<"Terminio de computar conos de probabilidad: ";
	////int ght;
	////std::cin>>ght;
	////////////////////////////////////////////////


	////countt = 0;
	//for( int tt=0; tt<_conesPru_prob.size(); tt++ )
	//{
	//	_conesPru_prob[tt].setOffset();
	//}

	//// Now calculate the span
	//double span = sqrt((double)hmax*hmax+(double)radius*radius);
	//double dx_dou = hmax/span;
	//double dy_dou = radius/span;
	//_intSpan = computeAngleIndex(dx_dou, dy_dou);
	//_voteDirec_prob = std::vector< pair< pair< int,int >, std::vector <int> > > (ntheta);
	//for( int tt=0; tt<ntheta; tt++ ) // 0 256
	//{
	//	pair<int,int> spann;
	//	if(tt-_intSpan<0)
	//	{
	//		spann.first = (tt-_intSpan+ntheta)%ntheta;
	//	}
	//	else
	//	{
	//		spann.first = (tt-_intSpan)%ntheta;
	//	}
	//	spann.second = (tt+_intSpan)%ntheta;

	//	// calculate how many points at disance 
	//	std::vector< int > vi_temp_1 (hmax-hmin+1,0); //hmax-hmin+1 integers of value 0
	//	for( int yy=spann.first; yy!=spann.second; yy=(yy+1)%ntheta )
	//	{
	//		for( int uu=0; uu<hmax-hmin+1; ++uu )
	//		{
	//			vi_temp_1[uu]  += _conesPru_prob[yy][uu].size();
	//		}
	//	}

	//	pair< pair< int,int >, std::vector<int> > spann_prob;
	//	spann_prob.first = spann;
	//	spann_prob.second = vi_temp_1;
	//	_voteDirec_prob[tt] = spann_prob;

	//}

	////for( int tt=0; tt<ntheta; tt++ ) // 0 256
	////{
	////	std::cout<<std::endl<<"Tam: "<<tt<<", ";
	////	for( int uu=0; uu<_voteDirec_prob[tt].second.size(); ++uu )
	////	{
	////		std::cout<<_voteDirec_prob[tt].second[uu]<<" ";
	////	}
	////}

	//std::cout<<"Cones computed";
	////std::cin>>stop11;

}


// ############################################################################################################################################################################
//int ftkVoting_3D::nextConeRadius(int h, int r)
//{
//	//  return  r-1;
//	if (h<=0 || r<0) {
//		cerr << "Error in ftkVoting_3D::nextConeRadius: null cone\n";
//		return 0;
//	}
//	return (int)((double) h * tan(atan((double)r/(double)h) * 0.7 )); // reduce the angle by 30% each time
//}

//// ############################################################################################################################################################################
void ftkVoting_3D::vote()
{

	std::cout << std::endl << "--> VOTING STARTED <--";

	std::cout << std::endl << "Memory Allocation .";
	_t1_begin=clock();
	// Allocate Image of sum of votes
	_votingSumVotes = VotingDirType_3D::New();
	VotingDirType_3D::IndexType start;
	start[0] = 0;
	start[1] = 0;
	start[2] = 0;
	VotingDirType_3D::SizeType size;
	size[0] = nx;
	size[1] = ny;
	size[2] = nz;
	VotingDirType_3D::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	_votingSumVotes->SetRegions( region );
	_votingSumVotes->Allocate();
	const VotingDirType_3D::PixelType ceros = 0;
	_votingSumVotes->FillBuffer( ceros );
	try{
		_votingSumVotes->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}

	// Allocate Image of Mask of votes
	_votingMaskVotes = VotingDirType_3D::New(); 
	_votingMaskVotes->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	_votingMaskVotes->Allocate();
	_votingMaskVotes->FillBuffer( ceros );
	try{
		_votingMaskVotes->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}


	//// Creates a vectorial image, that stores the histogram of voting form each pixel, too much memory, not used for now :(
	//ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
	//for( int oo=0; oo<256; oo++ )
	//{
	//	VotingDirPerType_scalar::Pointer image_1 = VotingDirPerType_scalar::New();
	//	VotingDirPerType_scalar::IndexType start_1;
	//	start_1[0] = 0;
	//	start_1[1] = 0;
	//	VotingDirPerType_scalar::SizeType size_1;
	//	size_1[0] = nx;
	//	size_1[1] = ny;
	//	VotingDirPerType_scalar::RegionType region_1;
	//	region_1.SetSize(size_1);
	//	region_1.SetIndex(start_1);
	//	image_1->SetRegions(region_1);
	//	image_1->Allocate();
	//	const VotingDirPerType_scalar::PixelType ceros_1 = 0;
	//	image_1->FillBuffer(ceros_1);
	//	image_1->Update();
	//	
	//	imageToVectorImageFilter->SetInput(oo, image_1);
	//}
	//imageToVectorImageFilter->Update();
	//VotingDirPerType::Pointer votingDirPerPixel = imageToVectorImageFilter->GetOutput();






	//std::vector < std::vector < int > > votingMaskVotes_dir = std::vector< std::vector< int > >(nx*ny);
	//for( int yr=0; yr<nx*ny; yr++ )
	//{
	//	std::vector< int > bindedir = std::vector< int >(256,0); //256 int of value 0
	//	votingMaskVotes_dir.push_back( bindedir );
	//}

	// Alocate Image of voting pixels
	VotingDirType_3D::Pointer imageOfVotingPixels = VotingDirType_3D::New(); 
	imageOfVotingPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageOfVotingPixels->Allocate();
	imageOfVotingPixels->FillBuffer( ceros );
	try{
		imageOfVotingPixels->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}

	// Allocate Imafe of cones Pixels
	VotingDirType_3D::Pointer imageOfConexPixels = VotingDirType_3D::New(); 
	imageOfConexPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageOfConexPixels->Allocate();
	imageOfConexPixels->FillBuffer( ceros );
	try{
		imageOfConexPixels->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}

	// Allocate Image of Derivative in X
	VotingDirType_3D::Pointer imageGradXPixels = VotingDirType_3D::New(); 
	imageGradXPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageGradXPixels->Allocate();
	imageGradXPixels->FillBuffer( ceros );
	try{
		imageGradXPixels->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}

	// Allocate Image of Derivative in Y
	VotingDirType_3D::Pointer imageGradYPixels = VotingDirType_3D::New(); 
	imageGradYPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageGradYPixels->Allocate();
	imageGradYPixels->FillBuffer( ceros );
	try{
		imageGradYPixels->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}

	// Allocate Image of Max Value
	VotingDirType_3D::Pointer imageMaxPixels = VotingDirType_3D::New(); 
	imageMaxPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageMaxPixels->Allocate();
	imageMaxPixels->FillBuffer( ceros );
	try{
		imageMaxPixels->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}

	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(_t1_end,_t1_begin)) << " s";

	VotingDirType_3D::PixelType * votingSumArray = _votingSumVotes->GetBufferPointer();
	VotingDirType_3D::PixelType * votingMaskArray = _votingMaskVotes->GetBufferPointer();
	VotingDirType_3D::PixelType * imageOfVotingPixelsArray = imageOfVotingPixels->GetBufferPointer();
	VotingDirType_3D::PixelType * imageOfConexPixelsArray = imageOfConexPixels->GetBufferPointer();
	VotingDirType_3D::PixelType * imageGradXPixelsArray = imageGradXPixels->GetBufferPointer();
	VotingDirType_3D::PixelType * imageGradYPixelsArray = imageGradYPixels->GetBufferPointer();
	VotingDirType_3D::PixelType * imageMaxPixelsArray = imageMaxPixels->GetBufferPointer();


	// Update the voting points to the new size (paddng)
	std::vector<VPoint3D>::iterator voting_points_begin = _voting_points.begin();
	std::vector<VPoint3D>::iterator voting_points_end = _voting_points.end();

	for(std::vector<VPoint3D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
		it->pos = it->x+nx*it->y+nx*ny*it->z; // Poscion donde esta x,y en nuestra recien creada _sum (Image)
	}

	std::cout << std::endl << "Update the voting points to the new size .";
	_t1_begin=clock();
	int temp = _voting_points_3D.size();
#pragma omp parallel for
	for( int i=0; i<temp; ++i )
	{
		_voting_points_3D.at(i).pos = _voting_points_3D.at(i).x + _voting_points_3D.at(i).y*nx + _voting_points_3D.at(i).z*nx*ny;
	}
	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVotingGlobal::diffclock(_t1_end,_t1_begin)) << " s";



	// One Vote just to test









	// PARALLEL
	clock_t begin33=clock();
#pragma omp parallel for
	for( int uu=0; uu<_NN_dir; ++uu )
	{
		for( int raddd2=0; raddd2<60; ++raddd2 )
		{

			for( int bin_cont2=0; bin_cont2 < _conesPru_3D_new[uu][raddd2].size(); ++bin_cont2 )
			{
				int x_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].x+100;
				int y_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].y+100;
				int z_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].z+100;
				votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
			}
		}
	}
	std::cout << std::endl << "Done drawing a crcle, Para, ";
	clock_t end33=clock();
	std::cout << "Time elapsed First: " << double(nftkVotingGlobal::diffclock(end33,begin33)) << " s";
	nftkVotingGlobal::stopProgram();

	// PARALLEL 2
	clock_t begin37=clock();
	for( int uu=0; uu<_NN_dir; ++uu )
	{
#pragma omp parallel for
		for( int raddd2=0; raddd2<60; ++raddd2 )
		{
			for( int bin_cont2=0; bin_cont2 < _conesPru_3D_new[uu][raddd2].size(); ++bin_cont2 )
			{
				int x_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].x+100;
				int y_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].y+100;
				int z_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].z+100;
				votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
			}
		}
	}
	std::cout << std::endl << "Done drawing a crcle, Para 2, ";
	clock_t end37=clock();
	std::cout << "Time elapsed First: " << double(nftkVotingGlobal::diffclock(end37,begin37)) << " s";
	nftkVotingGlobal::stopProgram();

	// PARALLEL 4
	clock_t begin39=clock();
	for( int uu=0; uu<_NN_dir; ++uu )
	{
#pragma omp parallel for
		for( int raddd2=0; raddd2<60; ++raddd2 )
		{
			for( int bin_cont2=0; bin_cont2 < _conesPru_3D_new[uu][raddd2].size(); ++bin_cont2 )
			{
				int x_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].x+100;
				int y_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].y+100;
				int z_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].z+100;
				votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
				if( (uu>=_NN_dir) || (raddd2>=60) || (bin_cont2>=_conesPru_3D_new[uu][raddd2].size()) )
				{
					std::cout << std::endl << "houston we have a problem: ";
				}
			}
		}
	}
	std::cout << std::endl << "Done drawing a crcle, Para 4, ";
	clock_t end39=clock();
	std::cout << "Time elapsed First: " << double(nftkVotingGlobal::diffclock(end39,begin39)) << " s";
	nftkVotingGlobal::stopProgram();

	// PARALLEL 3
	clock_t begin38=clock();
	for( int uu=0; uu<_NN_dir; ++uu )
	{
#pragma omp parallel for
		for( int raddd2=0; raddd2<60; ++raddd2 )
		{
			for( int bin_cont2=0; bin_cont2 < _conesPru_3D_new[uu][raddd2].size(); ++bin_cont2 )
			{
				int x_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].x+100;
				int y_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].y+100;
				int z_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].z+100;
				votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
				if( (uu>=_NN_dir) || (raddd2>=60) || (bin_cont2>=_conesPru_3D_new[uu][raddd2].size()) )
				{
					std::cout << std::endl << "houston we have a problem: ";
				}
			}
		}
	}
	std::cout << std::endl << "Done drawing a crcle, Para 3, ";
	clock_t end38=clock();
	std::cout << "Time elapsed First: " << double(nftkVotingGlobal::diffclock(end38,begin38)) << " s";
	nftkVotingGlobal::stopProgram();

	// NOT PARALLEL
	clock_t begin34=clock();
	//#pragma omp parallel for private(raddd2,bin_cont2)
	for( int uu=0; uu<_NN_dir; ++uu )
	{
		for( int raddd2=0; raddd2<60; ++raddd2 )
		{
			for(  int bin_cont2=0; bin_cont2 < _conesPru_3D_new[uu][raddd2].size(); ++bin_cont2 )
			{
				int x_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].x+100;
				int y_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].y+100;
				int z_posi = _conesPru_3D_new[uu][raddd2][bin_cont2].z+100;
				votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
			}
		}
	}
	std::cout << std::endl << "Done drawing a crcle, No Para, ";
	clock_t end34=clock();
	std::cout << "Time elapsed First: " << double(nftkVotingGlobal::diffclock(end34,begin34)) << " s";
	nftkVotingGlobal::stopProgram();



	////#pragma omp parallel for
	//for( int uu=0; uu<_NN_dir; ++uu )
	//{
	//	//std::cout<<std::endl<<"UU: "<<uu;
	//	for( unsigned int raddd=0; raddd<60; ++raddd )
	//	{
	//		//std::cout<<std::endl<<"aja 2";
	//		//std::cout<<std::endl<<"INTER: "<<_conesPru_3D_new[uu][raddd].size();
	//		for( unsigned int bin_cont=0; bin_cont < _conesPru_3D_new[uu][raddd].size(); ++bin_cont )
	//		{
	//			int x_posi = _conesPru_3D_new[uu][raddd][bin_cont].x+250;
	//			int y_posi = _conesPru_3D_new[uu][raddd][bin_cont].y+250;
	//			int z_posi = _conesPru_3D_new[uu][raddd][bin_cont].z+250;
	//			votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
	//		}
	//	}
	//	stringstream outas;
	//	outas<<uu;
	//	string sas = outas.str();
	//	string filenameCones_2as = "output\\out_new_votingSumVotes_cone"+sas+".tif";
	//	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//}


	int raddd;
	int bin_cont;

	// PARALLEL

	clock_t begin=clock();
	for( int oo=0; oo<1; ++oo)
	{
		//#pragma omp parallel for
		for( int uu=0; uu<200/*_NN_dir*/; ++uu )
		{
			//std::cout<<std::endl<<"UU_new: "<<oo<<" "<<uu;
			for( int angle_int = 0; angle_int<10; ++angle_int )
			{
				//std::cout<<std::endl<<"SIZE: "<<_voteDirec_3D_new[uu][angle_int].size();
				//// instead use vector iterator
				//				#pragma omp parallel for private(raddd,bin_cont)
				for( int vv = 0; vv<_voteDirec_3D_new[uu][angle_int].size(); ++vv )
				{
					//std::cout<<std::endl<<"VV: "<<vv;
					int temp_vv = _voteDirec_3D_new[uu][angle_int][vv];
					for( raddd=0; raddd<60; ++raddd )
					{
						//std::cout<<std::endl<<"INTER: "<<_conesPru_3D_new[temp_vv][raddd].size();
#pragma omp parallel for
						for( int bin_cont=0; bin_cont < _conesPru_3D_new[temp_vv][raddd].size(); ++bin_cont )
						{
							int x_posi = _conesPru_3D_new[temp_vv][raddd][bin_cont].x+250;
							int y_posi = _conesPru_3D_new[temp_vv][raddd][bin_cont].y+250;
							int z_posi = _conesPru_3D_new[temp_vv][raddd][bin_cont].z+250;
							votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
						}
					}
				}
				//stringstream outas;
				//outas<<angle_int;
				//string sas = outas.str();
				//string filenameCones_2as = "output\\out_new_votingSumVotes_cone"+sas+".tif";
				//if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
				//	std::cout<<std::endl<<"\tProblema escribiendo";
				//}
			}
		}
	}
	clock_t end=clock();
	std::cout << "Time elapsed Second Parallel: " << double(nftkVotingGlobal::diffclock(end,begin)) << " s";
	//nftkVotingGlobal::stopProgram();


	// NOT PARALLEL

	clock_t begin67=clock();
	for( int oo=0; oo<1; ++oo)
	{
		//#pragma omp parallel for
		for( int uu=0; uu<200/*_NN_dir*/; ++uu )
		{
			//std::cout<<std::endl<<"UU_new: "<<oo<<" "<<uu;
			for( int angle_int = 0; angle_int<10; ++angle_int )
			{
				//std::cout<<std::endl<<"SIZE: "<<_voteDirec_3D_new[uu][angle_int].size();
				//// instead use vector iterator
				for( int vv = 0; vv<_voteDirec_3D_new[uu][angle_int].size(); ++vv )
				{
					//std::cout<<std::endl<<"VV: "<<vv;
					int temp_vv = _voteDirec_3D_new[uu][angle_int][vv];
					for( int raddd=0; raddd<60; ++raddd )
					{
						//std::cout<<std::endl<<"INTER: "<<_conesPru_3D_new[temp_vv][raddd].size();
						for( int bin_cont=0; bin_cont < _conesPru_3D_new[temp_vv][raddd].size(); ++bin_cont )
						{
							int x_posi = _conesPru_3D_new[temp_vv][raddd][bin_cont].x+250;
							int y_posi = _conesPru_3D_new[temp_vv][raddd][bin_cont].y+250;
							int z_posi = _conesPru_3D_new[temp_vv][raddd][bin_cont].z+250;
							votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
						}
					}
				}
				//stringstream outas;
				//outas<<angle_int;
				//string sas = outas.str();
				//string filenameCones_2as = "output\\out_new_votingSumVotes_cone"+sas+".tif";
				//if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
				//	std::cout<<std::endl<<"\tProblema escribiendo";
				//}
			}
		}
	}
	clock_t end67=clock();
	std::cout << "Time elapsed Second Parallel: " << double(nftkVotingGlobal::diffclock(end67,begin67)) << " s";
	nftkVotingGlobal::stopProgram();








	//std::vector < std::vector < int > > temp_direc_1;
	//std::vector < int > temp_angle_2;
	//ftkCone3D  temp_direc_cone_1;
	//ftkBins3D  temp_dist_1;
	//int x_posi;
	//int y_posi;
	//int z_posi;


	clock_t begin3=clock();
	//#pragma omp parallel for
	for( int oo=0; oo<10; ++oo)
	{
		std::cout<<std::endl<<"OO: "<<oo;
		//#pragma omp parallel for
		for( int uu=0; uu<_NN_dir; ++uu )
		{
			//std::cout<<std::endl<<"UU_new2: "<<oo<<" "<<uu;
			//std::vector < std::vector < int > > temp_direc_1 = _voteDirec_3D_new[uu];
			std::vector < std::vector < int > > * temp_direc_1 = &(_voteDirec_3D_new[uu]);

			int vv;
			int raddd;
			int bin_cont;

			//#pragma omp parallel for private(vv,raddd,bin_cont)
			for( int angle_int = 0; angle_int<10; ++angle_int )
			{
				//std::cout<<std::endl<<"SIZE: "<<_voteDirec_3D_new[uu][angle_int].size();
				//// instead use vector iterator
				//std::vector < int > temp_angle_2 = temp_direc_1[angle_int];
				std::vector < int > * temp_angle_2 = &((*temp_direc_1)[angle_int]);
				//#pragma omp parallel for private(raddd,bin_cont)
#pragma omp parallel for
				for( int vv = 0; vv<(*temp_angle_2).size(); ++vv )
				{
					//std::cout<<std::endl<<"VV: "<<vv;
					//int temp_direc_cone_1 = temp_angle_2[vv];
					//ftkCone3D temp_direc_cone_1 = _conesPru_3D_new[temp_angle_2[vv]];
					ftkCone3D * temp_direc_cone_1 = &(_conesPru_3D_new[(*temp_angle_2)[vv]]);
					//#pragma omp parallel for private(bin_cont)
					for( int raddd=0; raddd<60; ++raddd )
					{
						//std::cout<<std::endl<<"INTER: "<<_conesPru_3D_new[temp_vv][raddd].size();
						//ftkBins3D temp_dist_1 = temp_direc_cone_1[raddd];
						ftkBins3D  * temp_dist_1 = &((*temp_direc_cone_1)[raddd]);
						//#pragma omp parallel for
						for( int bin_cont=0; bin_cont < (*temp_dist_1).size(); ++bin_cont )
						{
							//int x_posi = temp_dist_1[bin_cont].x+250;
							//int y_posi = temp_dist_1[bin_cont].y+250;
							//int z_posi = temp_dist_1[bin_cont].z+250;
							int x_posi = (*temp_dist_1)[bin_cont].x+250;
							int y_posi = (*temp_dist_1)[bin_cont].y+250;
							int z_posi = (*temp_dist_1)[bin_cont].z+250;
							//						#pragma omp critical
							//						{
							//int ban =0;
							//int pru_lock;
							//do{
							//	pru_lock = votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi];
							//	pru_lock = pru_lock +1;
							//	ban = ban+1;
							//}
							//while(pru_lock != votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi]+1);
							//if( ban==2 )
							//{std::cout<<std::endl<<"collition";
							//}


							//{
							//#pragma omp atomic
							votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
							//}
							//						}
						}
					}
				}
				//stringstream outas;
				//outas<<angle_int;
				//string sas = outas.str();
				//string filenameCones_2as = "output\\out_new_votingSumVotes_cone"+sas+".tif";
				//if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
				//	std::cout<<std::endl<<"\tProblema escribiendo";
				//}
			}
		}
	}

	clock_t end3=clock();
	std::cout << "Time elapsed Third: " << double(nftkVotingGlobal::diffclock(end3,begin3)) << " s";

	nftkVotingGlobal::stopProgram();




	clock_t begin2=clock();


	for( int oo=0; oo<100; ++oo)
	{
		//#pragma omp parallel for
		for( int uu=0; uu<_NN_dir; ++uu )
		{
			//std::cout<<std::endl<<"UU: "<<uu;
			for( int angle_int = 0; angle_int<10; ++angle_int )
			{
				//std::cout<<std::endl<<"SIZE: "<<_voteDirec_3D_new[uu][angle_int].size();
				//// instead use vector iterator
				for( int vv = 0; vv<_voteDirec_3D_new[uu][angle_int].size(); ++vv )
				{
					//std::cout<<std::endl<<"VV: "<<vv;
					int temp_vv = _voteDirec_3D_new[uu][angle_int][vv];

					std::vector< ftkBins3D >::iterator it_dist;
					for( it_dist = _conesPru_3D_new[temp_vv].begin(); it_dist != _conesPru_3D_new[temp_vv].end(); ++it_dist )
						//for( unsigned int raddd=0; raddd<60; ++raddd )
					{
						std::vector< ftkWPoint3D >::iterator it_bin_cont;

						for( it_bin_cont = (*it_dist).begin(); it_bin_cont!=(*it_dist).end(); ++it_bin_cont )
							//for( it_bin_cont = _conesPru_3D_new[temp_vv][raddd].begin(); it_bin_cont<_conesPru_3D_new[temp_vv][raddd].end(); ++it_bin_cont )
						{
							int x_posi = (*it_bin_cont).x+250;
							int y_posi = (*it_bin_cont).y+250;
							int z_posi = (*it_bin_cont).z+250;
							votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
						}
						////std::cout<<std::endl<<"INTER: "<<_conesPru_3D_new[temp_vv][raddd].size();
						//for( unsigned int bin_cont=0; bin_cont < _conesPru_3D_new[temp_vv][raddd].size(); ++bin_cont )1
						//{
						//	int x_posi = _conesPru_3D_new[temp_vv][raddd][bin_cont].x+250;
						//	int y_posi = _conesPru_3D_new[temp_vv][raddd][bin_cont].y+250;
						//	int z_posi = _conesPru_3D_new[temp_vv][raddd][bin_cont].z+250;
						//	votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
						//}
					}
				}
				//stringstream outas;
				//outas<<angle_int;
				//string sas = outas.str();
				//string filenameCones_2as = "output\\out_new_votingSumVotes_cone"+sas+".tif";
				//if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
				//	std::cout<<std::endl<<"\tProblema escribiendo";
				//}
			}
		}
	}

	clock_t end2=clock();
	std::cout << "Time elapsed Fourth: " << double(nftkVotingGlobal::diffclock(end2,begin2)) << " s";

	nftkVotingGlobal::stopProgram();

	std::cout<<std::endl<<"END";
	int rert;
	std::cin>>rert;













	std::string filenameCones_1 = "output\\out__votingSumVotes.tif";
	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_16 >(_votingSumVotes, filenameCones_1.c_str() )){
		std::cout<<std::endl<<"\tProblema escribiendo";
	}
	//#pragma omp parallel for
	//	for( int uu=20; uu<30; ++uu )
	//	{
	//		for( int uu2=10; uu2<40; ++uu2 )
	//		{
	//			//std::cout<<std::endl<<"aja 1";
	//			for( unsigned int raddd=0; raddd<60; ++raddd )
	//			{
	//				//std::cout<<std::endl<<"aja 2";
	//				for( unsigned int bin_cont=0; bin_cont < _conesPru_3D[uu][uu2][raddd].size(); ++bin_cont )
	//				{
	//					int x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+20;
	//					int y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+20;
	//					int z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
	//					votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
	//
	//					//x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+200;
	//					//y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+200;
	//					//z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
	//					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
	//
	//					//x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+250;
	//					//y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+250;
	//					//z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+250;
	//					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;
	//
	//					//std::cout<<std::endl<<"aja";
	//				}
	//			}
	//		}
	//	}

#pragma omp parallel for
	for( int uu=0; uu<32; ++uu )
	{
		for( int uu2=0; uu2<256; ++uu2 )
		{
			//std::cout<<std::endl<<"aja 1";
			for( unsigned int raddd=0; raddd<60; ++raddd )
			{
				//std::cout<<std::endl<<"aja 2";
				for( unsigned int bin_cont=0; bin_cont < _conesPru_3D[uu][uu2][raddd].size(); ++bin_cont )
				{
					//int x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+20;
					//int y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+20;
					//int z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

					int x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+200;
					int y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+200;
					int z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
					votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

					//x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+250;
					//y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+250;
					//z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+250;
					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

					//std::cout<<std::endl<<"aja";
				}
			}
		}
	}

	//#pragma omp parallel for
	for( int uu2=0; uu2<256; ++uu2 )
	{
		std::cout<<std::endl<<uu2;
		for( int uu=0; uu<32; ++uu )
		{
			//std::cout<<std::endl<<"UU: "<<uu;

			//std::cout<<std::endl<<"aja 1";
			for( unsigned int raddd=0; raddd<60; ++raddd )
			{
				//std::cout<<std::endl<<"aja 2";
				for( unsigned int bin_cont=0; bin_cont < _conesPru_3D[uu][uu2][raddd].size(); ++bin_cont )
				{
					//int x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+20;
					//int y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+20;
					//int z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

					//x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+200;
					//y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+200;
					//z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

					int x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+250;
					int y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+250;
					int z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+250;
					votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;



					//std::cout<<std::endl<<"aja";
				}
			}

		}
		std::stringstream outas;
		outas<<uu2;
		std::string sas = outas.str();
		std::string filenameCones_2as = "output\\out__votingSumVotes_cone"+sas+".tif";
		if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
			std::cout<<std::endl<<"\tProblema escribiendo";
		}

		for( int uu=0; uu<32; ++uu )
		{

			//std::cout<<std::endl<<"aja 1";
			for( unsigned int raddd=0; raddd<60; ++raddd )
			{
				//std::cout<<std::endl<<"aja 2";
				for( unsigned int bin_cont=0; bin_cont < _conesPru_3D[uu][uu2][raddd].size(); ++bin_cont )
				{
					//int x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+20;
					//int y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+20;
					//int z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

					//x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+200;
					//y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+200;
					//z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

					int x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+250;
					int y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+250;
					int z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+250;
					votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;



					//std::cout<<std::endl<<"aja";
				}
			}

		}
	}


	//#pragma omp parallel for
	for( int uu2=0; uu2<256; ++uu2 )
	{
		for( int uu=0; uu<32; ++uu )
		{
			std::cout<<std::endl<<"UU: "<<uu2<<" "<<uu;

			//std::cout<<std::endl<<"aja 1";
			for( unsigned int raddd=0; raddd<60; ++raddd )
			{
				//std::cout<<std::endl<<"aja 2";
				for( unsigned int bin_cont=0; bin_cont < _conesPru_3D[uu][uu2][raddd].size(); ++bin_cont )
				{
					//int x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+20;
					//int y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+20;
					//int z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

					//x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+200;
					//y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+200;
					//z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+20;
					//votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

					int x_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].x+100;
					int y_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].y+100;
					int z_posi = _conesPru_3D[uu][uu2][raddd][bin_cont].z+100;
					votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;



					//std::cout<<std::endl<<"aja";
				}
			}

		}
		std::stringstream outas;
		outas<<uu2;
		std::string sas = outas.str();
		std::string filenameCones_2as = "output\\z_out__votingSumVotes_cone"+sas+".tif";
		if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
			std::cout<<std::endl<<"\tProblema escribiendo";
		}

	}


	//#pragma omp parallel for
	for( int yx =0; yx<20; ++yx ){
		for( unsigned int yy =0; yy<20; ++yy ){
			for( unsigned int yz =0; yz<20; ++yz ){

				int x_posi = yx+300;
				int y_posi = yy+250;
				int z_posi = yz+250;
				votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi] = 1;

			}
		}
	}


	std::string filenameCones_2 = "output\\out__votingSumVotes_cone.tif";
	if( nftkVotingGlobal::writeImage< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2.c_str() )){
		std::cout<<std::endl<<"\tProblema escribiendo";
	}





	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//// Primer Voto
	//int nic=0;
	//int count = 0;
	//std::cout<<std::endl<<"Span"<<_intSpan;
	//
	//int ytt;
	//	for( int gg=0; gg<=_hmax-_hmin; gg++ )
	//	{
	//		//std::cout<<std::endl<<"Q pasa "<<gg;
	//		for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//		{
	//			////std::cout<<std::endl<<"\t	Q pasa "<<gg;
	//			votar(votingSumArray+it->pos,*it,it->angIndex,gg);
	//			votar(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//			////votar_dir(votingMaskVotes_dir,votingMaskArray+it->pos,*it,it->angIndex,gg,it->pos);

	//			//votar_prob(votingSumArray+it->pos,*it,it->angIndex,gg);
	//			//votar_prob(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//			
	//		}

	//		// Errase the mask (only to see step by step/by step) mas detallado todavia, ver los picos mas que todo
	//		memset(votingMaskArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//		
	//	}
	//	int rrr;
	//	std::cout<<std::endl<<"Fin primer voto: "<<ftkCone2D::contador_1;
	//	std::cin>>rrr;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////












	//
	//	// Primer Voto
	//	int nic=0;
	//	int count = 0;
	//	std::cout<<std::endl<<"Span"<<_intSpan;
	////	for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) { //it voting points
	////		//if ( count <256 ){
	////
	////		//std::cout<<std::endl<<"\tVoting pixel: "<<nic++<<", "<<it->pos;
	//////		_cones[it->angIndex].vote(votingSumArray+it->pos, *it);
	////		
	////			votar(votingSumArray+it->pos,*it,it->angIndex);
	////
	////
	////
	////		////
	////		//imageOfVotingPixelsArray[it->pos] = 1;
	////		////(*it).x = 1;
	////		////(*it).y = 0;
	////		////_cones[count].vote(imageOfConexPixelsArray+it->pos, *it);
	////		////std::cout<<std::endl<<"Epezando a votar "<<count;
	////		//votar(imageOfConexPixelsArray+it->pos,*it,it->angIndex);
	////		////_conesPru[count].vote(imageOfConexPixelsArray+it->pos, *it,count);
	////
	////
	////		//stringstream out;
	////		//out<<count;
	////		//string s = out.str();
	////		//string filename = "output\\cones\\out_ImageOfVotes_"+s+".jpg";
	////		//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(imageOfConexPixels, filename.c_str() )){
	////		//	std::cout<<std::endl<<"\tProblema escribiendo";
	////		//}
	////		//memset(imageOfConexPixelsArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	////
	////		count++;
	////		//
	////		//}
	////	}
	//
	//	
	//	int ytt;
	//		for( int gg=0; gg<=_hmax-_hmin; gg++ )
	//		{
	//			//std::cout<<std::endl<<"Q pasa "<<gg;
	//			for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//			{
	//				////std::cout<<std::endl<<"\t	Q pasa "<<gg;
	//				votar(votingSumArray+it->pos,*it,it->angIndex,gg);
	//				votar(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//				////votar_dir(votingMaskVotes_dir,votingMaskArray+it->pos,*it,it->angIndex,gg,it->pos);
	//
	//				//votar_prob(votingSumArray+it->pos,*it,it->angIndex,gg);
	//				//votar_prob(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//				
	//			}
	//			stringstream out4;
	//			out4<<_intSpan;
	//			string s4 = out4.str();
	//			stringstream out5;
	//			out5<<_hmax-_hmin-gg;
	//			string s5 = out5.str();
	//			////string filename3 = "output\\out_ImageOfMask_"+s3+".jpg";
	//			////if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
	//			////	std::cout<<std::endl<<"\tProblema escribiendo";
	//			////}
	//
	//			string filename5 = "output\\Sum_stepbystep\\out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename5.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename6 = "output\\Mask_stepbystep\\out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//			//string filename7 = "output\\Sum_stepbystep_mhd\\out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			//if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirType_3D >(_votingSumVotes, filename7.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//			//string filename8 = "output\\Mask_stepbystep_mhd\\out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			//if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirType_3D >(_votingMaskVotes, filename8.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//
	//			// Get the maximum value of step by step
	//			for( int re=0; re<npix; ++re )
	//			{
	//				if( votingMaskArray[re] > imageMaxPixelsArray[re] )
	//				{
	//					imageMaxPixelsArray[re] = votingMaskArray[re];
	//				}
	//			}
	//
	//
	//			string filename9_max = "output\\Max_stepbystep\\out_ImageOfMaxs_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(imageMaxPixels, filename9_max.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//
	//			// Errase the mask (only to see step by step/by step) mas detallado todavia, ver los picos mas que todo
	//			memset(votingMaskArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//
	//			//// The amount of memory required turns out to be to high :( what to do ? 
	//			//string filename9 = "output\\Direction_perpixel\\out_DirPerPixel_"+s4+"__"+s5+".mhd";
	//			//if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirPerType >(votingDirPerPixel, filename9.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//
	//			//string filename9 = "output\\Direction_perpixel\\out_DirPerPixel_"+s4+"__"+s5+".mhd";
	//			//writer->SetInput( votingDirPerPixel );
	//			//writer->SetFileName( filename9.c_str() );
	//			//writer->Update();
	//
	//			//VotingDirType_3D::IndexType acceso;
	//			//acceso[0] = 45;
	//			//acceso[1] = 345;
	//			//std::cout<<std::endl<<" S4: "<<s4<<", S5: "<<s5<<", x: 45, y: 345, "<<_votingSumVotes->GetPixel(acceso);
	//			//std::cout<<std::endl<<"W pasa";
	//			//std::cin>>ytt;
	//
	//
	//			
	//		}
	//		int rrr;
	//		std::cout<<std::endl<<"Fin primer voto: "<<ftkCone2D::contador_1;
	//		//std::cin>>rrr;
	//
	//	//stringstream out44;
	//	//out44<<_intSpan;
	//	//string s44 = out44.str();
	//	//string filename_dir = "output\\Directions\\Dri_vote_"+s44+".txt";
	//	//ofstream dirFilOut;
	//	//dirFilOut.open (filename_dir.c_str());
	//	//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//	//{
	//	//	dirFilOut << it->x <<"\t"<<it->y<<"\t"<<it->angIndex<<"\t"<<it->mag<<"\n";
	//	//}
	//		
	//	//stringstream out2;
	//	//out2<<_intSpan;
	//	//string s2 = out2.str();
	//	//string filename2 = "output\\out_ImageOfSums_"+s2+".jpg";
	//	//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename2.c_str() )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	////votingSumVotes
	//	////Store Derivative in X
	//	//const char* filename = "output\\out_ImageOfVotes.jpg";
	//	//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(imageOfVotingPixels, filename )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	//const char* filename2 = "output\\out_vote_1.jpg";
	//	//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename2 )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	int r = _radius;
	//
	//	while( _intSpan>=0 )
	//	{
	//		std::cout<<std::endl<<"Span"<<_intSpan;
	//		std::cout<<std::endl<<"Updating cones:";
	//		updateCones(); // Reduce _intSpan
	//		std::cout<<"\t done";
	//		
	//
	//		// Missing update direction
	//		int pru = _voting_points.size();
	//		std::vector<VPoint2D>::iterator it=voting_points_begin;
	//		std::cout<<std::endl<<"Updating direction:";
	//#pragma omp parallel for num_threads(16)
	//		for( int ty=0; ty<pru; ty++ ){
	//			updateDirection(*(it+ty));
	//		}
	//		std::cout<<"\t done";
	//		//stringstream out45;
	//		//out45<<_intSpan;
	//		//string s45 = out45.str();
	//		//string filename_dir2 = "output\\Directions\\Dri_vote_"+s45+".txt";
	//		//ofstream dirFilOut2;
	//		//dirFilOut2.open (filename_dir2.c_str());
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//		//{
	//		//	dirFilOut2 << it->x <<"\t"<<it->y<<"\t"<<it->angIndex<<"\t"<<it->mag<<"\n";
	//		//}
	//
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
	//				////getMaxDirection(*it);    
	//				//updateDirection(*it);
	//		//}   
	//
	//		//memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType_3D::PixelType));
	//		//memset(votingSumArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//
	//		//VOTING PERFECTO
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//		//{
	//		//	votar(votingSumArray+it->pos,*it,it->angIndex);
	//		//}
	//		//stringstream out3;
	//		//out3<<_intSpan;
	//		//string s3 = out3.str();
	//		////string filename3 = "output\\out_ImageOfMask_"+s3+".jpg";
	//		////if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
	//		////	std::cout<<std::endl<<"\tProblema escribiendo";
	//		////}
	//		//string filename4 = "output\\out_ImageOfSums_"+s3+".jpg";
	//		//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename4.c_str() )){
	//		//	std::cout<<std::endl<<"\tProblema escribiendo";
	//		//}
	//
	//		memset(votingMaskArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//		std::cout<<std::endl<<"Voting:";
	//		for( int gg=0; gg<=_hmax-_hmin; gg++ )
	//		{
	//			//std::cout<<std::endl<<"Q pasa "<<gg;
	//			for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//			{
	//				//std::cout<<std::endl<<"\t	Q pasa "<<gg;
	//				//votar(votingSumArray+it->pos,*it,it->angIndex,gg);
	//				votar(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//
	//				//votar_prob(votingSumArray+it->pos,*it,it->angIndex,gg);
	//				//votar_prob(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//			}
	//			stringstream out4;
	//			out4<<_intSpan;
	//			string s4 = out4.str();
	//			stringstream out5;
	//			out5<<_hmax-_hmin-gg;
	//			string s5 = out5.str();
	//
	//			string filename5 = "output\\Sum_stepbystep\\out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename5.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename6 = "output\\Mask_stepbystep\\out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//			string filename7 = "output\\Sum_stepbystep_mhd\\out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirType_3D >(_votingSumVotes, filename7.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename8 = "output\\Mask_stepbystep_mhd\\out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirType_3D >(_votingMaskVotes, filename8.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			
	//		}
	//		std::cout<<"\t done";
	//
	//		memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType_3D::PixelType));
	//
	//
	//
	//		if( _intSpan==0)
	//			break;
	//
	//	}
	//	/*memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType_3D::PixelType));*/
	//
	//	//while (r>= 0) 
	//	//{
	//	//	std::cout<<std::endl<<"	r: "<<r;
	//
	//	//	computeCones(_hmin, _hmax, r); // Eror pues los acabamos de calcular
	//
	//	//	for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
	//	//			getMaxDirection(*it);    
	//	//	}    
	//
	//	//	memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType_3D::PixelType));
	//	//	memset(votingSumArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//
	//
	//	//	 //NO ENTIENDO PARA QUE ES ESTO
	//	//	for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) 
	//	//	{
	//	//			_cones[it->angIndex].vote(votingSumArray+it->pos, votingMaskArray+it->pos, *it);
	//	//	}
	//
	//	//	if (r<=0)
	//	//		break;
	//
	//	//	r = nextConeRadius(_hmax, r);
	//
	//	//}
	//	////computeCones(_hmin, _hmax, 0);
	//
	//	//std::cout<<std::endl<<"Vote: End";
	//
	//
	////	_sum = doubleImage(nx, ny);
	////	_mask = doubleImage(nx, ny);
	////
	////	double* pmask = _mask.ptr();
	////	double* psum = _sum.ptr();
	////
	////	std::vector<VPoint2D>::iterator voting_points_begin = _voting_points.begin();
	////	std::vector<VPoint2D>::iterator voting_points_end = _voting_points.end();
	////
	////	for(std::vector<VPoint2D>::iterator it=voting_points_begin; 
	////		it!=voting_points_end; it++) {	
	////			it->pos = _sum.pos(it->x, it->y); // Poscion donde esta x,y en nuestra recien creada _sum (Image)
	////	}
	////
	////	computeCones(_hmin, _hmax, _radius);
	////
	////	for(std::vector<VPoint2D>::iterator it=voting_points_begin; 
	////		it!=voting_points_end; it++) {
	////			_cones[it->angIndex].vote(psum+it->pos, *it);
	////	}
	////
	////	Image* roi;
	////	int loop = 0;
	////	int r = _radius;
	////
	////	while (r>= 0) {
	////
	////		//    std::cout << "Radius = " << r<< std::endl;
	////		computeCones(_hmin, _hmax, r);
	////		// _reg_type = BrightReg;
	////
	////		if (_save_history) {
	////			roi = Roi::getSubImage(_sum, bw, bw, nx-bw-1, ny-bw-1);
	////			string vote_image_number = toString(r);
	////			if (r<10) vote_image_number = "0"+ vote_image_number;
	////			/* comment for saving 16bit ics file instead
	////			writePGM(_prefix+"vote_"+vote_image_number+".pgm",
	////			convertToByteImage(*roi, true));
	////			*/
	////			image_io::writeICS(_prefix+"vote_"+vote_image_number,*roi);
	////			delete roi;
	////		}
	////
	////		for(std::vector<VPoint2D>::iterator it=voting_points_begin; 
	////			it!=voting_points_end; it++) {	
	////				getMaxDirection(*it);    
	////		}    
	////
	////		memcpy(pmask, psum, npix*sizeof(double));
	////		memset(psum, 0,  npix*sizeof(double));
	////
	////		for(std::vector<VPoint2D>::iterator it=voting_points_begin; 
	////			it!=voting_points_end; it++) {
	////				_cones[it->angIndex].vote(psum+it->pos, pmask+it->pos, *it);
	////		}
	////
	////		if (r<=0)
	////			break;
	////
	////		r = nextConeRadius(_hmax, r);
	////
	////	}
	////	computeCones(_hmin, _hmax, 0);
	////}
	////



	//////////////////////////////////////////////////////////////////////////////////////
	// Copia Original 

	//	_votingSumVotes = VotingDirType_3D::New();
	//	VotingDirType_3D::IndexType start;
	//	start[0] = 0;
	//	start[1] = 0;
	//	VotingDirType_3D::SizeType size;
	//	size[0] = nx;
	//	size[1] = ny;
	//	VotingDirType_3D::RegionType region;
	//	region.SetSize( size );
	//	region.SetIndex( start );
	//	_votingSumVotes->SetRegions( region );
	//	_votingSumVotes->Allocate();
	//	const VotingDirType_3D::PixelType ceros = 0;
	//	_votingSumVotes->FillBuffer( ceros );
	//	_votingSumVotes->Update();
	//
	//	_votingMaskVotes = VotingDirType_3D::New(); 
	//	_votingMaskVotes->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	//	_votingMaskVotes->Allocate();
	//	_votingMaskVotes->FillBuffer( ceros );
	//	_votingMaskVotes->Update();
	//
	//
	//	//// Creates a vectorial image, that stores the histogram of voting form each pixel, too much memory, not used for now :(
	//	//ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
	//	//for( int oo=0; oo<256; oo++ )
	//	//{
	//	//	VotingDirPerType_scalar::Pointer image_1 = VotingDirPerType_scalar::New();
	//	//	VotingDirPerType_scalar::IndexType start_1;
	//	//	start_1[0] = 0;
	//	//	start_1[1] = 0;
	//	//	VotingDirPerType_scalar::SizeType size_1;
	//	//	size_1[0] = nx;
	//	//	size_1[1] = ny;
	//	//	VotingDirPerType_scalar::RegionType region_1;
	//	//	region_1.SetSize(size_1);
	//	//	region_1.SetIndex(start_1);
	//	//	image_1->SetRegions(region_1);
	//	//	image_1->Allocate();
	//	//	const VotingDirPerType_scalar::PixelType ceros_1 = 0;
	//	//	image_1->FillBuffer(ceros_1);
	//	//	image_1->Update();
	//	//	
	//	//	imageToVectorImageFilter->SetInput(oo, image_1);
	//	//}
	//	//imageToVectorImageFilter->Update();
	//	//VotingDirPerType::Pointer votingDirPerPixel = imageToVectorImageFilter->GetOutput();
	//
	//
	//
	//
	//
	//
	//	//std::vector < std::vector < int > > votingMaskVotes_dir = std::vector< std::vector< int > >(nx*ny);
	//	//for( int yr=0; yr<nx*ny; yr++ )
	//	//{
	//	//	std::vector< int > bindedir = std::vector< int >(256,0); //256 int of value 0
	//	//	votingMaskVotes_dir.push_back( bindedir );
	//	//}
	//
	//	//
	//	VotingDirType_3D::Pointer imageOfVotingPixels = VotingDirType_3D::New(); 
	//	imageOfVotingPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	//	imageOfVotingPixels->Allocate();
	//	imageOfVotingPixels->FillBuffer( ceros );
	//	imageOfVotingPixels->Update();
	//	//
	//
	//	//
	//	VotingDirType_3D::Pointer imageOfConexPixels = VotingDirType_3D::New(); 
	//	imageOfConexPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	//	imageOfConexPixels->Allocate();
	//	imageOfConexPixels->FillBuffer( ceros );
	//	imageOfConexPixels->Update();
	//	//
	//
	//	VotingDirType_3D::Pointer imageGradXPixels = VotingDirType_3D::New(); 
	//	imageGradXPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	//	imageGradXPixels->Allocate();
	//	imageGradXPixels->FillBuffer( ceros );
	//	imageGradXPixels->Update();
	//	//
	//
	//	VotingDirType_3D::Pointer imageGradYPixels = VotingDirType_3D::New(); 
	//	imageGradYPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	//	imageGradYPixels->Allocate();
	//	imageGradYPixels->FillBuffer( ceros );
	//	imageGradYPixels->Update();
	//	//
	//
	//	VotingDirType_3D::Pointer imageMaxPixels = VotingDirType_3D::New(); 
	//	imageMaxPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	//	imageMaxPixels->Allocate();
	//	imageMaxPixels->FillBuffer( ceros );
	//	imageMaxPixels->Update();
	//	//
	//
	//	VotingDirType_3D::PixelType * votingSumArray = _votingSumVotes->GetBufferPointer();
	//	VotingDirType_3D::PixelType * votingMaskArray = _votingMaskVotes->GetBufferPointer();
	//	VotingDirType_3D::PixelType * imageOfVotingPixelsArray = imageOfVotingPixels->GetBufferPointer();
	//	VotingDirType_3D::PixelType * imageOfConexPixelsArray = imageOfConexPixels->GetBufferPointer();
	//	VotingDirType_3D::PixelType * imageGradXPixelsArray = imageGradXPixels->GetBufferPointer();
	//	VotingDirType_3D::PixelType * imageGradYPixelsArray = imageGradYPixels->GetBufferPointer();
	//	VotingDirType_3D::PixelType * imageMaxPixelsArray = imageMaxPixels->GetBufferPointer();
	//	
	//
	//	std::vector<VPoint2D>::iterator voting_points_begin = _voting_points.begin();
	//	std::vector<VPoint2D>::iterator voting_points_end = _voting_points.end();
	//
	//	for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
	//		it->pos = it->x+nx*it->y; // Poscion donde esta x,y en nuestra recien creada _sum (Image)
	//	}
	//
	//	//std::cout<<std::endl<<"Vote: Conputecone: "<<"Hmin: "<<_hmin<<", Hmax: "<<_hmax;
	//	computeCones(_hmin, _hmax, _radius);
	//	computeCones_prob(_hmin, _hmax, _radius); // Just for one moment
	//
	//	// Primer Voto
	//	int nic=0;
	//	int count = 0;
	//	std::cout<<std::endl<<"Span"<<_intSpan;
	////	for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) { //it voting points
	////		//if ( count <256 ){
	////
	////		//std::cout<<std::endl<<"\tVoting pixel: "<<nic++<<", "<<it->pos;
	//////		_cones[it->angIndex].vote(votingSumArray+it->pos, *it);
	////		
	////			votar(votingSumArray+it->pos,*it,it->angIndex);
	////
	////
	////
	////		////
	////		//imageOfVotingPixelsArray[it->pos] = 1;
	////		////(*it).x = 1;
	////		////(*it).y = 0;
	////		////_cones[count].vote(imageOfConexPixelsArray+it->pos, *it);
	////		////std::cout<<std::endl<<"Epezando a votar "<<count;
	////		//votar(imageOfConexPixelsArray+it->pos,*it,it->angIndex);
	////		////_conesPru[count].vote(imageOfConexPixelsArray+it->pos, *it,count);
	////
	////
	////		//stringstream out;
	////		//out<<count;
	////		//string s = out.str();
	////		//string filename = "output\\cones\\out_ImageOfVotes_"+s+".jpg";
	////		//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(imageOfConexPixels, filename.c_str() )){
	////		//	std::cout<<std::endl<<"\tProblema escribiendo";
	////		//}
	////		//memset(imageOfConexPixelsArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	////
	////		count++;
	////		//
	////		//}
	////	}
	//
	//	
	//	int ytt;
	//		for( int gg=0; gg<=_hmax-_hmin; gg++ )
	//		{
	//			//std::cout<<std::endl<<"Q pasa "<<gg;
	//			for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//			{
	//				////std::cout<<std::endl<<"\t	Q pasa "<<gg;
	//				votar(votingSumArray+it->pos,*it,it->angIndex,gg);
	//				votar(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//				////votar_dir(votingMaskVotes_dir,votingMaskArray+it->pos,*it,it->angIndex,gg,it->pos);
	//
	//				//votar_prob(votingSumArray+it->pos,*it,it->angIndex,gg);
	//				//votar_prob(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//				
	//			}
	//			stringstream out4;
	//			out4<<_intSpan;
	//			string s4 = out4.str();
	//			stringstream out5;
	//			out5<<_hmax-_hmin-gg;
	//			string s5 = out5.str();
	//			////string filename3 = "output\\out_ImageOfMask_"+s3+".jpg";
	//			////if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
	//			////	std::cout<<std::endl<<"\tProblema escribiendo";
	//			////}
	//
	//			string filename5 = "output\\Sum_stepbystep\\out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename5.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename6 = "output\\Mask_stepbystep\\out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//			//string filename7 = "output\\Sum_stepbystep_mhd\\out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			//if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirType_3D >(_votingSumVotes, filename7.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//			//string filename8 = "output\\Mask_stepbystep_mhd\\out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			//if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirType_3D >(_votingMaskVotes, filename8.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//
	//			// Get the maximum value of step by step
	//			for( int re=0; re<npix; ++re )
	//			{
	//				if( votingMaskArray[re] > imageMaxPixelsArray[re] )
	//				{
	//					imageMaxPixelsArray[re] = votingMaskArray[re];
	//				}
	//			}
	//
	//
	//			string filename9_max = "output\\Max_stepbystep\\out_ImageOfMaxs_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(imageMaxPixels, filename9_max.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//
	//			// Errase the mask (only to see step by step/by step) mas detallado todavia, ver los picos mas que todo
	//			memset(votingMaskArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//
	//			//// The amount of memory required turns out to be to high :( what to do ? 
	//			//string filename9 = "output\\Direction_perpixel\\out_DirPerPixel_"+s4+"__"+s5+".mhd";
	//			//if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirPerType >(votingDirPerPixel, filename9.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//
	//			//string filename9 = "output\\Direction_perpixel\\out_DirPerPixel_"+s4+"__"+s5+".mhd";
	//			//writer->SetInput( votingDirPerPixel );
	//			//writer->SetFileName( filename9.c_str() );
	//			//writer->Update();
	//
	//			//VotingDirType_3D::IndexType acceso;
	//			//acceso[0] = 45;
	//			//acceso[1] = 345;
	//			//std::cout<<std::endl<<" S4: "<<s4<<", S5: "<<s5<<", x: 45, y: 345, "<<_votingSumVotes->GetPixel(acceso);
	//			//std::cout<<std::endl<<"W pasa";
	//			//std::cin>>ytt;
	//
	//
	//			
	//		}
	//		int rrr;
	//		std::cout<<std::endl<<"Fin primer voto: "<<ftkCone2D::contador_1;
	//		//std::cin>>rrr;
	//
	//	//stringstream out44;
	//	//out44<<_intSpan;
	//	//string s44 = out44.str();
	//	//string filename_dir = "output\\Directions\\Dri_vote_"+s44+".txt";
	//	//ofstream dirFilOut;
	//	//dirFilOut.open (filename_dir.c_str());
	//	//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//	//{
	//	//	dirFilOut << it->x <<"\t"<<it->y<<"\t"<<it->angIndex<<"\t"<<it->mag<<"\n";
	//	//}
	//		
	//	//stringstream out2;
	//	//out2<<_intSpan;
	//	//string s2 = out2.str();
	//	//string filename2 = "output\\out_ImageOfSums_"+s2+".jpg";
	//	//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename2.c_str() )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	////votingSumVotes
	//	////Store Derivative in X
	//	//const char* filename = "output\\out_ImageOfVotes.jpg";
	//	//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(imageOfVotingPixels, filename )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	//const char* filename2 = "output\\out_vote_1.jpg";
	//	//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename2 )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	int r = _radius;
	//
	//	while( _intSpan>=0 )
	//	{
	//		std::cout<<std::endl<<"Span"<<_intSpan;
	//		std::cout<<std::endl<<"Updating cones:";
	//		updateCones(); // Reduce _intSpan
	//		std::cout<<"\t done";
	//		
	//
	//		// Missing update direction
	//		int pru = _voting_points.size();
	//		std::vector<VPoint2D>::iterator it=voting_points_begin;
	//		std::cout<<std::endl<<"Updating direction:";
	//#pragma omp parallel for num_threads(16)
	//		for( int ty=0; ty<pru; ty++ ){
	//			updateDirection(*(it+ty));
	//		}
	//		std::cout<<"\t done";
	//		//stringstream out45;
	//		//out45<<_intSpan;
	//		//string s45 = out45.str();
	//		//string filename_dir2 = "output\\Directions\\Dri_vote_"+s45+".txt";
	//		//ofstream dirFilOut2;
	//		//dirFilOut2.open (filename_dir2.c_str());
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//		//{
	//		//	dirFilOut2 << it->x <<"\t"<<it->y<<"\t"<<it->angIndex<<"\t"<<it->mag<<"\n";
	//		//}
	//
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
	//				////getMaxDirection(*it);    
	//				//updateDirection(*it);
	//		//}   
	//
	//		//memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType_3D::PixelType));
	//		//memset(votingSumArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//
	//		//VOTING PERFECTO
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//		//{
	//		//	votar(votingSumArray+it->pos,*it,it->angIndex);
	//		//}
	//		//stringstream out3;
	//		//out3<<_intSpan;
	//		//string s3 = out3.str();
	//		////string filename3 = "output\\out_ImageOfMask_"+s3+".jpg";
	//		////if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
	//		////	std::cout<<std::endl<<"\tProblema escribiendo";
	//		////}
	//		//string filename4 = "output\\out_ImageOfSums_"+s3+".jpg";
	//		//if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename4.c_str() )){
	//		//	std::cout<<std::endl<<"\tProblema escribiendo";
	//		//}
	//
	//		memset(votingMaskArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//		std::cout<<std::endl<<"Voting:";
	//		for( int gg=0; gg<=_hmax-_hmin; gg++ )
	//		{
	//			//std::cout<<std::endl<<"Q pasa "<<gg;
	//			for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//			{
	//				//std::cout<<std::endl<<"\t	Q pasa "<<gg;
	//				//votar(votingSumArray+it->pos,*it,it->angIndex,gg);
	//				votar(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//
	//				//votar_prob(votingSumArray+it->pos,*it,it->angIndex,gg);
	//				//votar_prob(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//			}
	//			stringstream out4;
	//			out4<<_intSpan;
	//			string s4 = out4.str();
	//			stringstream out5;
	//			out5<<_hmax-_hmin-gg;
	//			string s5 = out5.str();
	//
	//			string filename5 = "output\\Sum_stepbystep\\out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename5.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename6 = "output\\Mask_stepbystep\\out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVotingGlobal::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//			string filename7 = "output\\Sum_stepbystep_mhd\\out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirType_3D >(_votingSumVotes, filename7.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename8 = "output\\Mask_stepbystep_mhd\\out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			if( nftkVotingGlobal::writeImage_mhdDouble< VotingDirType_3D >(_votingMaskVotes, filename8.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			
	//		}
	//		std::cout<<"\t done";
	//
	//		memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType_3D::PixelType));
	//
	//
	//
	//		if( _intSpan==0)
	//			break;
	//
	//	}
	//	/*memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType_3D::PixelType));*/
	//
	//	//while (r>= 0) 
	//	//{
	//	//	std::cout<<std::endl<<"	r: "<<r;
	//
	//	//	computeCones(_hmin, _hmax, r); // Eror pues los acabamos de calcular
	//
	//	//	for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
	//	//			getMaxDirection(*it);    
	//	//	}    
	//
	//	//	memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType_3D::PixelType));
	//	//	memset(votingSumArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//
	//
	//	//	 //NO ENTIENDO PARA QUE ES ESTO
	//	//	for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) 
	//	//	{
	//	//			_cones[it->angIndex].vote(votingSumArray+it->pos, votingMaskArray+it->pos, *it);
	//	//	}
	//
	//	//	if (r<=0)
	//	//		break;
	//
	//	//	r = nextConeRadius(_hmax, r);
	//
	//	//}
	//	////computeCones(_hmin, _hmax, 0);
	//
	//	//std::cout<<std::endl<<"Vote: End";
	//
	//
	////	_sum = doubleImage(nx, ny);
	////	_mask = doubleImage(nx, ny);
	////
	////	double* pmask = _mask.ptr();
	////	double* psum = _sum.ptr();
	////
	////	std::vector<VPoint2D>::iterator voting_points_begin = _voting_points.begin();
	////	std::vector<VPoint2D>::iterator voting_points_end = _voting_points.end();
	////
	////	for(std::vector<VPoint2D>::iterator it=voting_points_begin; 
	////		it!=voting_points_end; it++) {	
	////			it->pos = _sum.pos(it->x, it->y); // Poscion donde esta x,y en nuestra recien creada _sum (Image)
	////	}
	////
	////	computeCones(_hmin, _hmax, _radius);
	////
	////	for(std::vector<VPoint2D>::iterator it=voting_points_begin; 
	////		it!=voting_points_end; it++) {
	////			_cones[it->angIndex].vote(psum+it->pos, *it);
	////	}
	////
	////	Image* roi;
	////	int loop = 0;
	////	int r = _radius;
	////
	////	while (r>= 0) {
	////
	////		//    std::cout << "Radius = " << r<< std::endl;
	////		computeCones(_hmin, _hmax, r);
	////		// _reg_type = BrightReg;
	////
	////		if (_save_history) {
	////			roi = Roi::getSubImage(_sum, bw, bw, nx-bw-1, ny-bw-1);
	////			string vote_image_number = toString(r);
	////			if (r<10) vote_image_number = "0"+ vote_image_number;
	////			/* comment for saving 16bit ics file instead
	////			writePGM(_prefix+"vote_"+vote_image_number+".pgm",
	////			convertToByteImage(*roi, true));
	////			*/
	////			image_io::writeICS(_prefix+"vote_"+vote_image_number,*roi);
	////			delete roi;
	////		}
	////
	////		for(std::vector<VPoint2D>::iterator it=voting_points_begin; 
	////			it!=voting_points_end; it++) {	
	////				getMaxDirection(*it);    
	////		}    
	////
	////		memcpy(pmask, psum, npix*sizeof(double));
	////		memset(psum, 0,  npix*sizeof(double));
	////
	////		for(std::vector<VPoint2D>::iterator it=voting_points_begin; 
	////			it!=voting_points_end; it++) {
	////				_cones[it->angIndex].vote(psum+it->pos, pmask+it->pos, *it);
	////		}
	////
	////		if (r<=0)
	////			break;
	////
	////		r = nextConeRadius(_hmax, r);
	////
	////	}
	////	computeCones(_hmin, _hmax, 0);
	////}
	////
}

//// ############################################################################################################################################################################
void inline ftkVoting_3D::votar(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx)
{
	////if( angl_indx==1)
	////{
	////	std::cout<<std::endl<<_voteDirec[angl_indx].first<<" "<<_voteDirec[angl_indx].second<<", aca";
	////}
	////std::cout<<std::endl<<"\t angl_indx "<<angl_indx<<","<<_voteDirec[angl_indx].first<<","<<_voteDirec[angl_indx].second;
	//for( int tt=_voteDirec[angl_indx].first; tt!=_voteDirec[angl_indx].second; tt=(tt+1)%ntheta )
	//{
	//	_conesPru[tt].vote(p,vp);
	//}
}

//// ############################################################################################################################################################################
void inline ftkVoting_3D::votar(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx, int dist)
{
	////if( angl_indx==1)
	////{
	////	std::cout<<std::endl<<_voteDirec[angl_indx].first<<" "<<_voteDirec[angl_indx].second<<", aca";
	////}
	////std::cout<<std::endl<<"\t angl_indx "<<angl_indx<<","<<_voteDirec[angl_indx].first<<","<<_voteDirec[angl_indx].second;
	//for( int tt=_voteDirec[angl_indx].first; tt!=_voteDirec[angl_indx].second; tt=(tt+1)%ntheta )
	//{
	//	_conesPru[tt].vote(p,vp,dist);
	//}
}

//// ############################################################################################################################################################################
void inline ftkVoting_3D::votar_prob(VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx, int dist)
{
	////if( angl_indx==1)
	////{
	////	std::cout<<std::endl<<_voteDirec[angl_indx].first<<" "<<_voteDirec[angl_indx].second<<", aca";
	////}
	////std::cout<<std::endl<<"\t angl_indx "<<angl_indx<<","<<_voteDirec[angl_indx].first<<","<<_voteDirec[angl_indx].second;
	//for( int tt=_voteDirec_prob[angl_indx].first.first; tt!=_voteDirec_prob[angl_indx].first.second; tt=(tt+1)%ntheta )
	//{
	//	_conesPru_prob[tt].vote(p,vp,dist,_voteDirec_prob[angl_indx].second[dist]);
	//}
}

//// ############################################################################################################################################################################
void inline ftkVoting_3D::votar_dir(std::vector<std::vector<int> >& p_dir,VotingDirType_3D::PixelType * p, const VPoint3D& vp, int angl_indx, int dist, int& offset_1)
{
	////if( angl_indx==1)
	////{
	////	std::cout<<std::endl<<_voteDirec[angl_indx].first<<" "<<_voteDirec[angl_indx].second<<", aca";
	////}
	////std::cout<<std::endl<<"\t angl_indx "<<angl_indx<<","<<_voteDirec[angl_indx].first<<","<<_voteDirec[angl_indx].second;
	//for( int tt=_voteDirec[angl_indx].first; tt!=_voteDirec[angl_indx].second; tt=(tt+1)%ntheta )
	//{
	//	_conesPru[tt].vote_dir(p_dir,p,vp,dist,offset_1);
	//}
}

// ############################################################################################################################################################################
void ftkVoting_3D::updateCones()
{
	//if(_intSpan==0)
	//{
	//	int rr;
	//	std::cout<<std::endl<<"Error en updatecones";
	//	std::cin>>rr;
	//}

	//int reduction = 2;
	//while(_intSpan-reduction<0) //Si el span se reduce a menos de 0, entonces reduzca hasta que quede en cero
	//{
	//	reduction--;
	//}
	//_intSpan = _intSpan-reduction;
	//int temp;
	////_voteDirec = std::vector< pair< int,int > > (ntheta);
	//for( int tt=0; tt<ntheta; tt++ )
	//{
	//	_voteDirec[tt].first = (_voteDirec[tt].first+reduction)%ntheta;
	//	temp = _voteDirec[tt].second-reduction;
	//	if(temp<0)
	//	{
	//		_voteDirec[tt].second = (temp+ntheta)%ntheta;
	//	}
	//	else
	//	{
	//		_voteDirec[tt].second = (temp)%ntheta;
	//	}
	//}
}

// ############################################################################################################################################################################
int ftkVoting_3D::computeDirecIndex_3D(double dx, double dy, double dz) const
{
	double maxmax_dot = 0;
	int maxmax_pos = 0;
	double dotprod = 0;

	for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	{
		dotprod = _conesPru_3D_new_para[uu].dxx*dx+_conesPru_3D_new_para[uu].dyy*dy+_conesPru_3D_new_para[uu].dzz*dz;
		if ( maxmax_dot<dotprod)
		{
			maxmax_dot = dotprod;
			maxmax_pos = uu;
		}
	}
	return maxmax_pos;
}