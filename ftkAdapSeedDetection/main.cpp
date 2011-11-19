// ############################################################################################################################################################################
// ############################################################################################################################################################################

//#include "itkImage.h"
//#include "itkImageFileReader.h"
//#include "itkImageFileWriter.h"
//#include "itkRescaleIntensityImageFilter.h"
//
//#include "itkDerivativeImageFilter.h"
//#include "itkGradientRecursiveGaussianImageFilter.h"
//#include "itkConstantBoundaryCondition.h"
//#include "itkNeighborhoodAlgorithm.h"
//#include "itkShapedNeighborhoodIterator.h"
//#include "math.h"
//#include "itkSobelOperator.h"
//#include "itkNeighborhoodInnerProduct.h"
//#include "itkRealTimeClock.h"
//#include "itkGradientRecursiveGaussianImageFilter.h"
//
//// TEST
//#include <map>
//#include <utility> // make_pair

#ifdef _OPENMP
#include "omp.h"
#endif

//#include"ftkVoting.h"
//#include"ftkVotingGlobal.h"

#include"ftkVoting_3D.h"
#include"ftkVotingGlobal.h"


// For Extracting slide by slice image of the 3d Image
#include "itkImage.h"
#include "itkExtractImageFilter.h"

//using namespace nftkVotingGlobal;

//using namespace std;


// ############################################################################################################################################################################
// ############################################################################################################################################################################

int main( int argc, char * argv[] ){

	clock_t begin=clock();

	// User Parameters
	char* inputImageName = argv[1];
	char* folder = argv[2];
	//int hmax = atoi(argv[2]);

	// Initial Parameters
	//reg_type (have not underestand what is this for)
	int hmin = 9;
	int hmax = 70;//35;//77; //15 was ok, 25 just to test
	int radius = 70;//10;//77;
	double min_grad = 0.2;
	//threshold (for picking seed points, not necesary for now)
	double scale = 1.5; // Scale for computing the gradient using DoG
	//zc_only TRUE - voting points should also be the zero-crossing points of the image; FALSE - otherwise


	//// ###################################### 2D
	//// Set Up the Reader 2D
	//nftkVotingGlobal::InputImageType::Pointer inputImage = nftkVotingGlobal::readImage< nftkVotingGlobal::InputImageType, nftkVotingGlobal::InputImageType >( inputImageName );
	//inputImage->Update();

	//ftkVoting voteMain;
	//voteMain.setParams(hmin,hmax,radius,min_grad,scale);
	//voteMain.setPrefix("temp/");
	//voteMain.compute(inputImage); // DUDA COMO HACER PARA ENVIAR UN CONST POINTER, DESPUES DE QUE HE LEIDO LA IMAGEN COMO POINTER ??


	// ###################################### 3D
		// Set Up the Reader 2D
	
		// Input Image Type
	
		nftkVotingGlobal::InputImageType_3D::Pointer inputImage = nftkVotingGlobal::readImage_3D< nftkVotingGlobal::InputImageType_3D_16, nftkVotingGlobal::InputImageType_3D >( inputImageName );
		inputImage->Update();
	
		ftkVoting_3D voteMain_3D;
		voteMain_3D.setParams(hmin,hmax,radius,min_grad,scale);
		voteMain_3D.setPrefix("temp/");
		voteMain_3D.compute(inputImage); // DUDA COMO HACER PARA ENVIAR UN CONST POINTER, DESPUES DE QUE HE LEIDO LA IMAGEN COMO POINTER ??

	//// ###################################### 2D Slice by Slice

	//// Set Up the Reader 2D
	//nftkVotingGlobal::InputImageType_3D::Pointer inputImage_3D = nftkVotingGlobal::readImage_3D< nftkVotingGlobal::InputImageType_3D_16, nftkVotingGlobal::InputImageType_3D >( inputImageName );
	//inputImage_3D->Update();


	//double dxx = -0.1;

	//double dyy =0.3;
	//cout<<endl<<" ANGULO: "<<atan(dyy/dxx)<<" "<<atan(dxx/dyy);;

	//typedef itk::ExtractImageFilter< nftkVotingGlobal::InputImageType_3D, nftkVotingGlobal::InputImageType > FilterType;
	//FilterType::Pointer filter = FilterType::New();
	//filter->SetDirectionCollapseToSubmatrix();

	//filter->SetInput( inputImage_3D );


	//nftkVotingGlobal::InputImageType_3D::RegionType inputRegion = inputImage_3D->GetLargestPossibleRegion();
	//nftkVotingGlobal::InputImageType_3D::SizeType size = inputRegion.GetSize();
	//size[2] = 0;

	//nftkVotingGlobal::InputImageType_3D::IndexType start = inputRegion.GetIndex();
	//start[2] = 1;

	//nftkVotingGlobal::InputImageType_3D::RegionType desiredRegion;
	//desiredRegion.SetSize(  size  );
	//desiredRegion.SetIndex( start );

	//filter->SetExtractionRegion( desiredRegion );

	//try
	//{
	//	filter->Update();
	//}
	//catch( itk::ExceptionObject & err )
	//{
	//	std::cerr << "ExceptionObject caught !" << std::endl;
	//	std::cerr << err << std::endl;
	//	return EXIT_FAILURE;
	//}


	////nftkVotingGlobal::InputImageType_3D::Pointer votingBySlice = nftkVotingGlobal::InputImageType_3D::New(); //VotingDirType = double
	////votingBySlice->SetRegions( inputImage_3D->GetRequestedRegion()); // IMPORTANTE PARA CREAR UNA IMAGEN NUEVA EN BASE A UNA QUE YA EXISTE EN VEZ DE PONERME A LEER LOS TAMANOS Y LAS REGIONES
	////votingBySlice->Allocate();

	////nftkVotingGlobal::InputImageType_3D::Pointer votingBySlice = nftkVotingGlobal::readImage_3D< nftkVotingGlobal::InputImageType_3D_16, nftkVotingGlobal::InputImageType_3D >( inputImageName );
	////votingBySlice->Update();
	////int nxxX = votingBySlice->GetLargestPossibleRegion().GetSize()[0];
	////int nyxX = votingBySlice->GetLargestPossibleRegion().GetSize()[1];
	////int nzxX = votingBySlice->GetLargestPossibleRegion().GetSize()[2];
	////int npixX = nxxX*nyxX*nzxX;
	////nftkVotingGlobal::InputImageType_3D::PixelType * votingBySliceArray = votingBySlice->GetBufferPointer();
	////for( unsigned int tr=0; tr<npixX; ++tr ){
	////	votingBySliceArray[0] = 0;
	////}

	//nftkVotingGlobal::InputImageType_3D::Pointer votingBySlice = nftkVotingGlobal::InputImageType_3D::New();
	//nftkVotingGlobal::InputImageType_3D::IndexType start22;
	//start22[0] = 0;
	//start22[1] = 0;
	//start22[2] = 0;
	//nftkVotingGlobal::InputImageType_3D::SizeType size22;
	//double bw = sqrt((double)(radius*radius + hmax*hmax)) + 3; // Que es esto ??
	//const int bw2 = 2*bw;
	//size22[0] = inputImage_3D->GetLargestPossibleRegion().GetSize()[0]+bw2;
	//size22[1] = inputImage_3D->GetLargestPossibleRegion().GetSize()[1]+bw2;
	//size22[2] = inputImage_3D->GetLargestPossibleRegion().GetSize()[2];
	//nftkVotingGlobal::InputImageType_3D::RegionType region;
	//region.SetSize( size22 );
	//region.SetIndex( start22 );
	//votingBySlice->SetRegions( region );
	//votingBySlice->Allocate();
	//const nftkVotingGlobal::InputImageType_3D::PixelType ceros = 0;
	//votingBySlice->FillBuffer( ceros );
	//votingBySlice->Update();

	//nftkVotingGlobal::InputImageType_3D::Pointer votingBySliceLastbigSapan = nftkVotingGlobal::InputImageType_3D::New();
	//votingBySliceLastbigSapan->SetRegions( region );
	//votingBySliceLastbigSapan->Allocate();
	//votingBySliceLastbigSapan->FillBuffer( ceros );
	//votingBySliceLastbigSapan->Update();

	//nftkVotingGlobal::InputImageType_3D::Pointer votingBySliceLastbigSapanProb = nftkVotingGlobal::InputImageType_3D::New();
	//votingBySliceLastbigSapanProb->SetRegions( region );
	//votingBySliceLastbigSapanProb->Allocate();
	//votingBySliceLastbigSapanProb->FillBuffer( ceros );
	//votingBySliceLastbigSapanProb->Update();	




	////int nxxX = votingBySlice->GetLargestPossibleRegion().GetSize()[0];
	////int nyxX = votingBySlice->GetLargestPossibleRegion().GetSize()[1];
	////int nzxX = votingBySlice->GetLargestPossibleRegion().GetSize()[2];
	////int npixX = nxxX*nyxX*nzxX;
	////cout<<endl<<"APSDFASDF: "<<npixX;




	////#pragma omp parallel for
	//for( int uu=0; uu<275; ++uu )
	//{

	//ftkVoting voteMain;

	//voteMain.setOutBySlice( votingBySlice );
	//voteMain.setOutBySliceLastbigSapan( votingBySliceLastbigSapan );
	//voteMain.setOutBySliceLastbigSapanProb( votingBySliceLastbigSapanProb );

	//	cout<<endl<<endl<<"SLICE: "<<uu;
	//	stringstream out;
	//	out<<uu;
	//	string s = out.str();


	//	start[2] = uu;
	//	desiredRegion.SetIndex( start );
	//	filter->SetExtractionRegion( desiredRegion );
	//	filter->Update();

	//	voteMain.SetSlice(uu); // Tiene que empezar en 0
	//	voteMain.setParams(hmin,hmax,radius,min_grad,scale);
	//	voteMain.setPrefix(s);
	//	voteMain.compute(filter->GetOutput()); // DUDA COMO HACER PARA ENVIAR UN CONST POINTER, DESPUES DE QUE HE LEIDO LA IMAGEN COMO POINTER ??

	//}






















	//	char* inputImageName = argv[1];
	//	char* VotingMagName = "z_VMag.tif";//strcat(inputImageName, "VotMag");//argv[2];
	//	char* VotingImaName = "z_VIma.tif";//strcat(inputImageName, "VotIma");//argv[3];
	//	char* votingDirXName = "z_VDirX.tif";//strcat(inputImageName, "VotDirX");//argv[4];
	//	char* votingDirYName = "z_VDirY.tif";//strcat(inputImageName, "VotDirY");//argv[5];
	//	char* votingPixName = "z_VotPixel.tif";//strcat(inputImageName, "VotPix");//argv[6];
	//	double votingPixThre = atof(argv[2]);
	//	int rmin = atoi(argv[3]);
	//	int rmax = atoi(argv[4]);
	//	double delta = atof(argv[5])*3.1416/180;
	//	double maxiter = atof(argv[6]);
	//
	//	double deltaDecrea = (atof(argv[5])-1)/(maxiter-1); // Decrement linear and make the last iteration 1 degree
	//
	//	// Input Image Type
	//	typedef double InputPixelType;
	//	const unsigned int Dimension = 2;
	//	typedef itk::Image< InputPixelType, Dimension > InputImageType;
	//
	//	// Output Image Type
	//	typedef unsigned char OutputPixelType;
	//	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
	//
	//	// Voting Magnitude M
	//	typedef double VotingMagPixelType;
	//	typedef itk::Image< VotingMagPixelType, Dimension > VotingMagType;
	//
	//	// Voting Pixels S
	//	typedef bool VotingPixPixelType;
	//	typedef itk::Image< VotingPixPixelType, Dimension > VotingPixType;
	//
	//	// Voting Direction X and Y
	//	typedef double VotingDirPixelType;
	//	typedef itk::Image< VotingDirPixelType, Dimension > VotingDirType;
	//
	//	// Voting Area A
	//	typedef double VotingArePixelType;
	//	typedef itk::Image< VotingArePixelType, Dimension > VotingAreType;
	//
	//	// Voting Image V
	//	typedef double VotingImaPixelType;
	//	typedef itk::Image< VotingImaPixelType, Dimension > VotingImaType;
	//
	//	// Set Up the Reader
	//	InputImageType::Pointer inputImage = readImage< InputImageType, InputImageType >( inputImageName );
	//
	//	// Derivative XXXXXXXXXX
	//	VotingDirType::Pointer votingDirX = VotingDirType::New();
	//	votingDirX->SetRegions( inputImage->GetRequestedRegion() );
	//	votingDirX->Allocate();
	//
	//	itk::SobelOperator< VotingDirPixelType, 2 > sobelOperator;
	//	sobelOperator.SetDirection( 0 );
	//	sobelOperator.CreateDirectional();
	//	typedef itk::ImageRegionIterator< VotingDirType > IteratorType3;
	//	IteratorType3 out2( votingDirX, votingDirX->GetRequestedRegion() );
	//
	//	typedef itk::ConstNeighborhoodIterator< InputImageType > NeighborhoodIteratorType2;
	//	typedef itk::ImageRegionIterator< InputImageType > IteratorType2;
	//	NeighborhoodIteratorType2::RadiusType radius2 = sobelOperator.GetRadius();
	//	NeighborhoodIteratorType2 it2( radius2, inputImage, inputImage->GetRequestedRegion() );
	//	itk::NeighborhoodInnerProduct< InputImageType > innerProduct;
	//
	//	for (it2.GoToBegin(), out2.GoToBegin(); !it2.IsAtEnd(); ++it2, ++out2)
	//	{
	//		out2.Set( innerProduct( it2, sobelOperator ) );
	//	}
	//
	//	// Derivative YYYYYYYYYY
	//	VotingDirType::Pointer votingDirY = VotingDirType::New();
	//	votingDirY->SetRegions( inputImage->GetRequestedRegion() );
	//	votingDirY->Allocate();
	//
	//	itk::SobelOperator<VotingDirPixelType, 2> sobelOperatorY;
	//	sobelOperatorY.SetDirection( 1 );
	//	sobelOperatorY.CreateDirectional();
	//	typedef itk::ImageRegionIterator< VotingDirType > IteratorType4;
	//	IteratorType4 out3( votingDirY, votingDirY->GetRequestedRegion() );
	//
	//	typedef itk::ConstNeighborhoodIterator< InputImageType > NeighborhoodIteratorType3;
	//	typedef itk::ImageRegionIterator< InputImageType > IteratorType3;
	//	NeighborhoodIteratorType3::RadiusType radius3 = sobelOperatorY.GetRadius();
	//	NeighborhoodIteratorType3 it3( radius3, inputImage, inputImage->GetRequestedRegion() );
	//	itk::NeighborhoodInnerProduct< InputImageType > innerProduct2;
	//
	//	for (it3.GoToBegin(), out3.GoToBegin(); !it3.IsAtEnd(); ++it3, ++out3)
	//	{
	//		out3.Set( innerProduct2( it3, sobelOperatorY ) );
	//	}
	//
	//
	//	// Magnitude Image
	//	VotingMagType::Pointer VotingMagImage = VotingMagType::New();
	//	VotingMagImage->SetRegions( inputImage->GetRequestedRegion()); // IMPORTANTE PARA CREAR UNA IMAGEN NUEVA EN BASE A UNA QUE YA EXISTE EN VEZ DE PONERME A LEER LOS TAMANOS Y LAS REGIONES
	//	VotingMagImage->Allocate();
	//
	//	typedef itk::ImageRegionIteratorWithIndex< VotingMagType > ITVotingMag;
	//	ITVotingMag iVotingMag(VotingMagImage, VotingMagImage->GetLargestPossibleRegion() );
	//
	//	typedef itk::ImageRegionIteratorWithIndex< VotingDirType > ITVotingDir;
	//	ITVotingDir iVotingDirX(votingDirX, votingDirX->GetLargestPossibleRegion() );
	//	ITVotingDir iVotingDirY(votingDirY, votingDirY->GetLargestPossibleRegion() );
	//
	//	for ( iVotingMag.GoToBegin(); !iVotingMag.IsAtEnd(); ++iVotingMag, ++iVotingDirX, ++iVotingDirY ){
	//		iVotingMag.Set(sqrt( iVotingDirX.Get()*iVotingDirX.Get() + iVotingDirY.Get()*iVotingDirY.Get() ));
	//	}
	//
	//	// Voting Direction normalizada
	//	for ( iVotingMag.GoToBegin(), iVotingDirX.GoToBegin(), iVotingDirY.GoToBegin(); !iVotingMag.IsAtEnd(); ++iVotingMag, ++iVotingDirX, ++iVotingDirY){
	//		iVotingDirX.Set(iVotingDirX.Get()/iVotingMag.Get());
	//		iVotingDirY.Set(iVotingDirY.Get()/iVotingMag.Get());
	//	}
	//
	//
	//	// Votes
	//	VotingImaType::Pointer VotingImaImage = VotingImaType::New();
	//	VotingImaImage->SetRegions( inputImage->GetRequestedRegion() );
	//	VotingImaImage->Allocate();
	//	const VotingImaType::PixelType ceros = 0;
	//	VotingImaImage->FillBuffer( ceros );
	//	VotingImaImage->Update();
	//
	//	// Voting Pixels
	//	VotingPixType::Pointer VotingPixImage = VotingPixType::New();
	//	VotingPixImage->SetRegions( inputImage->GetRequestedRegion());
	//	VotingPixImage->Allocate();
	//	const VotingPixType::PixelType ceros2 = 0;
	//	VotingPixImage->FillBuffer( ceros2 );
	//	VotingPixImage->Update();
	//	VotingPixType::PixelType * votingPixArray = VotingPixImage->GetBufferPointer();
	//
	//	// Compute the pixels that actually work
	//	VotingMagType::PixelType * votingMagArray = VotingMagImage->GetBufferPointer();
	//	VotingMagType::SizeType votingMagSize = VotingMagImage->GetLargestPossibleRegion().GetSize();
	//	int votingMagTama = votingMagSize[0]*votingMagSize[1];
	//
	//
	//
	//
	//  ////Filter class is instantiated
	//  //typedef itk::GradientRecursiveGaussianImageFilter<InputImageType, OutputImageType> FilterType;
	// 
	//  //FilterType::Pointer filter = FilterType::New();
	// 
	//  ////sigma is specified in millimeters
	//  //filter->SetSigma( 1.5 );  
	//
	//
	//
	//
	//	pair< int, int > tempCoord;
	//	vector< pair< int, int > > votingPixelVector; // To store voting pixels
	//	int yCoord;
	//	for( int j = 0; j<votingMagSize[1]; ++j ){
	//		tempCoord.second = j;
	//		yCoord = j*votingMagSize[0];
	//		for( int i = 0; i<votingMagSize[0]; ++i ){
	//			if( votingMagArray[yCoord+i]>votingPixThre ){
	//				tempCoord.first = i;
	//				votingPixelVector.push_back(tempCoord);
	//				votingPixArray[yCoord+i] = 1;
	//			}
	//		}
	//	}
	//	cout<<endl<<"Number of votes "<<votingPixelVector.size();
	//
	//
	//	VotingImaType::PixelType * votingImaArray = VotingImaImage->GetBufferPointer();
	//	VotingDirType::PixelType * votingDirXArray = votingDirX->GetBufferPointer();
	//	VotingDirType::PixelType * votingDirYArray = votingDirY->GetBufferPointer();
	//	vector< pair< int, int > >::iterator itVotingPix;
	//
	//	int xCoord;
	//	int yCoord2;
	//	int xCoordVot;
	//	int yCoordVot;
	//	int xyCoord;
	//	//int yCoordVotMul;
	//	int xyCoordCenter;
	//
	//	double dis;
	//	double teta;
	//
	//	int re = 0;
	//
	//	// LISTO FUNCIONANDO
	//	// ------------------------------------------------------------------------------------
	//	//for( itVotingPix = votingPixelVector.begin(); itVotingPix != votingPixelVector.end(); ++itVotingPix )
	//	//{
	//	//	re++;
	//	//	//cout<<endl<<"vota.. ";
	//	//	xCoord = (*itVotingPix).first;
	//	//	yCoord2 = (*itVotingPix).second;
	//	//	xyCoordCenter = yCoord2*votingMagSize[0]+xCoord;
	//	//	//#pragma omp parallel for num_threads(16)
	//	//	for( int yy = -rmax; yy<=rmax; yy++ )
	//	//	{
	//	//		if( (0<=yCoord2+yy) && (yCoord2+yy<votingMagSize[1]) ){
	//	//			//#pragma omp parallel for num_threads(16)
	//	//			for( int xx = -rmax; xx<=rmax; xx++ )
	//	//			{
	//	//				xCoordVot = xCoord+xx;// Coordinate where the vote is going to be
	//	//				if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) ){
	//	//					xyCoord = (yCoord2+yy)*votingMagSize[0]+xCoordVot;
	//	//					dis = sqrt( (double)(xx*xx+yy*yy));
	//	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//	//					if( (rmin<=dis) && (dis<=rmax) && (teta<=delta/2) )
	//	//					{
	//	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//	//					}
	//	//				}
	//	//			}
	//	//		}
	//	//	}
	//	//}
	//	//cout<<endl<<re;
	//
	//
	//
	//	//int const1;
	//	//int const2;
	//
	//	////for( itVotingPix = votingPixelVector.begin(); itVotingPix != votingPixelVector.end(); ++itVotingPix )
	//	////#pragma omp parallel for shared(votingImaArray,votingDirXArray,votingDirYArray,votingMagArray) private(const1,const2,xCoord,yCoord2,xyCoordCenter)
	//	//for( int kk=0; kk<votingPixelVector.size(); kk++ )
	//	//{
	//	//	//xCoord = (*itVotingPix).first;
	//	//	//yCoord2 = (*itVotingPix).second;
	//
	//	//	xCoord = votingPixelVector[kk].first;
	//	//	yCoord2 = votingPixelVector[kk].second;
	//
	//	//	xyCoordCenter = yCoord2*votingMagSize[0]+xCoord;
	//	//	//#pragma omp parallel for num_threads(16) shared(votingImaArray,votingDirXArray,votingDirYArray,votingMagArray) private(const1,const2)
	//	//	for( int yy = -rmax; yy<=rmax; yy++ )
	//	//	{
	//	//		if( (0<=yCoord2+yy) && (yCoord2+yy<votingMagSize[1]) )
	//	//		{
	//	//			//#pragma omp parallel for num_threads(16)
	//	//			const1 = (int)sqrt((double)(rmax*rmax)-(double)(yy*yy));
	//	//			for( int xx = 0; xx<=const1; xx++ )
	//	//			{
	//	//				xCoordVot = xCoord+xx;// Coordinate where the vote is going to be
	//	//				if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) )
	//	//				{
	//	//					xyCoord = (yCoord2+yy)*votingMagSize[0]+xCoordVot;
	//	//					dis = sqrt( (double)(xx*xx+yy*yy));
	//	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//	//					//if( (rmin<=dis) && (dis<=rmax) && (teta<=delta/2) )
	//	//					if( (rmin<=dis) && (teta<=delta/2) )
	//	//					{
	//	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//	//					}
	//	//					//if( dis>rmax){
	//	//					//	cout<<endl<<"Algo fallo: "<<dis;
	//	//					//}
	//	//				}
	//	//			}
	//
	//	//			const2 = -(int)sqrt((double)(rmax*rmax)-(double)(yy*yy));
	//	//			//#pragma omp parallel for num_threads(16)
	//	//			for( int xx = const2; xx<0; xx++ )
	//	//			{
	//	//				xCoordVot = xCoord+xx;// Coordinate where the vote is going to be
	//	//				if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) )
	//	//				{
	//	//					xyCoord = (yCoord2+yy)*votingMagSize[0]+xCoordVot;
	//	//					dis = sqrt( (double)(xx*xx+yy*yy));
	//	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//	//					//if( (rmin<=dis) && (dis<=rmax) && (teta<=delta/2) )
	//	//					if( (rmin<=dis) && (teta<=delta/2) )
	//	//					{
	//	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//	//					}
	//	//					//if( dis>rmax){
	//	//					//	cout<<endl<<"Algo fallo";
	//	//					//}
	//	//				}
	//	//			}
	//
	//	//		}
	//	//	}
	//	//}
	//
	//
	//	int const1;
	//	int const2;
	//	int const3;
	//	int const4;
	//	int const5;
	//	int const6;
	//	int const7;
	//	int const8;
	//	int const9;
	//	int const10;
	//
	//	//pair< int, int > tempCoordCalc;
	//	//vector< pair< int, int > > votingPixelVectorIn;
	//	//vector< vector< pair< int, int > > > votingPixelVectorCalc; // To store voting pixels per pixel
	//
	//	//for( itVotingPix = votingPixelVector.begin(); itVotingPix != votingPixelVector.end(); ++itVotingPix )
	//	for( unsigned int uu=0; uu<maxiter; uu++ )
	//	{
	//		cout<<endl<<"Teration "<<uu+1<<" of "<<maxiter<<", delta: "<<delta;
	//		clock_t partialTimeBegin=clock();
	//#pragma omp parallel for shared(votingImaArray,votingDirXArray,votingDirYArray,votingMagArray) private(const1,const2,const3,const4,const5,const6,const7,const8,const9,const10,xCoord,yCoord2,xyCoordCenter,xCoordVot,yCoordVot)
	//		for( int kk=0; kk<votingPixelVector.size(); kk++ )
	//			//for( int kk=0; kk<1; kk++ )
	//		{
	//			//xCoord = (*itVotingPix).first;
	//			//yCoord2 = (*itVotingPix).second;
	//
	//			xCoord = votingPixelVector[kk].first;
	//			yCoord2 = votingPixelVector[kk].second;
	//
	//			xyCoordCenter = yCoord2*votingMagSize[0]+xCoord;
	//
	//			//if( (0<=xCoord-rmax) && (xCoord+rmax<votingMagSize[0]) && (0<=yCoord2-rmax) && (yCoord2+rmax<votingMagSize[1]) )
	//			//{
	//
	//			const3 = -floor(sqrt( double(rmax*rmax/2) ));
	//			const4 = -ceil(sqrt( double(rmin*rmin/2) ));
	//
	//			//cout<<endl<<const3<<" "<<const4;
	//			for( int yy = const3; yy<=const4; yy++ )
	//			{
	//				const5 = -yy+1;
	//				const6 = floor(sqrt(double(rmax*rmax-yy*yy)));
	//				//cout<<endl<<const5<<" "<<const6;
	//				for( int xx = const5; xx<=const6; xx++ )
	//				{
	//					dis = sqrt( (double)(xx*xx+yy*yy));
	//
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+xx;
	//						yCoordVot = yCoord2+yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-xx;
	//						yCoordVot = yCoord2+yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)(-yy))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-xx;
	//						yCoordVot = yCoord2-yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)(-yy))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+xx;
	//						yCoordVot = yCoord2-yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)yy+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+yy;
	//						yCoordVot = yCoord2+xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-yy)+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-yy;
	//						yCoordVot = yCoord2+xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-yy)+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-yy;
	//						yCoordVot = yCoord2-xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)yy+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+yy;
	//						yCoordVot = yCoord2-xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//
	//
	//
	//					if( (dis<rmin) || (dis>rmax) )
	//					{
	//						cout<<endl<<"falla 333";
	//					}
	//				}
	//			}
	//
	//			const7 = const4+1;
	//			const8 = -1;
	//			//cout<<endl<<const7<<" "<<const8;
	//			for( int yy = const7; yy<=const8; yy++ )
	//			{
	//				const9 = ceil(sqrt(double(rmin*rmin-yy*yy)));
	//				const10 = floor(sqrt(double(rmax*rmax-yy*yy)));
	//				//cout<<endl<<const9<<" "<<const10;
	//				for( int xx = const9; xx<=const10; xx++ )
	//				{
	//					dis = sqrt( (double)(xx*xx+yy*yy));
	//
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+xx;
	//						yCoordVot = yCoord2+yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-xx;
	//						yCoordVot = yCoord2+yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)(-yy))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-xx;
	//						yCoordVot = yCoord2-yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)(-yy))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+xx;
	//						yCoordVot = yCoord2-yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)yy+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+yy;
	//						yCoordVot = yCoord2+xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-yy)+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-yy;
	//						yCoordVot = yCoord2+xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-yy)+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-yy;
	//						yCoordVot = yCoord2-xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)yy+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+yy;
	//						yCoordVot = yCoord2-xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//						}
	//					}
	//
	//					if( (dis<rmin) || (dis>rmax) )
	//					{
	//						cout<<endl<<"falla 333";
	//					}
	//				}
	//			}
	//
	//			// x = xx, y=0
	//			// x = -xx, y=0
	//			// x = 0, y=xx
	//			// x = 0, y=-xx
	//			// x = xx, y=xx
	//			// x = -xx, y=xx
	//			// x = xx, y=-xx
	//			// x = -xx, y=-xx
	//			for( int xx = rmin ;xx<=rmax; xx++ )
	//			{
	//				//dot product
	//				//xx = 0, yy = xx
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord+xx;
	//					yCoordVot = yCoord2+xx;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//						xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//					}
	//				}
	//
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx)/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord+xx;
	//					yCoordVot = yCoord2;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) )
	//					{
	//						xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx))/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord-xx;
	//					yCoordVot = yCoord2;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) )
	//					{
	//						xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//					}
	//				}
	//				teta = abs(acos( (votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord;
	//					yCoordVot = yCoord2+xx;
	//					if( (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//						xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//					}
	//				}
	//				teta = abs(acos( (votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord;
	//					yCoordVot = yCoord2-xx;
	//					if( (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//						xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord+xx;
	//					yCoordVot = yCoord2+xx;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//						xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord-xx;
	//					yCoordVot = yCoord2+xx;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//						xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord+xx;
	//					yCoordVot = yCoord2-xx;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//						xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord-xx;
	//					yCoordVot = yCoord2-xx;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//						xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//						votingImaArray[xyCoord] = votingImaArray[xyCoord]+votingMagArray[xyCoordCenter];
	//					}
	//				}
	//			}
	//		}
	//
	//// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//		// Exact copy of the lines above with change in the internal calculation of each finarl if
	//
	//		double maximImaVot = 0;
	//		int xCoordmaximImaVot = 0;
	//		int yCoordmaximImaVot = 0;
	//		double temp = 0;
	//
	//#pragma omp parallel for shared(votingImaArray,votingDirXArray,votingDirYArray,votingMagArray) private(const1,const2,const3,const4,const5,const6,const7,const8,const9,const10,xCoord,yCoord2,xyCoordCenter,xCoordVot,yCoordVot,xCoordmaximImaVot,yCoordmaximImaVot,maximImaVot,temp)
	//		for( int kk=0; kk<votingPixelVector.size(); kk++ )
	//			//for( int kk=0; kk<1; kk++ )
	//		{
	//			maximImaVot = 0;
	//			//xCoord = (*itVotingPix).first;
	//			//yCoord2 = (*itVotingPix).second;
	//
	//			xCoord = votingPixelVector[kk].first;
	//			yCoord2 = votingPixelVector[kk].second;
	//
	//			xyCoordCenter = yCoord2*votingMagSize[0]+xCoord;
	//
	//			//if( (0<=xCoord-rmax) && (xCoord+rmax<votingMagSize[0]) && (0<=yCoord2-rmax) && (yCoord2+rmax<votingMagSize[1]) )
	//			//{
	//
	//			const3 = -floor(sqrt( double(rmax*rmax/2) ));
	//			const4 = -ceil(sqrt( double(rmin*rmin/2) ));
	//
	//			//cout<<endl<<const3<<" "<<const4;
	//			for( int yy = const3; yy<=const4; yy++ )
	//			{
	//				const5 = -yy+1;
	//				const6 = floor(sqrt(double(rmax*rmax-yy*yy)));
	//				//cout<<endl<<const5<<" "<<const6;
	//				for( int xx = const5; xx<=const6; xx++ )
	//				{
	//					dis = sqrt( (double)(xx*xx+yy*yy));
	//
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+xx;
	//						yCoordVot = yCoord2+yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = xx;
	//								yCoordmaximImaVot = yy;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-xx;
	//						yCoordVot = yCoord2+yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -xx;
	//								yCoordmaximImaVot = yy;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)(-yy))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-xx;
	//						yCoordVot = yCoord2-yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -xx;
	//								yCoordmaximImaVot = -yy;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)(-yy))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+xx;
	//						yCoordVot = yCoord2-yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = xx;
	//								yCoordmaximImaVot = -yy;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)yy+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+yy;
	//						yCoordVot = yCoord2+xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = yy;
	//								yCoordmaximImaVot = xx;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-yy)+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-yy;
	//						yCoordVot = yCoord2+xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -yy;
	//								yCoordmaximImaVot = xx;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-yy)+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-yy;
	//						yCoordVot = yCoord2-xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -yy;
	//								yCoordmaximImaVot = -xx;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)yy+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+yy;
	//						yCoordVot = yCoord2-xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = yy;
	//								yCoordmaximImaVot = -xx;
	//							}
	//						}
	//					}
	//
	//
	//
	//
	//					if( (dis<rmin) || (dis>rmax) )
	//					{
	//						cout<<endl<<"falla 333";
	//					}
	//				}
	//			}
	//
	//			const7 = const4+1;
	//			const8 = -1;
	//			//cout<<endl<<const7<<" "<<const8;
	//			for( int yy = const7; yy<=const8; yy++ )
	//			{
	//				const9 = ceil(sqrt(double(rmin*rmin-yy*yy)));
	//				const10 = floor(sqrt(double(rmax*rmax-yy*yy)));
	//				//cout<<endl<<const9<<" "<<const10;
	//				for( int xx = const9; xx<=const10; xx++ )
	//				{
	//					dis = sqrt( (double)(xx*xx+yy*yy));
	//
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+xx;
	//						yCoordVot = yCoord2+yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = xx;
	//								yCoordmaximImaVot = yy;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)yy)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-xx;
	//						yCoordVot = yCoord2+yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -xx;
	//								yCoordmaximImaVot = yy;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)(-yy))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-xx;
	//						yCoordVot = yCoord2-yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -xx;
	//								yCoordmaximImaVot = -yy;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)(-yy))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+xx;
	//						yCoordVot = yCoord2-yy;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = xx;
	//								yCoordmaximImaVot = -yy;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)yy+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+yy;
	//						yCoordVot = yCoord2+xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = yy;
	//								yCoordmaximImaVot = xx;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-yy)+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-yy;
	//						yCoordVot = yCoord2+xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -yy;
	//								yCoordmaximImaVot = xx;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-yy)+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord-yy;
	//						yCoordVot = yCoord2-xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -yy;
	//								yCoordmaximImaVot = -xx;
	//							}
	//						}
	//					}
	//
	//					teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)yy+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//					if( (teta<=delta/2) )
	//					{
	//						xCoordVot = xCoord+yy;
	//						yCoordVot = yCoord2-xx;
	//						if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//						{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = yy;
	//								yCoordmaximImaVot = -xx;
	//							}
	//						}
	//					}
	//
	//					if( (dis<rmin) || (dis>rmax) )
	//					{
	//						cout<<endl<<"falla 333";
	//					}
	//				}
	//			}
	//
	//			// x = xx, y=0
	//			// x = -xx, y=0
	//			// x = 0, y=xx
	//			// x = 0, y=-xx
	//			// x = xx, y=xx
	//			// x = -xx, y=xx
	//			// x = xx, y=-xx
	//			// x = -xx, y=-xx
	//			for( int xx = rmin ;xx<=rmax; xx++ )
	//			{
	//				//dot product
	//				//xx = 0, yy = xx
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx)/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord+xx;
	//					yCoordVot = yCoord2;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) )
	//					{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = xx;
	//								yCoordmaximImaVot = 0;
	//							}
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx))/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord-xx;
	//					yCoordVot = yCoord2;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) )
	//					{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -xx;
	//								yCoordmaximImaVot = 0;
	//							}
	//					}
	//				}
	//				teta = abs(acos( (votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord;
	//					yCoordVot = yCoord2+xx;
	//					if( (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = 0;
	//								yCoordmaximImaVot = xx;
	//							}
	//					}
	//				}
	//				teta = abs(acos( (votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord;
	//					yCoordVot = yCoord2-xx;
	//					if( (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = 0;
	//								yCoordmaximImaVot = -xx;
	//							}
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord+xx;
	//					yCoordVot = yCoord2+xx;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = xx;
	//								yCoordmaximImaVot = xx;
	//							}
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)xx)/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord-xx;
	//					yCoordVot = yCoord2+xx;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -xx;
	//								yCoordmaximImaVot = xx;
	//							}
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)xx+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord+xx;
	//					yCoordVot = yCoord2-xx;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = xx;
	//								yCoordmaximImaVot = -xx;
	//							}
	//					}
	//				}
	//				teta = abs(acos( (votingDirXArray[xyCoordCenter]*(double)(-xx)+votingDirYArray[xyCoordCenter]*(double)(-xx))/dis ));
	//				if( (teta<=delta/2) )
	//				{
	//					xCoordVot = xCoord-xx;
	//					yCoordVot = yCoord2-xx;
	//					if( (0<=xCoordVot) && (xCoordVot<votingMagSize[0]) && (0<=yCoordVot) && (yCoordVot<votingMagSize[1]) )
	//					{
	//							xyCoord = (yCoordVot)*votingMagSize[0]+xCoordVot;
	//							if( votingImaArray[xyCoord] > maximImaVot )
	//							{ 
	//								maximImaVot = votingImaArray[xyCoord];
	//								xCoordmaximImaVot = -xx;
	//								yCoordmaximImaVot = -xx;
	//							}
	//					}
	//				}
	//			}
	//			temp = sqrt( double(xCoordmaximImaVot*xCoordmaximImaVot + yCoordmaximImaVot*yCoordmaximImaVot));
	//			votingDirXArray[xyCoordCenter] = ((double)xCoordmaximImaVot)/temp;
	//			votingDirYArray[xyCoordCenter] = ((double)yCoordmaximImaVot)/temp;
	//		}
	//
	//
	//		delta = delta -deltaDecrea*3.1416/180;
	//	//// --------------------------------------------------------------------------------------------------------------------------------------------------
	//
	//	// Store Voting Image
	//	if( writeImage< VotingImaType, OutputImageType >(VotingImaImage, VotingImaName )){
	//		return 1;
	//	}
	//	clock_t partialTimeEnd=clock();
	//	cout << endl << "	Time this Iteration: " << double(diffclock(partialTimeEnd,partialTimeBegin)) << " s"<< endl;
	//
	//	
	//
	//
	//
	//
	//	} // Fin delas iteraciones
	//
	//
	//
	//
	//
	//	//// --------------------------------------------------------------------------------------------------------------------------------------------------
	//
	//	// Store Voting Image
	//	if( writeImage< VotingImaType, OutputImageType >(VotingImaImage, VotingImaName )){
	//		return 1;
	//	}
	//
	//	// --------------------------------------------------------------------------------------------------------------------------------------------------
	//
	//	//Store Derivative in X
	//	if( writeImage< VotingDirType, OutputImageType >(votingDirX, votingDirXName )){
	//		return 1;
	//	}
	//
	//	// --------------------------------------------------------------------------------------------------------------------------------------------------
	//
	//	//// Store Derivative in Y
	//	if( writeImage< VotingDirType, OutputImageType >(votingDirY, votingDirYName )){
	//		return 1;
	//	}
	//
	//	// --------------------------------------------------------------------------------------------------------------------------------------------------
	//
	//	// Store Magnitude fo the Derivative
	//	if( writeImage< VotingMagType, OutputImageType >(VotingMagImage, VotingMagName )){
	//		return 1;
	//	}
	//
	//	//// --------------------------------------------------------------------------------------------------------------------------------------------------
	//
	//	// Store Magnitude fo the Derivative
	//	if( writeImage< VotingPixType, OutputImageType >(VotingPixImage, votingPixName )){
	//		return 1;
	//	}




	clock_t end=clock();
	cout << "Time elapsed: " << double(nftkVotingGlobal::diffclock(end,begin)) << " s";
	return 1;
};

// ############################################################################################################################################################################
// ############################################################################################################################################################################