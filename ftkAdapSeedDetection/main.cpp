// ############################################################################################################################################################################
// ############################################################################################################################################################################

// Main functions
//#include"ftkVoting.h"
#include"ftkVoting_3D.h"
#include"ftkVotingGlobal.h"

// For Extracting slide by slice image of the 3d Image
#include "itkImage.h"
#include "itkExtractImageFilter.h"

//using namespace nftkVot;



// ############################################################################################################################################################################
// ############################################################################################################################################################################

int main( int argc, char * argv[] ){

	clock_t begin=clock();

	// User Parameters
	char* inputImageName = argv[1];
	char* folder = argv[2];
	//int hmax = atoi(argv[2]);

	// Initial Parameters
	int hmin = 1;
	int hmax = 25;//35;//77; //15 was ok, 25 just to test
	int radius = 0.7*hmax;//25;//10;//77;
	double min_grad = 0.1;
	//threshold (for picking seed points, not necesary for now)
	double scale = 1.5; // Scale for computing the gradient using Gradient
	//zc_only TRUE - voting points should also be the zero-crossing points of the image; FALSE - otherwise


	//// ###################################### 2D
	//// Set Up the Reader 2D
	//nftkVot::InputImageType::Pointer inputImage = nftkVot::readImage< nftkVot::InputImageType, nftkVot::InputImageType >( inputImageName );

	//ftkVoting voteMain;
	//voteMain.setParams(hmin,hmax,radius,min_grad,scale);
	//voteMain.setPrefix("temp/");
	//voteMain.compute(inputImage); // DUDA COMO HACER PARA ENVIAR UN CONST POINTER, DESPUES DE QUE HE LEIDO LA IMAGEN COMO POINTER ??


	// ###################################### 3D
		// Set Up the Reader 3D
	
		// Input Image Type
	
		nftkVot::InputImageType_3D::Pointer inputImage = nftkVot::readImage_3D< nftkVot::InputImageType_3D_16, nftkVot::InputImageType_3D >( inputImageName );
	
		ftkVoting_3D voteMain_3D;
		voteMain_3D.setParams(hmin,hmax,radius,min_grad,scale);
		voteMain_3D.setPrefix("temp/");
		voteMain_3D.compute(inputImage);

	//// ###################################### 2D Slice by Slice

	//// Set Up the Reader 2D
	//nftkVot::InputImageType_3D::Pointer inputImage_3D = nftkVot::readImage_3D< nftkVot::InputImageType_3D_16, nftkVot::InputImageType_3D >( inputImageName );
	//inputImage_3D->Update();


	//double dxx = -0.1;

	//double dyy =0.3;
	//cout<<endl<<" ANGULO: "<<atan(dyy/dxx)<<" "<<atan(dxx/dyy);;

	//typedef itk::ExtractImageFilter< nftkVot::InputImageType_3D, nftkVot::InputImageType > FilterType;
	//FilterType::Pointer filter = FilterType::New();
	//filter->SetDirectionCollapseToSubmatrix();

	//filter->SetInput( inputImage_3D );


	//nftkVot::InputImageType_3D::RegionType inputRegion = inputImage_3D->GetLargestPossibleRegion();
	//nftkVot::InputImageType_3D::SizeType size = inputRegion.GetSize();
	//size[2] = 0;

	//nftkVot::InputImageType_3D::IndexType start = inputRegion.GetIndex();
	//start[2] = 1;

	//nftkVot::InputImageType_3D::RegionType desiredRegion;
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


	////nftkVot::InputImageType_3D::Pointer votingBySlice = nftkVot::InputImageType_3D::New(); //VotingDirType = double
	////votingBySlice->SetRegions( inputImage_3D->GetRequestedRegion()); // IMPORTANTE PARA CREAR UNA IMAGEN NUEVA EN BASE A UNA QUE YA EXISTE EN VEZ DE PONERME A LEER LOS TAMANOS Y LAS REGIONES
	////votingBySlice->Allocate();

	////nftkVot::InputImageType_3D::Pointer votingBySlice = nftkVot::readImage_3D< nftkVot::InputImageType_3D_16, nftkVot::InputImageType_3D >( inputImageName );
	////votingBySlice->Update();
	////int nxxX = votingBySlice->GetLargestPossibleRegion().GetSize()[0];
	////int nyxX = votingBySlice->GetLargestPossibleRegion().GetSize()[1];
	////int nzxX = votingBySlice->GetLargestPossibleRegion().GetSize()[2];
	////int npixX = nxxX*nyxX*nzxX;
	////nftkVot::InputImageType_3D::PixelType * votingBySliceArray = votingBySlice->GetBufferPointer();
	////for( unsigned int tr=0; tr<npixX; ++tr ){
	////	votingBySliceArray[0] = 0;
	////}

	//nftkVot::InputImageType_3D::Pointer votingBySlice = nftkVot::InputImageType_3D::New();
	//nftkVot::InputImageType_3D::IndexType start22;
	//start22[0] = 0;
	//start22[1] = 0;
	//start22[2] = 0;
	//nftkVot::InputImageType_3D::SizeType size22;
	//double bw = sqrt((double)(radius*radius + hmax*hmax)) + 3; // Que es esto ??
	//const int bw2 = 2*bw;
	//size22[0] = inputImage_3D->GetLargestPossibleRegion().GetSize()[0]+bw2;
	//size22[1] = inputImage_3D->GetLargestPossibleRegion().GetSize()[1]+bw2;
	//size22[2] = inputImage_3D->GetLargestPossibleRegion().GetSize()[2];
	//nftkVot::InputImageType_3D::RegionType region;
	//region.SetSize( size22 );
	//region.SetIndex( start22 );
	//votingBySlice->SetRegions( region );
	//votingBySlice->Allocate();
	//const nftkVot::InputImageType_3D::PixelType ceros = 0;
	//votingBySlice->FillBuffer( ceros );
	//votingBySlice->Update();

	//nftkVot::InputImageType_3D::Pointer votingBySliceLastbigSapan = nftkVot::InputImageType_3D::New();
	//votingBySliceLastbigSapan->SetRegions( region );
	//votingBySliceLastbigSapan->Allocate();
	//votingBySliceLastbigSapan->FillBuffer( ceros );
	//votingBySliceLastbigSapan->Update();

	//nftkVot::InputImageType_3D::Pointer votingBySliceLastbigSapanProb = nftkVot::InputImageType_3D::New();
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






	clock_t end=clock();
	std::cout << "Time elapsed: " << double(nftkVot::diffclock(end,begin)) << " s";

	std::cout << std::endl;
	return 0;
};

// ############################################################################################################################################################################
// ############################################################################################################################################################################