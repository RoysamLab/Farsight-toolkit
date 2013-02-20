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
const double ftkVoting_3D::epsilon((double)0.0001);
const double ftkVoting_3D::pi = (double)3.1415927;

// ############################################################################################################################################################################
ftkVoting_3D::ftkVoting_3D()
{
	// Initial Parameters
	//reg_type (have not underestand what is this for)
	_hmin = 1;
	_hmax = 20;
	_radius = 4;
	_min_grad = 0.005;
	//_threshold (for picking seed points, not necesary for now)
	_scale = 1.5; // Scale for computing the gradient using DoG

	// Default Parameters for quantizing the direction of voting
	_NN_dir = 500;
	_sigmaG = 1.5;
	_Z_factor = 2;
	_maxSpan = 10;
	_intSpan = _maxSpan; // Start with the maximum span
}

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
void ftkVoting_3D::setParams(	int hmin, int hmax,	int radius,	double min_grad, double scale )
{
	_hmin = hmin;
	_hmax = hmax;
	_radius = radius;
	_min_grad = min_grad;
	_scale = scale;
}

// ############################################################################################################################################################################
void ftkVoting_3D::compute(nftkVot::InputImageType_3D::Pointer inputImage)
{
	_inputImage = inputImage;

	nx = inputImage->GetLargestPossibleRegion().GetSize()[0];
	ny = inputImage->GetLargestPossibleRegion().GetSize()[1];
	nz = inputImage->GetLargestPossibleRegion().GetSize()[2];
	npix = nx*ny*nz;
	_npix_org = npix;

	_nx_org = nx;
	_ny_org = ny;
	_nz_org = nz;

	std::cout << std::endl << std::endl << "--> COMPUTING DERIVATIVES <--";

	// Derivatives using Gaussian
	typedef itk::RecursiveGaussianImageFilter<nftkVot::InputImageType_3D,nftkVot::InputImageType_3D > ShortFilterType;
	typedef itk::ImageDuplicator< nftkVot::InputImageType_3D > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();

	// Calculate the derivative in X
	ShortFilterType::Pointer dx = ShortFilterType::New();
	dx->SetDirection( 0 ); //dx works along x
	dx->SetSigma( _sigmaG );
	dx->SetFirstOrder();
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
	nftkVot::InputImageType_3D::Pointer votingDirX_3D = duplicator->GetOutput(); //first order derivative along x
	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";


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
	nftkVot::InputImageType_3D::Pointer votingDirY_3D = duplicator->GetOutput();
	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";


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
	nftkVot::InputImageType_3D::Pointer votingDirZ_3D = duplicator->GetOutput();
	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";


	// Pointers for the X,Y and Z directions
	nftkVot::InputImageType_3D::PixelType * votingDirX_3DArray = votingDirX_3D->GetBufferPointer();
	nftkVot::InputImageType_3D::PixelType * votingDirY_3DArray = votingDirY_3D->GetBufferPointer();
	nftkVot::InputImageType_3D::PixelType * votingDirZ_3DArray = votingDirZ_3D->GetBufferPointer();

	// Gradient magnitud
	// Duplicate image X (to create a copy, the gradient magnitude)
	typedef itk::ImageDuplicator< nftkVot::InputImageType_3D > DuplicatorType_2;
	DuplicatorType_2::Pointer duplicator_2 = DuplicatorType_2::New();
	duplicator_2->SetInputImage( votingDirX_3D );
	try{
		duplicator_2->Update();}
	catch( itk::ExceptionObject & err ){
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
	}
	nftkVot::InputImageType_3D::Pointer votingDirXYZ_3D = duplicator_2->GetOutput(); //first order derivative along z

	// Pointer fot the XYZ image
	nftkVot::InputImageType_3D::PixelType * votingDirXYZ_3DArray = votingDirXYZ_3D->GetBufferPointer();

	// Gradient magnitude and normalize gradient X, Y and Z
#pragma omp parallel for
	for( int tt=0; tt<npix; ++tt ){
		votingDirXYZ_3DArray[tt] = sqrt(votingDirX_3DArray[tt]*votingDirX_3DArray[tt]+votingDirY_3DArray[tt]*votingDirY_3DArray[tt]+votingDirZ_3DArray[tt]*votingDirZ_3DArray[tt]*_Z_factor);
		if(votingDirXYZ_3DArray[tt]>epsilon){
			votingDirX_3DArray[tt] = votingDirX_3DArray[tt]/votingDirXYZ_3DArray[tt];
			votingDirY_3DArray[tt] = votingDirY_3DArray[tt]/votingDirXYZ_3DArray[tt];
			votingDirZ_3DArray[tt] = votingDirZ_3DArray[tt]*_Z_factor/votingDirXYZ_3DArray[tt];
		}
	}

	// Calculate the Maximum value of the gradient
	typedef itk::MinimumMaximumImageCalculator < nftkVot::InputImageType_3D > ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
	imageCalculatorFilter->SetImage( votingDirXYZ_3D);
	imageCalculatorFilter->Compute();
	nftkVot::InputImageType_3D::PixelType maxVotVal = imageCalculatorFilter->GetMaximum();

	// Scale the values of the gradient by the maximum value of the gradiente, to ensure that the maximum value is 1
	std::cout << std::endl << "Scalling the values of the gradiet (0 to 1) -> ";
	_t1_begin=clock();
#pragma omp parallel for
	for(int i=0; i<npix; i++) {
		if (votingDirXYZ_3DArray[i]<epsilon)
		{
			votingDirXYZ_3DArray[i] = 0;
		}
		else
		{
			votingDirXYZ_3DArray[i] = votingDirXYZ_3DArray[i]/maxVotVal;
		}
	}
	std::cout << "Done.";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";


	//// Store the derivative
	//std::cout << std::endl << "Store the derivative";
	//_t1_begin=clock();
	//
	//	// Save the derivatives
	//	string filenameDerivativeX = "output/out_ftkDerivativeX.tif";
	//	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_16 >(votingDirX_3D, filenameDerivativeX.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//	string filenameDerivativeY = "output/out_ftkDerivativeY.tif";
	//	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_16 >(votingDirY_3D, filenameDerivativeY.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//	string filenameDerivativeZ = "output/out_ftkDerivativeZ.tif";
	//	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_16 >(votingDirZ_3D, filenameDerivativeZ.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//	string filenameDerivativeXYZ = "output/out_ftkDerivativeXYZ.tif";
	//	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_16 >(votingDirXYZ_3D, filenameDerivativeXYZ.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//std::cout << " -> Done.";
	//_t1_end=clock();
	//std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";


	std::cout << std::endl << "Now we are going to compute the cones";
	computeCones_3D(_hmin, _hmax, _radius);
	std::cout << " -> Done.";
	//computeCones_3D_prob(_hmin, _hmax, _radius); // Just for one moment

	// Padding
	bw = sqrt((double)(_radius*_radius + _hmax*_hmax)) + 3; // Que es esto ??
	const int bw2 = 2*bw;
	_bw = bw;

	// Calculate the voting poitns
	std::cout << std::endl << "Calculate Voting Points 3D -> ";
	_t1_begin=clock();
// This can be parallel
	VPoint3D vp_3D;
	for(int z = 0; z < nz; ++z)
	{
		for(int y = 0; y < ny; ++y)
		{
			for(int x = 0; x < nx; ++x)
			{
				int k = x+nx*y+nx*ny*z;
				if (votingDirXYZ_3DArray[k] > _min_grad) // This is the original one, only vote if the mag of gradient is greater than epsilon (and gra_min) previously thresholded
				{
					vp_3D.x = x + bw;
					vp_3D.y = y + bw;
					vp_3D.z = z + bw;
					vp_3D.mag = votingDirXYZ_3DArray[k]; 
					vp_3D.posOri = k; // This is the original position in the original image
					vp_3D.scale = _hmax;
					_voting_points_3D.push_back(vp_3D);
				}
			}
		}
	}
	std::cout << "Done.";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";

	// Calculate the direction of the voting poitns
	std::cout << std::endl << "Quantize Voting Points Directions in 3D -> ";
	_t1_begin=clock();
	int limit = _voting_points_3D.size();
#pragma omp parallel for
	for( int i = 0; i<limit; ++i )\
	{
		int k = _voting_points_3D[i].posOri;
		_voting_points_3D[i].direc_vote = computeDirecIndex_3D(votingDirX_3DArray[k], votingDirY_3DArray[k], votingDirZ_3DArray[k]); 
	}
	std::cout << "Done.";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";


	// Add the padding to the image size
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

}

//// ############################################################################################################################################################################
void ftkVoting_3D::computeCones_3D(int hmin, int hmax, int radius)
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

	//int z1_new;
	//int R_new;
	//double R_dou_new;
	//int R_quan_new;
	//double x_nor_new, y_nor_new, z_nor_new;

	//ftkWPoint3D wp_new;

#pragma omp parallel for
	for( int xx=-hmax; xx<=hmax; ++xx )
	{

	int z1_new;
	int R_new;
	double R_dou_new;
	int R_quan_new;
	double x_nor_new, y_nor_new, z_nor_new;

	ftkWPoint3D wp_new;

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
#pragma omp critical
					{
					_conesPru_3D_new[maxmax_pos][R_quan_new].push_back(wp_new);
					}
				}
			}
		}
	}

	_voteDirec_3D_new = std::vector< std::vector< std::vector< int > > >(_NN_dir,std::vector< std::vector <int> >(_maxSpan));
	double max_dotpro_delta = atan((double)_radius/(double)_hmax)/_maxSpan; // Max span

	//INEFFICIENT, 
	//#pragma omp parallel for 
	for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	{
		for( int uu2=0; uu2<_NN_dir; ++uu2 ) // All possible directions
		{
			double dotproduct = acos(_conesPru_3D_new[uu].dxx*_conesPru_3D_new[uu2].dxx+_conesPru_3D_new[uu].dyy*_conesPru_3D_new[uu2].dyy+_conesPru_3D_new[uu].dzz*_conesPru_3D_new[uu2].dzz);

			int angle = ceil(dotproduct/max_dotpro_delta);
			if( (0<=angle) && (angle<_maxSpan) )
			{
				_voteDirec_3D_new[uu][angle].push_back(uu2);
				//!!!
				//std::cout<<std::endl<<angle;
			}
		}
	}

	std::cout << " -> Done";
	clock_t end_1=clock();
	std::cout << "Time elapsed: " << double(nftkVot::diffclock(end_1,begin_1)) << " s";


	// NOT QUANTIZED

	std::cout << std::endl << "Computing the cones without qua";
	clock_t begin_3=clock();

	_conesPru_3D_new_notqu = std::vector < ftkCone3D_new >(_NN_dir); // Todos los conos posibles

	double dlong_para_notqu = pi*(3-sqrt(5.0));
	double dz_para_notqu = 2.0/_NN_dir;
	double long__para_notqu = 0;
	double z_para_notqu = 1-dz_para_notqu/2;

	for( int uu=0; uu<_NN_dir; uu++ ) // Phi
	{
		double r_para_notqu = sqrt(1-z_para_notqu*z_para_notqu);
		_conesPru_3D_new_notqu[uu].dxx = cos(long__para_notqu)*r_para_notqu;
		_conesPru_3D_new_notqu[uu].dyy = sin(long__para_notqu)*r_para_notqu;
		_conesPru_3D_new_notqu[uu].dzz = z_para_notqu;
		z_para_notqu    = z_para_notqu - dz_para_notqu;
		long__para_notqu = long__para_notqu + dlong_para_notqu;
	}

#pragma omp parallel for
	for( int xx=-hmax; xx<=hmax; ++xx )
	{

	int z1_new;
	int R_new;
	double R_dou_new;
	int R_quan_new;
	double x_nor_new, y_nor_new, z_nor_new;

	ftkWPoint3D wp_new;

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

					ftkWPoint3D wp_new;
					wp_new.x = xx;
					wp_new.y = yy;
					wp_new.z = zz;
					wp_new.w = 1; // In case of deciding to put some weight

					double maxmax_dot = 0;
					int maxmax_pos = 0;

					for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
					{
						double dotprod = _conesPru_3D_new_notqu[uu].dxx*x_nor_new+_conesPru_3D_new_notqu[uu].dyy*y_nor_new+_conesPru_3D_new_notqu[uu].dzz*z_nor_new;
						if ( maxmax_dot<dotprod)
						{
							maxmax_dot = dotprod;
							maxmax_pos = uu;
						}
					}
#pragma omp critical
					{
					_conesPru_3D_new_notqu[maxmax_pos].push_back(wp_new);
					}
				}
			}
		}
	}

	_voteDirec_3D_new_notqu = std::vector< std::vector< std::vector< int > > >(_NN_dir,std::vector< std::vector <int> >(_maxSpan));
	double max_dotpro_delta_para_notqu = atan((double)_radius/(double)_hmax)/_maxSpan; // Max span

	//INEFFICIENT, 
	//#pragma omp parallel for 
	for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	{
		for( int uu2=0; uu2<_NN_dir; ++uu2 ) // All possible directions
		{
			double dotproduct = acos(_conesPru_3D_new_notqu[uu].dxx*_conesPru_3D_new_notqu[uu2].dxx+_conesPru_3D_new_notqu[uu].dyy*_conesPru_3D_new_notqu[uu2].dyy+_conesPru_3D_new_notqu[uu].dzz*_conesPru_3D_new_notqu[uu2].dzz);

			int angle = ceil(dotproduct/max_dotpro_delta_para_notqu);
			if( (0<=angle) && (angle<_maxSpan) )
			{
				_voteDirec_3D_new_notqu[uu][angle].push_back(uu2);
			}
		}
	}
	std::cout << " -> Done. ";
	clock_t end_3=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(end_3,begin_3)) << " s";
}



//// ############################################################################################################################################################################
void ftkVoting_3D::vote()
{

	std::cout << std::endl << std::endl << "--> VOTING STARTED <--";

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
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";

	VotingDirType_3D::PixelType * votingSumArray = _votingSumVotes->GetBufferPointer();
	VotingDirType_3D::PixelType * votingMaskArray = _votingMaskVotes->GetBufferPointer();
	VotingDirType_3D::PixelType * imageOfVotingPixelsArray = imageOfVotingPixels->GetBufferPointer();
	VotingDirType_3D::PixelType * imageOfConexPixelsArray = imageOfConexPixels->GetBufferPointer();
	VotingDirType_3D::PixelType * imageGradXPixelsArray = imageGradXPixels->GetBufferPointer();
	VotingDirType_3D::PixelType * imageGradYPixelsArray = imageGradYPixels->GetBufferPointer();
	VotingDirType_3D::PixelType * imageMaxPixelsArray = imageMaxPixels->GetBufferPointer();


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
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";


	// Update the voting cones to the new size (paddng)
	// FALTA REVIAR ESTE PEDASO
#pragma omp parallel for
	for( int uu=0; uu<_NN_dir; uu++ )
	{
		_conesPru_3D_new[uu].setOffset();
	}

	// Update the voting cones to the new size (paddng)
	// FALTA REVIAR ESTE PEDASO
#pragma omp parallel for
	for( int uu=0; uu<_NN_dir; uu++ )
	{
		_conesPru_3D_new_notqu[uu].setOffset();
	}


	// Store the image of voting points
	std::cout << std::endl << "Store the image of voting points .";
	_t1_begin=clock();
	int temp34 = _voting_points_3D.size();
#pragma omp parallel for
	for( int i=0; i<temp34; ++i )
	{
		imageOfVotingPixelsArray[_voting_points_3D.at(i).pos] = 1;
	}
	std::string filenameCones_votpixe = "output/Vot_Pixels.tif";
	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(imageOfVotingPixels, filenameCones_votpixe.c_str() )){
		std::cout<<std::endl<<"\tProblema escribiendo";
	}
	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";


	// Store the original image with padding
	std::cout << std::endl << "Store the original image with padding (only necessary to show results, disable otherwise).";
	_t1_begin=clock();
	nftkVot::InputImageType_3D::PixelType * inputImageArray = _inputImage->GetBufferPointer();
#pragma omp parallel for
	for( int xx = 0; xx<_nx_org; ++xx )
	{
		for( int yy = 0; yy<_ny_org; ++yy )
		{
			for( int zz = 0; zz<_nz_org; ++zz )
			{
				votingSumArray[bw+bw*(2*bw+_nx_org)+bw*(2*bw+_nx_org)*(2*bw+_ny_org)+ xx+yy*(_nx_org+2*bw)+zz*(_nx_org+2*bw)*(_ny_org+2*bw)] = inputImageArray[xx+yy*_nx_org+zz*_nx_org*_ny_org];
			}
		}
	}
	std::string filenameOrig = "output/Vot_OriginalImage.tif";
	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameOrig.c_str() )){
		std::cout<<std::endl<<"\tProblema escribiendo";
	}
	memset(votingSumArray, 0,  npix*sizeof(VotingDirType_3D::PixelType)); // Errase the 
	std::cout << " -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";



	//// One Vote just to test
	//std::cout << std::endl << "One Vote just to test.";
	//_t1_begin=clock();

	//int temp2 = _voting_points_3D.size();
	//for( int i=0; i<temp2; ++i )
	//{
	//	std::cout<<std::endl<<"Vote: "<<temp2-i;
	//	int dir_vote = _voting_points_3D.at(i).direc_vote;
	//	for( int angle_int = 0; angle_int<_maxSpan; ++angle_int )
	//	{
	//		//std::cout<<std::endl<<"SIZE: "<<_voteDirec_3D_new[uu][angle_int].size();
	//		//// instead use vector iterator
	//		for( int vv = 0; vv<_voteDirec_3D_new_para[dir_vote][angle_int].size(); ++vv )
	//		{
	//			//std::cout<<std::endl<<"VV: "<<vv;
	//			int temp_vv = _voteDirec_3D_new_para[dir_vote][angle_int][vv];
	//			for( int raddd=0; raddd<_conesPru_3D_new_para[temp_vv].size(); ++raddd )
	//			{
	//				//std::cout<<std::endl<<"INTER: "<<_conesPru_3D_new[temp_vv][raddd].size();
	//				for( int bin_cont=0; bin_cont < _conesPru_3D_new_para[temp_vv][raddd].size(); ++bin_cont )
	//				{
	//					int x_posi = _conesPru_3D_new_para[temp_vv][raddd][bin_cont].x;
	//					int y_posi = _conesPru_3D_new_para[temp_vv][raddd][bin_cont].y;
	//					int z_posi = _conesPru_3D_new_para[temp_vv][raddd][bin_cont].z;
	//					votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi+_voting_points_3D.at(i).pos] = votingSumArray[x_posi+nx*y_posi+nx*ny*z_posi+_voting_points_3D.at(i).pos] + _voting_points_3D.at(i).mag;
	//				}
	//			}
	//		}
	//	}
	//}
	//std::string filenameCones_2as = "output/oneVote.tif";
	//if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//}
	//std::cout << " -> Done. ";
	//_t1_end=clock();
	//std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";
	//nftkVot::stopProgram();


	// Let's Vote
	std::cout << std::endl << "Let's Vote.";
	_t1_begin=clock();
for( int spanofvote = 10; spanofvote>0; --spanofvote )
{
	std::cout << std::endl << "\tVote (version 3)" << _intSpan - spanofvote + 1;
	_t2_begin=clock();
	int temp2 = _voting_points_3D.size();
#pragma omp parallel for
	for( int i=0; i<temp2; ++i ) // All the voting points
	{
		//// VERSION 1
		//VotingDirType_3D::PixelType * votingMaskArray_pos = _votingMaskVotes->GetBufferPointer()+_voting_points_3D.at(i).pos;

		//int dir_vote = _voting_points_3D.at(i).direc_vote;
		////#pragma omp parallel for
		//for( int angle_int = 0; angle_int<spanofvote; ++angle_int ) // The different angles (1-10)
		//{
		//	for( int vv = 0; vv<_voteDirec_3D_new_notqu[dir_vote][angle_int].size(); ++vv )	// All the bins at a given dot product distance from the given direction (dir_vote)
		//	{
		//		int temp_vv = _voteDirec_3D_new_notqu[dir_vote][angle_int][vv];
		//		for( int bin_cont=0; bin_cont < _conesPru_3D_new_notqu[temp_vv].size(); ++bin_cont )
		//		{
		//			votingMaskArray_pos[_conesPru_3D_new_notqu[temp_vv][bin_cont].off] = votingMaskArray_pos[_conesPru_3D_new_notqu[temp_vv][bin_cont].off] + _voting_points_3D.at(i).mag;
		//		}
		//	}
		//}


// 		// VERSION 2
// 		VotingDirType_3D::PixelType * votingMaskArray_pos = _votingMaskVotes->GetBufferPointer()+_voting_points_3D.at(i).pos;
// 
// 		int dir_vote = _voting_points_3D.at(i).direc_vote;
// 		//#pragma omp parallel for
// 		for( int angle_int = 0; angle_int<spanofvote; ++angle_int ) // The different angles (1-10)
// 		{
// 			std::vector< int >::iterator it_dir_3D;
// 			for( it_dir_3D = _voteDirec_3D_new_notqu[dir_vote][angle_int].begin(); it_dir_3D != _voteDirec_3D_new_notqu[dir_vote][angle_int].end(); ++it_dir_3D )
// 			{
// 				ftkCone3D_new::iterator it_cone_3D;
// 				for( it_cone_3D = _conesPru_3D_new_notqu[(*it_dir_3D)].begin(); it_cone_3D != _conesPru_3D_new_notqu[(*it_dir_3D)].end(); ++it_cone_3D )
// 				{
// 					votingMaskArray_pos[(*it_cone_3D).off] = votingMaskArray_pos[(*it_cone_3D).off] + _voting_points_3D.at(i).mag;
// 
// 				}
// 			}
// 		}
		
// 		// VERSION 3
// 		VotingDirType_3D::PixelType * votingMaskArray_pos = _votingMaskVotes->GetBufferPointer()+_voting_points_3D.at(i).pos;
// 
// 		int dir_vote = _voting_points_3D.at(i).direc_vote;
// 		//#pragma omp parallel for
// 		for( int angle_int = 0; angle_int<spanofvote; ++angle_int ) // The different angles (1-10)
// 		{
// 			std::vector< int >::iterator it_dir_3D;
// 			for( it_dir_3D = _voteDirec_3D_new[dir_vote][angle_int].begin(); it_dir_3D != _voteDirec_3D_new[dir_vote][angle_int].end(); ++it_dir_3D )
// 			{
// 				ftkCone3D::iterator it_cone_3D;
// 				ftkBins3D::iterator it_bin_3D;
// 				for( it_cone_3D = _conesPru_3D_new[(*it_dir_3D)].begin(); it_cone_3D != _conesPru_3D_new[(*it_dir_3D)].end(); ++it_cone_3D )
// 				{
// 					for( it_bin_3D = (*it_cone_3D).begin(); it_bin_3D != (*it_cone_3D).end(); ++it_bin_3D )
// 					{
// 						votingMaskArray_pos[(*it_bin_3D).off] = votingMaskArray_pos[(*it_bin_3D).off] + _voting_points_3D.at(i).mag;
// 					}
// 
// 				}
// 			}
// 		}

		// VERSION 1
		VotingDirType_3D::PixelType * votingMaskArray_pos = _votingMaskVotes->GetBufferPointer()+_voting_points_3D.at(i).pos;
		
		ftkBins3D * poin_bin;

		int dir_vote = _voting_points_3D.at(i).direc_vote;
		//#pragma omp parallel for
		for( int angle_int = 0; angle_int<spanofvote; ++angle_int ) // The different angles (1-10)
		{
			for( int vv = 0; vv<_voteDirec_3D_new[dir_vote][angle_int].size(); ++vv )	// All the bins at a given dot product distance from the given direction (dir_vote)
			{
				int temp_vv = _voteDirec_3D_new[dir_vote][angle_int][vv];
				for( int bin_cont=0; bin_cont < _voting_points_3D.at(i).scale/*_conesPru_3D_new[temp_vv].size()*/; ++bin_cont )
				{
					poin_bin =  &_conesPru_3D_new[temp_vv][bin_cont];
					
					for( int votes=0; votes < (*poin_bin).size(); ++votes )
						votingMaskArray_pos[(*poin_bin)[votes].off] = votingMaskArray_pos[(*poin_bin)[votes].off] + _voting_points_3D.at(i).mag;
				}
			}
		}

	}
	memcpy(votingSumArray, votingMaskArray, npix*sizeof(VotingDirType_3D::PixelType)); // !!! This should be sum = sum + mask
	memset(votingMaskArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	std::cout << " -> Done. ";
	_t2_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t2_end,_t2_begin)) << " s";



	// Update the direction of the votes
	// !!! ??? OJO QUE AHOR EL UPDATE SE ESTA HACIENDO CON LA SUMA, NO CON LA MASCARA
	std::cout << std::endl << "\tLet's Update Directions.";
	_t3_begin=clock();
	int temp3 = _voting_points_3D.size();
#pragma omp parallel for
	for( int i=0; i<temp3; ++i ) // All the voting points
	{
		updateDirection_3D(_voting_points_3D.at(i));
	}
	std::cout << " -> Done. ";
	_t3_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t3_end,_t3_begin)) << " s";

	// Writing results
	std::cout << std::endl << "\tWriting results.";
	_t4_begin=clock();
	std::stringstream outas;
	outas << spanofvote;
	std::string sas = outas.str();
	std::string filenameCones_2as = "output/"+sas+"_Vote_sum.tif";
	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
		std::cout<<std::endl<<"\tProblema escribiendo";
	}
	std::string filenameCones_2as_mask = "output/"+sas+"_Vote_mask.tif";
	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingMaskVotes, filenameCones_2as_mask.c_str() )){
		std::cout<<std::endl<<"\tProblema escribiendo";
	}
	std::cout << " -> Done. ";
	_t4_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t4_end,_t4_begin)) << " s";
	//nftkVot::stopProgram();

}

	std::cout << std::endl << "\t -> Done. ";
	_t1_end=clock();
	std::cout << "Time: " << double(nftkVot::diffclock(_t1_end,_t1_begin)) << " s";
	nftkVot::stopProgram();
















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
	std::cout << "Time elapsed First: " << double(nftkVot::diffclock(end33,begin33)) << " s";
	nftkVot::stopProgram();

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
	std::cout << "Time elapsed First: " << double(nftkVot::diffclock(end37,begin37)) << " s";
	nftkVot::stopProgram();

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
	std::cout << "Time elapsed First: " << double(nftkVot::diffclock(end39,begin39)) << " s";
	nftkVot::stopProgram();

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
	std::cout << "Time elapsed First: " << double(nftkVot::diffclock(end38,begin38)) << " s";
	nftkVot::stopProgram();

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
	std::cout << "Time elapsed First: " << double(nftkVot::diffclock(end34,begin34)) << " s";
	nftkVot::stopProgram();



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
	//	string filenameCones_2as = "output/out_new_votingSumVotes_cone"+sas+".tif";
	//	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
	//		std::cout<<std::endl<<"\tProblema escribiendo";
	//	}
	//}



	// PARALLEL

	clock_t begin=clock();
	for( int oo=0; oo<1; ++oo)
	{
		//#pragma omp parallel for
		for( int uu=0; uu<200/*_NN_dir*/; ++uu )
		{
			//std::cout<<std::endl<<"UU_new: "<<oo<<" "<<uu;
			for( int angle_int = 0; angle_int<_maxSpan; ++angle_int )
			{
				//std::cout<<std::endl<<"SIZE: "<<_voteDirec_3D_new[uu][angle_int].size();
				//// instead use vector iterator
				//				#pragma omp parallel for private(raddd,bin_cont)
				for( int vv = 0; vv<_voteDirec_3D_new[uu][angle_int].size(); ++vv )
				{
					//std::cout<<std::endl<<"VV: "<<vv;
					int temp_vv = _voteDirec_3D_new[uu][angle_int][vv];
					for( int raddd=0; raddd<60; ++raddd )
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
				//string filenameCones_2as = "output/out_new_votingSumVotes_cone"+sas+".tif";
				//if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
				//	std::cout<<std::endl<<"\tProblema escribiendo";
				//}
			}
		}
	}
	clock_t end=clock();
	std::cout << "Time elapsed Second Parallel: " << double(nftkVot::diffclock(end,begin)) << " s";
	//nftkVot::stopProgram();


	// NOT PARALLEL

	clock_t begin67=clock();
	for( int oo=0; oo<1; ++oo)
	{
		//#pragma omp parallel for
		for( int uu=0; uu<200/*_NN_dir*/; ++uu )
		{
			//std::cout<<std::endl<<"UU_new: "<<oo<<" "<<uu;
			for( int angle_int = 0; angle_int<_maxSpan; ++angle_int )
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
				//string filenameCones_2as = "output/out_new_votingSumVotes_cone"+sas+".tif";
				//if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
				//	std::cout<<std::endl<<"\tProblema escribiendo";
				//}
			}
		}
	}
	clock_t end67=clock();
	std::cout << "Time elapsed Second Parallel: " << double(nftkVot::diffclock(end67,begin67)) << " s";
	nftkVot::stopProgram();








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

			//#pragma omp parallel for private(vv,raddd,bin_cont)
			for( int angle_int = 0; angle_int<_maxSpan; ++angle_int )
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
				//string filenameCones_2as = "output/out_new_votingSumVotes_cone"+sas+".tif";
				//if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
				//	std::cout<<std::endl<<"\tProblema escribiendo";
				//}
			}
		}
	}

	clock_t end3=clock();
	std::cout << "Time elapsed Third: " << double(nftkVot::diffclock(end3,begin3)) << " s";

	nftkVot::stopProgram();




	clock_t begin2=clock();


	for( int oo=0; oo<100; ++oo)
	{
		//#pragma omp parallel for
		for( int uu=0; uu<_NN_dir; ++uu )
		{
			//std::cout<<std::endl<<"UU: "<<uu;
			for( int angle_int = 0; angle_int<_maxSpan; ++angle_int )
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
				//string filenameCones_2as = "output/out_new_votingSumVotes_cone"+sas+".tif";
				//if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
				//	std::cout<<std::endl<<"\tProblema escribiendo";
				//}
			}
		}
	}

	clock_t end2=clock();
	std::cout << "Time elapsed Fourth: " << double(nftkVot::diffclock(end2,begin2)) << " s";

	nftkVot::stopProgram();

	std::cout<<std::endl<<"END";
	int rert;
	std::cin>>rert;













	std::string filenameCones_1 = "output/out__votingSumVotes.tif";
	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_16 >(_votingSumVotes, filenameCones_1.c_str() )){
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
		std::string filenameCones_2as = "output/out__votingSumVotes_cone"+sas+".tif";
		if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
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
		std::string filenameCones_2as = "output/z_out__votingSumVotes_cone"+sas+".tif";
		if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2as.c_str() )){
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


	std::string filenameCones_2 = "output/out__votingSumVotes_cone.tif";
	if( nftkVot::writeImage< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_8 >(_votingSumVotes, filenameCones_2.c_str() )){
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
	////		//string filename = "output/cones/out_ImageOfVotes_"+s+".jpg";
	////		//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(imageOfConexPixels, filename.c_str() )){
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
	//			////string filename3 = "output/out_ImageOfMask_"+s3+".jpg";
	//			////if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
	//			////	std::cout<<std::endl<<"\tProblema escribiendo";
	//			////}
	//
	//			string filename5 = "output/Sum_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename5.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename6 = "output/Mask_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//			//string filename7 = "output/Sum_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			//if( nftkVot::writeImage_mhdDouble< VotingDirType_3D >(_votingSumVotes, filename7.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//			//string filename8 = "output/Mask_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			//if( nftkVot::writeImage_mhdDouble< VotingDirType_3D >(_votingMaskVotes, filename8.c_str() )){
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
	//			string filename9_max = "output/Max_stepbystep/out_ImageOfMaxs_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(imageMaxPixels, filename9_max.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//
	//			// Errase the mask (only to see step by step/by step) mas detallado todavia, ver los picos mas que todo
	//			memset(votingMaskArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//
	//			//// The amount of memory required turns out to be to high :( what to do ? 
	//			//string filename9 = "output/Direction_perpixel/out_DirPerPixel_"+s4+"__"+s5+".mhd";
	//			//if( nftkVot::writeImage_mhdDouble< VotingDirPerType >(votingDirPerPixel, filename9.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//
	//			//string filename9 = "output/Direction_perpixel/out_DirPerPixel_"+s4+"__"+s5+".mhd";
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
	//	//string filename_dir = "output/Directions/Dri_vote_"+s44+".txt";
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
	//	//string filename2 = "output/out_ImageOfSums_"+s2+".jpg";
	//	//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename2.c_str() )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	////votingSumVotes
	//	////Store Derivative in X
	//	//const char* filename = "output/out_ImageOfVotes.jpg";
	//	//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(imageOfVotingPixels, filename )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	//const char* filename2 = "output/out_vote_1.jpg";
	//	//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename2 )){
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
	//			computeCones_3D(*(it+ty));
	//		}
	//		std::cout<<"\t done";
	//		//stringstream out45;
	//		//out45<<_intSpan;
	//		//string s45 = out45.str();
	//		//string filename_dir2 = "output/Directions/Dri_vote_"+s45+".txt";
	//		//ofstream dirFilOut2;
	//		//dirFilOut2.open (filename_dir2.c_str());
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//		//{
	//		//	dirFilOut2 << it->x <<"\t"<<it->y<<"\t"<<it->angIndex<<"\t"<<it->mag<<"\n";
	//		//}
	//
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
	//				////getMaxDirection(*it);    
	//				//computeCones_3D(*it);
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
	//		////string filename3 = "output/out_ImageOfMask_"+s3+".jpg";
	//		////if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
	//		////	std::cout<<std::endl<<"\tProblema escribiendo";
	//		////}
	//		//string filename4 = "output/out_ImageOfSums_"+s3+".jpg";
	//		//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename4.c_str() )){
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
	//			string filename5 = "output/Sum_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename5.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename6 = "output/Mask_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//			string filename7 = "output/Sum_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			if( nftkVot::writeImage_mhdDouble< VotingDirType_3D >(_votingSumVotes, filename7.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename8 = "output/Mask_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			if( nftkVot::writeImage_mhdDouble< VotingDirType_3D >(_votingMaskVotes, filename8.c_str() )){
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
	//	//	computeCones_3D(_hmin, _hmax, r); // Eror pues los acabamos de calcular
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
	//	////computeCones_3D(_hmin, _hmax, 0);
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
	////	computeCones_3D(_hmin, _hmax, _radius);
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
	////		computeCones_3D(_hmin, _hmax, r);
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
	////	computeCones_3D(_hmin, _hmax, 0);
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
	//	computeCones_3D(_hmin, _hmax, _radius);
	//	computeCones_3D_3D_prob(_hmin, _hmax, _radius); // Just for one moment
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
	////		//string filename = "output/cones/out_ImageOfVotes_"+s+".jpg";
	////		//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(imageOfConexPixels, filename.c_str() )){
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
	//			////string filename3 = "output/out_ImageOfMask_"+s3+".jpg";
	//			////if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
	//			////	std::cout<<std::endl<<"\tProblema escribiendo";
	//			////}
	//
	//			string filename5 = "output/Sum_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename5.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename6 = "output/Mask_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//			//string filename7 = "output/Sum_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			//if( nftkVot::writeImage_mhdDouble< VotingDirType_3D >(_votingSumVotes, filename7.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//			//string filename8 = "output/Mask_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			//if( nftkVot::writeImage_mhdDouble< VotingDirType_3D >(_votingMaskVotes, filename8.c_str() )){
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
	//			string filename9_max = "output/Max_stepbystep/out_ImageOfMaxs_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(imageMaxPixels, filename9_max.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//
	//			// Errase the mask (only to see step by step/by step) mas detallado todavia, ver los picos mas que todo
	//			memset(votingMaskArray, 0,  npix*sizeof(VotingDirType_3D::PixelType));
	//
	//			//// The amount of memory required turns out to be to high :( what to do ? 
	//			//string filename9 = "output/Direction_perpixel/out_DirPerPixel_"+s4+"__"+s5+".mhd";
	//			//if( nftkVot::writeImage_mhdDouble< VotingDirPerType >(votingDirPerPixel, filename9.c_str() )){
	//			//	std::cout<<std::endl<<"\tProblema escribiendo";
	//			//}
	//
	//			//string filename9 = "output/Direction_perpixel/out_DirPerPixel_"+s4+"__"+s5+".mhd";
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
	//	//string filename_dir = "output/Directions/Dri_vote_"+s44+".txt";
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
	//	//string filename2 = "output/out_ImageOfSums_"+s2+".jpg";
	//	//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename2.c_str() )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	////votingSumVotes
	//	////Store Derivative in X
	//	//const char* filename = "output/out_ImageOfVotes.jpg";
	//	//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(imageOfVotingPixels, filename )){
	//	//	std::cout<<std::endl<<"\tProblema escribiendo";
	//	//}
	//
	//	//const char* filename2 = "output/out_vote_1.jpg";
	//	//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename2 )){
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
	//			computeCones_3D(*(it+ty));
	//		}
	//		std::cout<<"\t done";
	//		//stringstream out45;
	//		//out45<<_intSpan;
	//		//string s45 = out45.str();
	//		//string filename_dir2 = "output/Directions/Dri_vote_"+s45+".txt";
	//		//ofstream dirFilOut2;
	//		//dirFilOut2.open (filename_dir2.c_str());
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//		//{
	//		//	dirFilOut2 << it->x <<"\t"<<it->y<<"\t"<<it->angIndex<<"\t"<<it->mag<<"\n";
	//		//}
	//
	//		//for(std::vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
	//				////getMaxDirection(*it);    
	//				//computeCones_3D(*it);
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
	//		////string filename3 = "output/out_ImageOfMask_"+s3+".jpg";
	//		////if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
	//		////	std::cout<<std::endl<<"\tProblema escribiendo";
	//		////}
	//		//string filename4 = "output/out_ImageOfSums_"+s3+".jpg";
	//		//if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename4.c_str() )){
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
	//			string filename5 = "output/Sum_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingSumVotes, filename5.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename6 = "output/Mask_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//			if( nftkVot::writeImage< VotingDirType_3D, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//
	//			string filename7 = "output/Sum_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			if( nftkVot::writeImage_mhdDouble< VotingDirType_3D >(_votingSumVotes, filename7.c_str() )){
	//				std::cout<<std::endl<<"\tProblema escribiendo";
	//			}
	//			string filename8 = "output/Mask_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//			if( nftkVot::writeImage_mhdDouble< VotingDirType_3D >(_votingMaskVotes, filename8.c_str() )){
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
	//	//	computeCones_3D(_hmin, _hmax, r); // Eror pues los acabamos de calcular
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
	//	////computeCones_3D(_hmin, _hmax, 0);
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
	////	computeCones_3D(_hmin, _hmax, _radius);
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
	////		computeCones_3D(_hmin, _hmax, r);
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
	////	computeCones_3D(_hmin, _hmax, 0);
	////}
	////
}

// ############################################################################################################################################################################
int ftkVoting_3D::computeDirecIndex_3D(double dx, double dy, double dz) const
{
	double maxmax_dot = 0;
	int maxmax_pos = 0;
	double dotprod = 0;

	for( int uu=0; uu<_NN_dir; ++uu ) // All possible directions
	{
		dotprod = _conesPru_3D_new[uu].dxx*dx+_conesPru_3D_new[uu].dyy*dy+_conesPru_3D_new[uu].dzz*dz;
		if ( maxmax_dot<dotprod)
		{
			maxmax_dot = dotprod;
			maxmax_pos = uu;
		}
	}
	return maxmax_pos;
}

// ############################################################################################################################################################################
void ftkVoting_3D::updateCones()
{

}

//// ############################################################################################################################################################################
void inline ftkVoting_3D::updateDirection_3D(VPoint3D& vp)
{
	VotingDirType_3D::PixelType * votingSumArray_pos = _votingSumVotes->GetBufferPointer() + vp.pos;

	double maxx = -100;
	double posMax = 0;

	int temp_vv_max = 0; // Maxima directino 
	int scale_max = _hmax;
	
	ftkBins3D * poin_bin;

	int dir_vote = vp.direc_vote;
	for( int angle_int = 0; angle_int<_intSpan; ++angle_int ) // The different angles (1-10)
	{
// 		for( int vv = 0; vv<_voteDirec_3D_new_notqu[dir_vote][angle_int].size(); ++vv )	// All the bins at a given dot product distance from the given direction (dir_vote)
// 		{
// 			int temp_vv = _voteDirec_3D_new_notqu[dir_vote][angle_int][vv];
// 			for( int bin_cont=0; bin_cont < _conesPru_3D_new_notqu[temp_vv].size(); ++bin_cont )
// 			{
// 				posMax = votingSumArray_pos[_conesPru_3D_new_notqu[temp_vv][bin_cont].off];
// 				if( posMax > maxx )
// 				{
// 					maxx = posMax;
// 					temp_vv_max = temp_vv;
// 				}
// 
// 			}
// 		}
		
		for( int vv = 0; vv<_voteDirec_3D_new[dir_vote][angle_int].size(); ++vv )	// All the bins at a given dot product distance from the given direction (dir_vote)
		{
			int temp_vv = _voteDirec_3D_new[dir_vote][angle_int][vv];
			for( int bin_cont=0; bin_cont < vp.scale/*_conesPru_3D_new[temp_vv].size()*/; ++bin_cont )
			{
				poin_bin =  &_conesPru_3D_new[temp_vv][bin_cont];
				
				for( int votes=0; votes < (*poin_bin).size(); ++votes )
				{
					posMax = votingSumArray_pos[(*poin_bin)[votes].off] = votingSumArray_pos[(*poin_bin)[votes].off];
					if( posMax > maxx )
					{
						maxx = posMax;
						temp_vv_max = temp_vv;
						scale_max = bin_cont;
					}
				}
			}
		}
	}
	
	//std::cout<<std::endl<<scale_max;

	// If there is not maximum continue with the same gradient ??? !!! (es decir no mas grande que un determinado threshdol)
	vp.direc_vote = temp_vv_max;
	
 	//// If this is less than the maximum otherwise do not change
 	//if(scale_max+(int)((double)scale_max/4)<=_hmax)
 	//{
 	//	vp.scale = scale_max+(int)((double)scale_max/4); // ??? Give the third part (option relajada)
 	//}

 	// If this is less than the maximum otherwise do not change
 	if(vp.scale = scale_max + 2<=_hmax)
 	{
 		vp.scale = scale_max + 2;
 	}
	
	
	

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