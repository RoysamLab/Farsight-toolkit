// ############################################################################################################################################################################
#include "ftkVoting.h"


#include <cmath>

int ftkWPoint2D::nx = 0;
int ftkWPoint2D::ny = 0;

int ftkCone2D::nx = 0;
int ftkCone2D::ny = 0;
int ftkCone2D::contador_1 = 0;

int ftkBins2D::nx = 0;
int ftkBins2D::ny = 0;

// Static constants
const double ftkVoting::epsilon((double)0.001);
const double ftkVoting::pi = (double)3.1415927;

// ############################################################################################################################################################################
ftkVoting::ftkVoting(){

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

}

//// ############################################################################################################################################################################
void ftkVoting::setParams(	int hmin, int hmax,	int radius,	double min_grad, double scale ){
	_hmin = hmin;
	_hmax = hmax;
	_radius = radius;
	_min_grad = min_grad;
	_scale = scale; // Scale for computing the gradient using DoG
}

// ############################################################################################################################################################################
void ftkVoting::compute(nftkVot::InputImageType::Pointer inputImage)
{
	nx = inputImage->GetLargestPossibleRegion().GetSize()[0];
	ny = inputImage->GetLargestPossibleRegion().GetSize()[1];
	npix = nx*ny;

	_voting_points = vector<VPoint2D>(); // This are the points that actually vote
	_voting_points.reserve(npix/2);		// To speed up, reserve some space (on average what will be used)


	// Derivative XXXXXXXXXX
	VotingDirType::Pointer votingDirX = VotingDirType::New();
	votingDirX->SetRegions( inputImage->GetRequestedRegion() );
	votingDirX->Allocate();

	itk::SobelOperator< VotingDirPixelType, 2 > sobelOperator;
	sobelOperator.SetDirection( 0 );
	sobelOperator.CreateDirectional();
	typedef itk::ImageRegionIterator< VotingDirType > IteratorType3;
	IteratorType3 out2( votingDirX, votingDirX->GetRequestedRegion() );

	typedef itk::ConstNeighborhoodIterator< nftkVot::InputImageType > NeighborhoodIteratorType2;
	typedef itk::ImageRegionIterator< nftkVot::InputImageType > IteratorType2;
	NeighborhoodIteratorType2::RadiusType radius2 = sobelOperator.GetRadius();
	NeighborhoodIteratorType2 it2( radius2, inputImage, inputImage->GetRequestedRegion() );
	itk::NeighborhoodInnerProduct< nftkVot::InputImageType > innerProduct;

	for (it2.GoToBegin(), out2.GoToBegin(); !it2.IsAtEnd(); ++it2, ++out2)
	{
		out2.Set( innerProduct( it2, sobelOperator ) );
	}

	// Derivative YYYYYYYYYY
	VotingDirType::Pointer votingDirY = VotingDirType::New();
	votingDirY->SetRegions( inputImage->GetRequestedRegion() );
	votingDirY->Allocate();

	itk::SobelOperator<VotingDirPixelType, 2> sobelOperatorY;
	sobelOperatorY.SetDirection( 1 );
	sobelOperatorY.CreateDirectional();
	typedef itk::ImageRegionIterator< VotingDirType > IteratorType4;
	IteratorType4 out3( votingDirY, votingDirY->GetRequestedRegion() );

	typedef itk::ConstNeighborhoodIterator< nftkVot::InputImageType > NeighborhoodIteratorType3;
	typedef itk::ImageRegionIterator< nftkVot::InputImageType > IteratorType3;
	NeighborhoodIteratorType3::RadiusType radius3 = sobelOperatorY.GetRadius();
	NeighborhoodIteratorType3 it3( radius3, inputImage, inputImage->GetRequestedRegion() );
	itk::NeighborhoodInnerProduct< nftkVot::InputImageType > innerProduct2;

	for (it3.GoToBegin(), out3.GoToBegin(); !it3.IsAtEnd(); ++it3, ++out3)
	{
		out3.Set( innerProduct2( it3, sobelOperatorY ) );
	}

	// Magnitude Image
	VotingDirType::Pointer votingMagImage = VotingDirType::New(); //VotingDirType = double
	votingMagImage->SetRegions( inputImage->GetRequestedRegion()); // IMPORTANTE PARA CREAR UNA IMAGEN NUEVA EN BASE A UNA QUE YA EXISTE EN VEZ DE PONERME A LEER LOS TAMANOS Y LAS REGIONES
	votingMagImage->Allocate();

	typedef itk::ImageRegionIteratorWithIndex< VotingDirType > ITVotingMag;
	ITVotingMag iVotingMag(votingMagImage, votingMagImage->GetLargestPossibleRegion() );

	typedef itk::ImageRegionIteratorWithIndex< VotingDirType > ITVotingDir;
	ITVotingDir iVotingDirX(votingDirX, votingDirX->GetLargestPossibleRegion() );
	ITVotingDir iVotingDirY(votingDirY, votingDirY->GetLargestPossibleRegion() );

	for ( iVotingMag.GoToBegin(),iVotingDirX.GoToBegin(),iVotingDirY.GoToBegin(); !iVotingMag.IsAtEnd(); ++iVotingMag, ++iVotingDirX, ++iVotingDirY ){
		iVotingMag.Set(sqrt( iVotingDirX.Get()*iVotingDirX.Get() + iVotingDirY.Get()*iVotingDirY.Get() ));
		if(iVotingMag.Get()>epsilon){
			iVotingDirX.Set(iVotingDirX.Get()/iVotingMag.Get());
			iVotingDirY.Set(iVotingDirY.Get()/iVotingMag.Get());
		}
	}


	// Magnitude Image Binary (just to store in file the voting pixels)
	VotingDirType::Pointer votingMagImage_bin = VotingDirType::New();
	votingMagImage_bin->SetRegions( inputImage->GetRequestedRegion() );
	votingMagImage_bin->Allocate();

	// Canny edge detection
	typedef itk::CannyEdgeDetectionImageFilter <VotingDirType, VotingDirType> CannyEdgeDetectionImageFilterType;
	CannyEdgeDetectionImageFilterType::Pointer cannyFilter = CannyEdgeDetectionImageFilterType::New();
	cannyFilter->SetInput(inputImage);
	cannyFilter->SetVariance( 2.5 );
	cannyFilter->SetUpperThreshold( 0.03/*0.0238*/ );
	cannyFilter->SetLowerThreshold( 0.005/*0.0175*/ );
	cannyFilter->Update();

	// Save the canny edge result
	string filenameCanny = "output/out_ImageOfCanny.jpg";
	if( nftkVot::writeImage< VotingDirType, OutputImageType >(cannyFilter->GetOutput(), filenameCanny.c_str() )){
		cout<<endl<<"\tProblema escribiendo";
	}



	// Arrays of data
	nftkVot::InputImageType::PixelType * votingImaArray = inputImage->GetBufferPointer();
	VotingDirType::PixelType * votingDirXArray = votingDirX->GetBufferPointer();
	VotingDirType::PixelType * votingDirYArray = votingDirY->GetBufferPointer();
	VotingDirType::PixelType * votingMagArray = votingMagImage->GetBufferPointer();
	VotingDirType::PixelType * votingMagArray_bin = votingMagImage_bin->GetBufferPointer();
	VotingDirType::PixelType * votingCannyArray = cannyFilter->GetOutput()->GetBufferPointer();


	//double maxgrad = votingMagImag

	typedef itk::MinimumMaximumImageCalculator < VotingDirType > ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New ();
	imageCalculatorFilter->SetImage( votingMagImage );
	imageCalculatorFilter->Compute();

	VotingDirType::PixelType maxVotVal = imageCalculatorFilter->GetMaximum();


	cout<<endl<<"Max Grad: "<<maxVotVal<<", Thresh: "<<_min_grad;

	//_OriginalScale = maxVotVal;

	// Scale the values of the gradient by the maximum value of the gradiente, to ensure that the maximum value is 1
	for(int i=0; i<npix; i++) {
		if (votingMagArray[i]</*0.001*/_min_grad)
		{
			votingMagArray[i] = 0;
		}
		else
		{
			//votingMagArray[i] /= maxVotVal;
		}
	}


	// Testing to store the resulting voting image base on the gradient
	// Save the canny edge result
	string filenameGradVot = "output/out_ImageOfGradVot.jpg";
	if( nftkVot::writeImage< VotingDirType, OutputImageType >(votingMagImage, filenameGradVot.c_str() )){
		cout<<endl<<"\tProblema escribiendo";
	}
	// Put 1 if mag array is different to zero
	for(int i=0; i<npix; i++) {
		votingMagArray_bin[i] = 0;
		if (votingMagArray[i]!=0)
		{
			votingMagArray_bin[i] = 1;
		}
	}
	string filenameGradVot_bin = "output/out_ImageOfGradVot_bin.jpg";
	if( nftkVot::writeImage< VotingDirType, OutputImageType >(votingMagImage_bin, filenameGradVot_bin.c_str() )){
		cout<<endl<<"\tProblema escribiendo";
	}
	// Put 1 if mag array intersect with canny edge is di
	for(int i=0; i<npix; i++) {
		votingMagArray_bin[i] = 0;
		if (votingMagArray[i]!=0 && votingCannyArray[i]!=0)
		{
			votingMagArray_bin[i] = 1;
		}
	}
	string filenameGradVotInterCanny = "output/out_ImageOfGradVotInterCanny.jpg";
	if( nftkVot::writeImage< VotingDirType, OutputImageType >(votingMagImage_bin, filenameGradVotInterCanny.c_str() )){
		cout<<endl<<"\tProblema escribiendo";
	}






	// Padding
	bw = sqrt((double)(_radius*_radius + _hmax*_hmax)) + 3; // Que es esto ??
	const int bw2 = 2*bw;


	// Calculate the voting poitns
	//int stop1;
	int indx;
	VPoint2D vp;
	int x, y, k;
	for(y = 0; y < ny; y++)
	{
		for(x = 0; x < nx; x++)
		{
			k = x+nx*y;
			if (votingMagArray[k] > epsilon) // This is the original one, only vote if the mag of gradient is greater than epsilon (and gra_min) previously thresholded
				//if(votingCannyArray[k]!=0) // This is new, use the canny edge detection to vote form this values (IT DOES NOT WORK TO USE THE CANNY EDGE DETECTION, INCREIBLE)
			{
				if ((indx=computeAngleIndex(votingDirXArray[k], votingDirYArray[k])) >= 0 ) // NO ESTAN INCLUYENDO LOS CEROS PORQUE ??
				{
					vp.x = x + bw;
					vp.y = y + bw;
					vp.mag = votingMagArray[k]; 
					vp.angIndex = indx;
					_voting_points.push_back(vp);
				}	
			}
		}
	}
	//cin>>stop1;









	nx += bw2;
	ny += bw2;
	npix = nx*ny;

	ftkWPoint2D::setImageSize(nx,ny);
	ftkBins2D::setImageSize(nx,ny);
	ftkCone2D::setImageSize(nx,ny);

	vote();

	// Switch back to the original size
	nx -= bw2;
	ny -= bw2;
	npix = nx*ny;

	//_votingVotes = VotingDirType::New();
	//_votingVotes->SetRegions( inputImage->GetRequestedRegion() );
	//_votingVotes->Allocate();

	// Copy sum votes to 


}


//void ftkVoting::compute(const Image& I)
//{
//	//cout<<"Executing local version of circular voting "<<" Thresh :"<<_threshold<<endl;
//	nx = I.nx();
//	ny = I.ny();
//	npix = nx*ny;
//
//	_voting_points = vector<VPoint2D>();
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
//			//cout << "#Voting points = " << _voting_points.size() << endl;
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
//				//cout<<"Voting Threshold "<<_threshold<<" max Value "<<maxval<<" thresh * maxval:" <<_threshold * maxval<<endl;
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
//				vector<RegionFeature> rf = regionFeatures(label, nc, _vote);
//				_vote.reset(0);
//				_center_weights.clear();
//				for(vector<RegionFeature>::iterator it=rf.begin(); it!=rf.end(); it++) {
//					if (_vote.inImage(it->xc, it->yc)) {
//						_vote(it->xc, it->yc) = it->mean;
//						_centers.push_back(Point( it->xc, it->yc ));
//						_center_weights.push_back(it->mean);
//					}   
//				}
//
//				for(vector<VPoint2D>::iterator it=_voting_points.begin(); it!=_voting_points.end(); it++) {
//					it->x -= bw;
//					it->y -= bw;
//					it->xc -= bw;
//					it->yc -= bw;
//				}
//}
//
// ############################################################################################################################################################################
ftkVoting::~ftkVoting()
{
}

// ############################################################################################################################################################################
void ftkVoting::setPrefix(const string& p)
{
	_prefix = p;
}

// ############################################################################################################################################################################
int ftkVoting::computeAngleIndex(double dx, double dy) const 
{
	double a = acos(dx);
	if (dy<0) {
		a = 2*pi-a;
	}

	//int indx_theta = nftkVot::round_double(a/delta_theta); // De por si esta divisino no veo como pueda ser mas de 255
	//return max(0, min(ntheta-1, indx_theta)); // Se asegura que el maximo devuelto sea 255 (256 nunca va a salir de aca)

	int indx_theta = floor(a/delta_theta); // De por si esta divisino no veo como pueda ser mas de 255
	if(indx_theta<0)
	{int rr;
	cout<<endl<<"Error en compute angle: "<<indx_theta<<" "<<dx<<" "<<dy;
	cin>>rr;
	}
	return max(0, min(ntheta-1, indx_theta));
}

// ############################################################################################################################################################################
void inline ftkCone2D::vote(VotingDirType::PixelType * p, const VPoint2D& vp)
{
	//int stop2;
	iterator from = begin(); //Pointer to Bin
	iterator to = end();     
	for(iterator it = from; it != to; it++) 
	{
		for(ftkBins2D::iterator it1=it->begin(); it1!=it->end(); ++it1) {
			p[it1->off] += vp.mag * it1->w;
			//p[it1->off] += 1; //Numer of votes, not using magnitude
			//contador_1++;
		}      
	}
}

//// ############################################################################################################################################################################
//void inline ftkCone2D::vote_dir(VotingDirType::PixelType * p, const VPoint2D& vp) //Vota y guarda direcciones de los votos
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
void inline ftkCone2D::vote(VotingDirType::PixelType * p, const VPoint2D& vp, int& dist)
{

	iterator binRequired = begin();
	ftkBins2D::iterator from = binRequired[dist].begin();
	ftkBins2D::iterator to = binRequired[dist].end();

	for(ftkBins2D::iterator it = from; it != to; it++) 
	{
		p[it->off] += vp.mag * it->w; 
		//p[it->off] += 1; // Number of votes, not using magnitude
		//contador_1++;
	}
}


// I HAVE TO FINISH THIS ONE, VOTES WITH PROPABILITY AND WITHOUT DISTANCE (FASTER VERSION)
//// ############################################################################################################################################################################
//void inline ftkCone2D::vote(VotingDirType::PixelType * p, const VPoint2D& vp) // Incluye la probabilidad
//{
//
//	iterator binRequired = begin();
//	ftkBins2D::iterator from = binRequired[dist].begin();
//	ftkBins2D::iterator to = binRequired[dist].end();
//
//	for(ftkBins2D::iterator it = from; it != to; it++) 
//	{
//		//p[it->off] += vp.mag;// * (1/(double)mag)/*it->w*/; 
//		p[it->off] += vp.mag * (1/(double)mag); 
//		//p[it->off] += (1/(double)mag); 
//		//p[it->off] += (1/(double)mag)/*it->w*/; 
//		//cout<<" "<<mag;
//		//p[it->off] += 1; // Number of votes, not using magnitude
//		contador_1++;
//	}
//
//	//int stop2;
//	iterator from = begin(); //Pointer to Bin
//	iterator to = end();     
//	for(iterator it = from; it != to; it++) 
//	{
//		for(ftkBins2D::iterator it1=it->begin(); it1!=it->end(); ++it1) {
//			p[it1->off] += vp.mag * (1/(double)_voteDirec_prob[].second[dist]);
//			//p[it1->off] += 1; //Numer of votes, not using magnitude
//			contador_1++;
//		}      
//	}
//}

// ############################################################################################################################################################################
void inline ftkCone2D::vote(VotingDirType::PixelType * p, const VPoint2D& vp, int& dist, int &mag) // Incluye la probabilidad
{

	iterator binRequired = begin();
	ftkBins2D::iterator from = binRequired[dist].begin();
	ftkBins2D::iterator to = binRequired[dist].end();

	for(ftkBins2D::iterator it = from; it != to; it++) 
	{
		//p[it->off] += vp.mag;// * (1/(double)mag)/*it->w*/; 
		p[it->off] += vp.mag * (1/(double)mag); 
		//p[it->off] += (1/(double)mag); 
		//p[it->off] += (1/(double)mag)/*it->w*/; 
		//cout<<" "<<mag;
		//p[it->off] += 1; // Number of votes, not using magnitude
		//contador_1++;
	}
}

// ############################################################################################################################################################################
void inline ftkCone2D::vote_dir(vector<vector<int> >& p_dir, VotingDirType::PixelType * p, const VPoint2D& vp, int& dist, int& offset_1)
{

	iterator binRequired = begin();
	ftkBins2D::iterator from = binRequired[dist].begin();
	ftkBins2D::iterator to = binRequired[dist].end();

	for(ftkBins2D::iterator it = from; it != to; it++) 
	{
		p[it->off] += vp.mag * it->w;
		p_dir.at(it->off+offset_1).at(vp.angIndex) += 1; // angIndex 0 255
		//p_dir[it->off][vp.angIndex] += 1;
		//p[it->off] += 1; // Number of votes, not using magnitude
		//contador_1++;
	}
}

//// ############################################################################################################################################################################
void inline ftkVoting::updateDirection(VPoint2D& vp)
{
	ftkBins2D::iterator max_point;
	double maxvote = 0;
	for( int tt=_voteDirec[vp.angIndex].first; tt!=_voteDirec[vp.angIndex].second; tt=(tt+1)%ntheta )
	{
		ftkCone2D::iterator from = _conesPru[tt].begin();
		ftkCone2D::iterator to = _conesPru[tt].end();
		//_conesPru[tt].vote(p,vp);

		VotingDirType::PixelType * votingSumArray = _votingSumVotes->GetBufferPointer()+vp.pos;
		for(ftkCone2D::iterator it = from; it != to; it++) 
		{
			for(ftkBins2D::iterator it1=it->begin(); it1!=it->end(); it1++)
			{
				if ( votingSumArray[it1->off]>maxvote) 
				{
					maxvote = max(maxvote, votingSumArray[it1->off]);
					max_point = it1;
				}    
			}
		} 
	}

	if (maxvote<epsilon)
		return;

	double dx = max_point->x;
	double dy = max_point->y;

	double r = sqrt(dx*dx+dy*dy);
	if (r>epsilon) {
		vp.xc = vp.x + max_point->x;
		vp.yc = vp.y + max_point->y;
		vp.angIndex = computeAngleIndex(dx/r, dy/r);
		//cout<<"\n"<<vp.angIndex;
	}
}


//// ############################################################################################################################################################################
void inline ftkVoting::updateDirection_prob(VPoint2D& vp)
{
	ftkBins2D::iterator max_point;
	double maxvote = 0;
	for( int tt=_voteDirec[vp.angIndex].first; tt!=_voteDirec[vp.angIndex].second; tt=(tt+1)%ntheta ) // Update direction for each voting point
	{
		ftkCone2D::iterator from = _conesPru[tt].begin();
		ftkCone2D::iterator to = _conesPru[tt].end();
		//_conesPru[tt].vote(p,vp);

		VotingDirType::PixelType * votingSumArray = _votingSumVotes->GetBufferPointer()+vp.pos;
		for(ftkCone2D::iterator it = from; it != to; it++)  // Iterate to the cones
		{
			for(ftkBins2D::iterator it1=it->begin(); it1!=it->end(); it1++) // Iterate to the bins of each cone
			{
				if ( votingSumArray[it1->off]>maxvote) 
				{
					maxvote = max(maxvote, votingSumArray[it1->off]);
					max_point = it1;
				}    
			}
		} 
	}

	if (maxvote<epsilon)
		return;

	double dx = max_point->x;
	double dy = max_point->y;

	double r = sqrt(dx*dx+dy*dy);
	if (r>epsilon) {
		vp.xc = vp.x + max_point->x;
		vp.yc = vp.y + max_point->y;
		vp.angIndex = computeAngleIndex(dx/r, dy/r);
		//cout<<"\n"<<vp.angIndex;
	}
}


//// ############################################################################################################################################################################
void ftkVoting::computeCones(int hmin, int hmax, int radius)
{

	int stop11;

	_conesPru = vector<ftkCone2D>(ntheta); // Todos los conos posibles
	for( int uu=0; uu<ntheta; uu++ )
	{
		for( int uuu=0; uuu<hmax-hmin+1; uuu++ )
		{
			ftkBins2D bin;
			_conesPru[uu].push_back(bin);
		}
	}

	int y1;
	int R;
	double R_dou;
	int R_quan;
	int ang_quan;
	double x_nor, y_nor;
	double x_nor2, y_nor2;
	int countt=0;

	ftkWPoint2D wp;
#pragma omp parallel for num_threads(8) private(y1)
	for( int xx=1; xx<=hmax; xx++ )
	{
		y1 = (int)floor((double)sqrt((double)hmax*hmax-(double)xx*xx));
#pragma omp parallel for num_threads(8) private(wp,R_dou,R,R_quan,x_nor,y_nor,ang_quan)
		for( int yy=1; yy<=y1; yy++ )
		{
			R_dou = (double)sqrt((double)yy*yy+(double)xx*xx);
			R = (int)ceil(R_dou);
			if( R>= hmin) // Can be done more efficiently but for now is ok the speed
			{

				R_quan = R-hmin;
				//cout<<endl<<R_quan;
				x_nor = ((double)xx)/R_dou;
				y_nor = ((double)yy)/R_dou;
				ang_quan = computeAngleIndex(x_nor, y_nor);
				//ftkWPoint2D wp;
				wp.x = xx;
				wp.y = yy;
				wp.w = 1; // In case of deciding to put some weight
				//cout<<endl<<ang_quan<<" "<<R_quan;
				_conesPru[ang_quan][R_quan].push_back(wp);

				x_nor2 = -x_nor;
				y_nor2 = y_nor;
				ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 64, + 128 etc...
				wp.x = -xx;
				wp.y = yy;
				wp.w = 1; // In case of deciding to put some weight
				_conesPru[ang_quan][R_quan].push_back(wp);

				x_nor2 = x_nor;
				y_nor2 = -y_nor;
				ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 64, + 128 etc...
				wp.x = xx;
				wp.y = -yy;
				wp.w = 1; // In case of deciding to put some weight
				_conesPru[ang_quan][R_quan].push_back(wp);

				x_nor2 = -x_nor;
				y_nor2 = -y_nor;
				ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 64, + 128 etc...
				wp.x = -xx;
				wp.y = -yy;
				wp.w = 1; // In case of deciding to put some weight
				_conesPru[ang_quan][R_quan].push_back(wp);

				//cout<<endl<<"coutn: "<<countt++;

				//cout<<endl<<" "<<R_quan<<" "<<ang_quan<<" "<<R_dou<<"-> "<<xx<<","<<yy<<"|"<<-xx<<" "<<yy<<"|"<<xx<<" "<<-yy<<"|"<<-xx<<" "<<-yy;
			}
		}
	}
	//ftkWPoint2D wp;
#pragma omp parallel for num_threads(16) private(wp,R_dou,R_quan,x_nor,y_nor,ang_quan)
	for( int xx=hmin; xx<=hmax; xx++ )
	{

		R_dou = xx;
		R_quan = R_dou-hmin;
		x_nor = ((double)xx)/R_dou;
		y_nor = 0;
		ang_quan = computeAngleIndex(x_nor, y_nor);
		//ftkWPoint2D wp;
		wp.x = xx;
		wp.y = 0;
		wp.w = 1; // In case of deciding to put some weight
		_conesPru[ang_quan][R_quan].push_back(wp);

		x_nor = -((double)xx)/R_dou;
		y_nor = 0;
		ang_quan = computeAngleIndex(x_nor, y_nor);
		wp.x = -xx;
		wp.y = 0;
		wp.w = 1; // In case of deciding to put some weight
		_conesPru[ang_quan][R_quan].push_back(wp);

		x_nor = 0;
		y_nor = ((double)xx)/R_dou;
		ang_quan = computeAngleIndex(x_nor, y_nor);
		wp.x = 0;
		wp.y = xx;
		wp.w = 1; // In case of deciding to put some weight
		_conesPru[ang_quan][R_quan].push_back(wp);

		x_nor = 0;
		y_nor = -((double)xx)/R_dou;
		ang_quan = computeAngleIndex(x_nor, y_nor);
		wp.x = 0;
		wp.y = -xx;
		wp.w = 1; // In case of deciding to put some weight
		_conesPru[ang_quan][R_quan].push_back(wp);


		//cout<<endl<<"coutn: "<<countt++;

		//cout<<endl<<" "<<R_quan<<" "<<ang_quan<<" "<<R_dou<<"-> "<<xx;

	}


	//cout<<"Terminio: ";
	//int ght;
	//cin>>ght;


	//countt = 0;
	for( int tt=0; tt<_conesPru.size(); tt++ )
	{
		_conesPru[tt].setOffset();
	}

	// Now calculate the span
	double span = sqrt((double)hmax*hmax+(double)radius*radius);
	double dx_dou = hmax/span;
	double dy_dou = radius/span;
	_intSpan = computeAngleIndex(dx_dou, dy_dou);
	_voteDirec = vector< pair< int,int > > (ntheta);
	// This copy is for a test on the voting, just to use the last state of the voting directions (updated gradients) and perform actually a vote in the biggest span
	_voteDirec_copy = vector< pair< int,int > > (ntheta);
	for( int tt=0; tt<ntheta; tt++ )
	{
		pair<int,int> spann;
		if(tt-_intSpan<0)
		{
			spann.first = (tt-_intSpan+ntheta)%ntheta;
		}
		else
		{
			spann.first = (tt-_intSpan)%ntheta;
		}
		//cout<<endl<<"\t"<<spann.first;
		spann.second = (tt+_intSpan)%ntheta;
		_voteDirec[tt] = spann;
	}
	// Copy the initial votes direction (span)
	_voteDirec_copy = _voteDirec;



	cout<<"Cones computed";
	//cin>>stop11;

}


//// ############################################################################################################################################################################
void ftkVoting::computeCones_prob(int hmin, int hmax, int radius)
{

	int stop11;

	_conesPru_prob = vector<ftkCone2D>(ntheta); // Todos los conos posibles
	for( int uu=0; uu<ntheta; uu++ )
	{
		for( int uuu=0; uuu<hmax-hmin+1; uuu++ )
		{
			ftkBins2D bin;
			_conesPru_prob[uu].push_back(bin);
		}
	}

	int y1;
	int R;
	double R_dou;
	int R_quan;
	int ang_quan;
	double x_nor, y_nor;
	double x_nor2, y_nor2;
	int countt=0;

	ftkWPoint2D wp;
	//#pragma omp parallel for num_threads(8) private(y1)
	for( int xx=1; xx<=hmax; xx++ )
	{
		y1 = (int)floor((double)sqrt((double)hmax*hmax-(double)xx*xx));
		//#pragma omp parallel for num_threads(8) private(wp,R_dou,R,R_quan,x_nor,y_nor,ang_quan)
		for( int yy=1; yy<=y1; yy++ )
		{
			R_dou = (double)sqrt((double)yy*yy+(double)xx*xx);
			R = (int)ceil(R_dou);
			if( R>= hmin) // Can be done more efficiently but for now is ok the speed
			{

				R_quan = R-hmin;
				//cout<<endl<<R_quan;
				x_nor = ((double)xx)/R_dou;
				y_nor = ((double)yy)/R_dou;
				ang_quan = computeAngleIndex(x_nor, y_nor);
				//ftkWPoint2D wp;
				wp.x = xx;
				wp.y = yy;
				wp.w = 1; // In case of deciding to put some weight
				//cout<<endl<<ang_quan<<" "<<R_quan;
				_conesPru_prob[ang_quan][R_quan].push_back(wp);

				x_nor2 = -x_nor;
				y_nor2 = y_nor;
				ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 64, + 128 etc...
				wp.x = -xx;
				wp.y = yy;
				wp.w = 1; // In case of deciding to put some weight
				_conesPru_prob[ang_quan][R_quan].push_back(wp);

				x_nor2 = x_nor;
				y_nor2 = -y_nor;
				ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 64, + 128 etc...
				wp.x = xx;
				wp.y = -yy;
				wp.w = 1; // In case of deciding to put some weight
				_conesPru_prob[ang_quan][R_quan].push_back(wp);

				x_nor2 = -x_nor;
				y_nor2 = -y_nor;
				ang_quan = computeAngleIndex(x_nor2, y_nor2); // It is like ang_quan + 64, + 128 etc...
				wp.x = -xx;
				wp.y = -yy;
				wp.w = 1; // In case of deciding to put some weight
				_conesPru_prob[ang_quan][R_quan].push_back(wp);

				//cout<<endl<<"coutn: "<<countt++<<" BBB";

				//cout<<endl<<" "<<R<<" "<<R_quan<<" "<<ang_quan<<" "<<R_dou<<"-> "<<xx<<","<<yy<<"|"<<-xx<<" "<<yy<<"|"<<xx<<" "<<-yy<<"|"<<-xx<<" "<<-yy<<" ACA";
			}
		}
	}
	//cout<<endl<<"gggg";
	//ftkWPoint2D wp;
	//#pragma omp parallel for num_threads(16) private(wp,R_dou,R_quan,x_nor,y_nor,ang_quan)
	for( int xx=hmin; xx<=hmax; xx++ )
	{

		R_dou = xx;
		R_quan = R_dou-hmin;
		x_nor = ((double)xx)/R_dou;
		y_nor = 0;
		ang_quan = computeAngleIndex(x_nor, y_nor);
		//ftkWPoint2D wp;
		wp.x = xx;
		wp.y = 0;
		wp.w = 1; // In case of deciding to put some weight
		_conesPru_prob[ang_quan][R_quan].push_back(wp);

		x_nor = -((double)xx)/R_dou;
		y_nor = 0;
		ang_quan = computeAngleIndex(x_nor, y_nor);
		wp.x = -xx;
		wp.y = 0;
		wp.w = 1; // In case of deciding to put some weight
		_conesPru_prob[ang_quan][R_quan].push_back(wp);

		x_nor = 0;
		y_nor = ((double)xx)/R_dou;
		ang_quan = computeAngleIndex(x_nor, y_nor);
		wp.x = 0;
		wp.y = xx;
		wp.w = 1; // In case of deciding to put some weight
		_conesPru_prob[ang_quan][R_quan].push_back(wp);

		x_nor = 0;
		y_nor = -((double)xx)/R_dou;
		ang_quan = computeAngleIndex(x_nor, y_nor);
		wp.x = 0;
		wp.y = -xx;
		wp.w = 1; // In case of deciding to put some weight
		_conesPru_prob[ang_quan][R_quan].push_back(wp);

		//cout<<endl<<"coutn: "<<countt++;

		//cout<<endl<<" "<<R_quan<<" "<<ang_quan<<" "<<R_dou<<"-> "<<xx;

	}


	cout<<endl<<"Computing probability";

	//for( int tt=0; tt<_conesPru_prob.size(); tt++ )
	//{
	//	cout<<endl<<"tamanos: "<<tt<<" "<<_conesPru_prob[tt].size()<<" -> ";
	//	for( int uu=0; uu<_conesPru_prob[tt].size(); ++uu )
	//	{
	//		cout<<" "<<_conesPru_prob[tt][uu].size();
	//	}
	//	
	//}

	//////////////////////////////////////////////
	//cout<<"Terminio de computar conos de probabilidad: ";
	//int ght;
	//cin>>ght;
	//////////////////////////////////////////////


	//countt = 0;
	for( int tt=0; tt<_conesPru_prob.size(); tt++ )
	{
		_conesPru_prob[tt].setOffset();
	}

	// Now calculate the span
	double span = sqrt((double)hmax*hmax+(double)radius*radius);
	double dx_dou = hmax/span;
	double dy_dou = radius/span;
	_intSpan = computeAngleIndex(dx_dou, dy_dou);
	_voteDirec_prob = vector< pair< pair< int,int >, vector <int> > > (ntheta);
	// This copy is for a test on the voting, just to use the last state of the voting directions (updated gradients) and perform actually a vote in the biggest span
	_voteDirec_prob_copy = vector< pair< pair< int,int >, vector <int> > > (ntheta);
	for( int tt=0; tt<ntheta; tt++ ) // 0 256
	{
		pair<int,int> spann;
		if(tt-_intSpan<0)
		{
			spann.first = (tt-_intSpan+ntheta)%ntheta;
		}
		else
		{
			spann.first = (tt-_intSpan)%ntheta;
		}
		spann.second = (tt+_intSpan)%ntheta;

		// calculate how many points at disance 
		vector< int > vi_temp_1 (hmax-hmin+1,0); //hmax-hmin+1 integers of value 0
		for( int yy=spann.first; yy!=spann.second; yy=(yy+1)%ntheta )
		{
			for( int uu=0; uu<hmax-hmin+1; ++uu )
			{
				vi_temp_1[uu]  += _conesPru_prob[yy][uu].size();
			}
		}

		pair< pair< int,int >, vector<int> > spann_prob;
		spann_prob.first = spann;
		spann_prob.second = vi_temp_1;
		_voteDirec_prob[tt] = spann_prob;
	}

	// Just to test yosef method at the end when we have the best directions
	_voteDirec_prob_copy = _voteDirec_prob;


	//for( int tt=0; tt<ntheta; tt++ ) // 0 256
	//{
	//	cout<<endl<<"Tam: "<<tt<<", ";
	//	for( int uu=0; uu<_voteDirec_prob[tt].second.size(); ++uu )
	//	{
	//		cout<<_voteDirec_prob[tt].second[uu]<<" ";
	//	}
	//}

	cout<<"Cones computed";
	//cin>>stop11;

}


// ############################################################################################################################################################################
int ftkVoting::nextConeRadius(int h, int r)
{
	//  return  r-1;
	if (h<=0 || r<0) {
		cerr << "Error in ftkVoting::nextConeRadius: null cone\n";
		return 0;
	}
	return (int)((double) h * tan(atan((double)r/(double)h) * 0.7 )); // reduce the angle by 30% each time
}

//// ############################################################################################################################################################################
void ftkVoting::vote()
{
	_votingSumVotes = VotingDirType::New();
	VotingDirType::IndexType start;
	start[0] = 0;
	start[1] = 0;
	VotingDirType::SizeType size;
	size[0] = nx;
	size[1] = ny;
	VotingDirType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	_votingSumVotes->SetRegions( region );
	_votingSumVotes->Allocate();
	const VotingDirType::PixelType ceros = 0;
	_votingSumVotes->FillBuffer( ceros );
	_votingSumVotes->Update();

	_votingMaskVotes = VotingDirType::New(); 
	_votingMaskVotes->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	_votingMaskVotes->Allocate();
	_votingMaskVotes->FillBuffer( ceros );
	_votingMaskVotes->Update();


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






	//vector < vector < int > > votingMaskVotes_dir = vector< vector< int > >(nx*ny);
	//for( int yr=0; yr<nx*ny; yr++ )
	//{
	//	vector< int > bindedir = vector< int >(256,0); //256 int of value 0
	//	votingMaskVotes_dir.push_back( bindedir );
	//}

	//
	VotingDirType::Pointer imageOfVotingPixels = VotingDirType::New(); 
	imageOfVotingPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageOfVotingPixels->Allocate();
	imageOfVotingPixels->FillBuffer( ceros );
	imageOfVotingPixels->Update();
	//

	//
	VotingDirType::Pointer imageOfConexPixels = VotingDirType::New(); 
	imageOfConexPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageOfConexPixels->Allocate();
	imageOfConexPixels->FillBuffer( ceros );
	imageOfConexPixels->Update();
	//

	VotingDirType::Pointer imageGradXPixels = VotingDirType::New(); 
	imageGradXPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageGradXPixels->Allocate();
	imageGradXPixels->FillBuffer( ceros );
	imageGradXPixels->Update();
	//

	VotingDirType::Pointer imageGradYPixels = VotingDirType::New(); 
	imageGradYPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageGradYPixels->Allocate();
	imageGradYPixels->FillBuffer( ceros );
	imageGradYPixels->Update();
	//

	VotingDirType::Pointer imageMaxPixels = VotingDirType::New(); 
	imageMaxPixels->SetRegions( region );// Puede generar problemas depronto toque declara otro size y otro start
	imageMaxPixels->Allocate();
	imageMaxPixels->FillBuffer( ceros );
	imageMaxPixels->Update();
	//

	VotingDirType::PixelType * votingSumArray = _votingSumVotes->GetBufferPointer();
	VotingDirType::PixelType * votingMaskArray = _votingMaskVotes->GetBufferPointer();
	VotingDirType::PixelType * imageOfVotingPixelsArray = imageOfVotingPixels->GetBufferPointer();
	VotingDirType::PixelType * imageOfConexPixelsArray = imageOfConexPixels->GetBufferPointer();
	VotingDirType::PixelType * imageGradXPixelsArray = imageGradXPixels->GetBufferPointer();
	VotingDirType::PixelType * imageGradYPixelsArray = imageGradYPixels->GetBufferPointer();
	VotingDirType::PixelType * imageMaxPixelsArray = imageMaxPixels->GetBufferPointer();

	// Update the voting points to the new size (paddng)
	vector<VPoint2D>::iterator voting_points_begin = _voting_points.begin();
	vector<VPoint2D>::iterator voting_points_end = _voting_points.end();

	for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
		it->pos = it->x+nx*it->y; // Poscion donde esta x,y en nuestra recien creada _sum (Image)
	}

	//cout<<endl<<"Vote: Conputecone: "<<"Hmin: "<<_hmin<<", Hmax: "<<_hmax;
	computeCones(_hmin, _hmax, _radius);
	computeCones_prob(_hmin, _hmax, _radius); // Just for one moment

	// Primer Voto
	int nic=0;
	int count = 0;
	cout<<endl<<"Span"<<_intSpan;
	//	for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) { //it voting points
	//		//if ( count <256 ){
	//
	//		//cout<<endl<<"\tVoting pixel: "<<nic++<<", "<<it->pos;
	////		_cones[it->angIndex].vote(votingSumArray+it->pos, *it);
	//		
	//			votar(votingSumArray+it->pos,*it,it->angIndex);
	//
	//
	//
	//		////
	//		//imageOfVotingPixelsArray[it->pos] = 1;
	//		////(*it).x = 1;
	//		////(*it).y = 0;
	//		////_cones[count].vote(imageOfConexPixelsArray+it->pos, *it);
	//		////cout<<endl<<"Epezando a votar "<<count;
	//		//votar(imageOfConexPixelsArray+it->pos,*it,it->angIndex);
	//		////_conesPru[count].vote(imageOfConexPixelsArray+it->pos, *it,count);
	//
	//
	//		//stringstream out;
	//		//out<<count;
	//		//string s = out.str();
	//		//string filename = "output/cones/out_ImageOfVotes_"+s+".jpg";
	//		//if( nftkVot::writeImage< VotingDirType, OutputImageType >(imageOfConexPixels, filename.c_str() )){
	//		//	cout<<endl<<"\tProblema escribiendo";
	//		//}
	//		//memset(imageOfConexPixelsArray, 0,  npix*sizeof(VotingDirType::PixelType));
	//
	//		count++;
	//		//
	//		//}
	//	}


	int ytt;
	////// -------------------------------------------------------------------------------------------
	////// 1 of 2 || Step by step
	//for( int gg=0; gg<=_hmax-_hmin; gg++ )
	//{
	//	//cout<<endl<<"Q pasa "<<gg;
	//	for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//	{
	//		////cout<<endl<<"\t	Q pasa "<<gg;
	//		votar(votingSumArray+it->pos,*it,it->angIndex,gg);
	//		votar(votingMaskArray+it->pos,*it,it->angIndex,gg);
	//		////votar_dir(votingMaskVotes_dir,votingMaskArray+it->pos,*it,it->angIndex,gg,it->pos);

	//		//votar_prob(votingSumArray+it->pos,*it,it->angIndex,gg);
	//		//votar_prob(votingMaskArray+it->pos,*it,it->angIndex,gg);

	//	}
	//	stringstream out4;
	//	out4<<_intSpan;
	//	string s4 = out4.str();
	//	stringstream out5;
	//	out5<<_hmax-_hmin-gg;
	//	string s5 = out5.str();
	//	////string filename3 = "output/out_ImageOfMask_"+s3+".jpg";
	//	////if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
	//	////	cout<<endl<<"\tProblema escribiendo";
	//	////}

	//	string filename5 = "output/Sum_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//	if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingSumVotes, filename5.c_str() )){
	//		cout<<endl<<"\tProblema escribiendo";
	//	}
	//	string filename6 = "output/Mask_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
	//	if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
	//		cout<<endl<<"\tProblema escribiendo";
	//	}

	//	//string filename7 = "output/Sum_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//	//if( nftkVot::writeImage_mhdDouble< VotingDirType >(_votingSumVotes, filename7.c_str() )){
	//	//	cout<<endl<<"\tProblema escribiendo";
	//	//}
	//	//string filename8 = "output/Mask_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
	//	//if( nftkVot::writeImage_mhdDouble< VotingDirType >(_votingMaskVotes, filename8.c_str() )){
	//	//	cout<<endl<<"\tProblema escribiendo";
	//	//}

	//	// Get the maximum value of step by step
	//	for( int re=0; re<npix; ++re )
	//	{
	//		if( votingMaskArray[re] > imageMaxPixelsArray[re] )
	//		{
	//			imageMaxPixelsArray[re] = votingMaskArray[re];
	//		}
	//	}


	//	string filename9_max = "output/Max_stepbystep/out_ImageOfMaxs_"+s4+"__"+s5+".jpg";
	//	if( nftkVot::writeImage< VotingDirType, OutputImageType >(imageMaxPixels, filename9_max.c_str() )){
	//		cout<<endl<<"\tProblema escribiendo";
	//	}




	//	// Errase the mask (only to see step by step/by step) mas detallado todavia, ver los picos mas que todo
	//	memset(votingMaskArray, 0,  npix*sizeof(VotingDirType::PixelType));

	//	//// The amount of memory required turns out to be to high :( what to do ? 
	//	//string filename9 = "output/Direction_perpixel/out_DirPerPixel_"+s4+"__"+s5+".mhd";
	//	//if( nftkVot::writeImage_mhdDouble< VotingDirPerType >(votingDirPerPixel, filename9.c_str() )){
	//	//	cout<<endl<<"\tProblema escribiendo";
	//	//}

	//	//string filename9 = "output/Direction_perpixel/out_DirPerPixel_"+s4+"__"+s5+".mhd";
	//	//writer->SetInput( votingDirPerPixel );
	//	//writer->SetFileName( filename9.c_str() );
	//	//writer->Update();

	//	//VotingDirType::IndexType acceso;
	//	//acceso[0] = 45;
	//	//acceso[1] = 345;
	//	//cout<<endl<<" S4: "<<s4<<", S5: "<<s5<<", x: 45, y: 345, "<<_votingSumVotes->GetPixel(acceso);
	//	//cout<<endl<<"W pasa";
	//	//cin>>ytt;
	//}
	// -------------------------------------------------------------------------------------------
	// 2 of 2 ||  Voting without step by step
	// VOTING PERFECTO
	for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	{
		votar(votingMaskArray+it->pos,*it,it->angIndex);
	}
	//memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType::PixelType)); // Mask = Sum;
	memcpy(votingSumArray, votingMaskArray, npix*sizeof(VotingDirType::PixelType)); // Sum = Mask;







	cout<<endl<<"Fin primer voto: "<<ftkCone2D::contador_1;
	int rrr;
	//cin>>rrr;

	//stringstream out44;
	//out44<<_intSpan;
	//string s44 = out44.str();
	//string filename_dir = "output/Directions/Dri_vote_"+s44+".txt";
	//ofstream dirFilOut;
	//dirFilOut.open (filename_dir.c_str());
	//for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	//{
	//	dirFilOut << it->x <<"\t"<<it->y<<"\t"<<it->angIndex<<"\t"<<it->mag<<"\n";
	//}

	//stringstream out2;
	//out2<<_intSpan;
	//string s2 = out2.str();
	//string filename2 = "output/out_ImageOfSums_"+s2+".jpg";
	//if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingSumVotes, filename2.c_str() )){
	//	cout<<endl<<"\tProblema escribiendo";
	//}

	////votingSumVotes
	////Store Derivative in X
	//const char* filename = "output/out_ImageOfVotes.jpg";
	//if( nftkVot::writeImage< VotingDirType, OutputImageType >(imageOfVotingPixels, filename )){
	//	cout<<endl<<"\tProblema escribiendo";
	//}

	//const char* filename2 = "output/out_vote_1.jpg";
	//if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingSumVotes, filename2 )){
	//	cout<<endl<<"\tProblema escribiendo";
	//}

	int r = _radius;

	while( _intSpan>=0 )
	{
		cout<<endl<<"Span"<<_intSpan;
		cout<<endl<<"Updating cones:";
		updateCones(); // Reduce _intSpan
		cout<<"\t done";


		// Missing update direction
		int pru = _voting_points.size();
		vector<VPoint2D>::iterator it=voting_points_begin;
		cout<<endl<<"Updating direction:";
#pragma omp parallel for num_threads(16)
		for( int ty=0; ty<pru; ty++ ){
			updateDirection(*(it+ty)); // The Update is with the sum image
		}
		cout<<"\t done";
		//stringstream out45;
		//out45<<_intSpan;
		//string s45 = out45.str();
		//string filename_dir2 = "output/Directions/Dri_vote_"+s45+".txt";
		//ofstream dirFilOut2;
		//dirFilOut2.open (filename_dir2.c_str());
		//for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
		//{
		//	dirFilOut2 << it->x <<"\t"<<it->y<<"\t"<<it->angIndex<<"\t"<<it->mag<<"\n";
		//}

		//for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
		////getMaxDirection(*it);    
		//updateDirection(*it);
		//}   

		//memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType::PixelType));
		//memset(votingSumArray, 0,  npix*sizeof(VotingDirType::PixelType));

		//VOTING PERFECTO
		//for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
		//{
		//	votar(votingSumArray+it->pos,*it,it->angIndex);
		//}
		//stringstream out3;
		//out3<<_intSpan;
		//string s3 = out3.str();
		////string filename3 = "output/out_ImageOfMask_"+s3+".jpg";
		////if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingMaskVotes, filename3.c_str() )){
		////	cout<<endl<<"\tProblema escribiendo";
		////}
		//string filename4 = "output/out_ImageOfSums_"+s3+".jpg";
		//if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingSumVotes, filename4.c_str() )){
		//	cout<<endl<<"\tProblema escribiendo";
		//}

		memset(votingMaskArray, 0,  npix*sizeof(VotingDirType::PixelType));
		cout<<endl<<"Voting:";

		//// -------------------------------------------------------------------------------------------
		//// 1 of 2 || Step by step
		//for( int gg=0; gg<=_hmax-_hmin; gg++ )
		//{
		//	//cout<<endl<<"Q pasa "<<gg;
		//	for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
		//	{
		//		//cout<<endl<<"\t	Q pasa "<<gg;
		//		//votar(votingSumArray+it->pos,*it,it->angIndex,gg);
		//		votar(votingMaskArray+it->pos,*it,it->angIndex,gg);

		//		//votar_prob(votingSumArray+it->pos,*it,it->angIndex,gg);
		//		//votar_prob(votingMaskArray+it->pos,*it,it->angIndex,gg);
		//	}
		//	stringstream out4;
		//	out4<<_intSpan;
		//	string s4 = out4.str();
		//	stringstream out5;
		//	out5<<_hmax-_hmin-gg;
		//	string s5 = out5.str();

		//	string filename5 = "output/Sum_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
		//	if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingSumVotes, filename5.c_str() )){
		//		cout<<endl<<"\tProblema escribiendo";
		//	}
		//	string filename6 = "output/Mask_stepbystep/out_ImageOfSums_"+s4+"__"+s5+".jpg";
		//	if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingMaskVotes, filename6.c_str() )){
		//		cout<<endl<<"\tProblema escribiendo";
		//	}

		//	//string filename7 = "output/Sum_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
		//	//if( nftkVot::writeImage_mhdDouble< VotingDirType >(_votingSumVotes, filename7.c_str() )){
		//	//	cout<<endl<<"\tProblema escribiendo";
		//	//}
		//	//string filename8 = "output/Mask_stepbystep_mhd/out_ImageOfSums_"+s4+"__"+s5+".mhd";
		//	//if( nftkVot::writeImage_mhdDouble< VotingDirType >(_votingMaskVotes, filename8.c_str() )){
		//	//	cout<<endl<<"\tProblema escribiendo";
		//	//}
		//}
		//// -------------------------------------------------------------------------------------------
		//// 2 of 2 ||  Voting without step by step
		//// VOTING PERFECTO
		//for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
		//{
		//	votar(votingMaskArray+it->pos,*it,it->angIndex);
		//}
		// ------------------------------------------------------------------
		// 3 of 2 ||  Voting without step by step
		// VOTING PERFECTO
		for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
		{
			//votar(votingMaskArray+it->pos,*it,it->angIndex);
			for( int tt=_voteDirec[it->angIndex].first; tt!=_voteDirec[it->angIndex].second; tt=(tt+1)%ntheta )
			{
				_conesPru[tt].vote(votingMaskArray+it->pos,*it);
			}
		}
		// ------------------------------------------------------------------
		cout<<"\t voting done";
		//memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType::PixelType)); // Mask = Sum;
		memcpy(votingSumArray, votingMaskArray, npix*sizeof(VotingDirType::PixelType)); // Mask = Sum;




		// TEST TO GENERATE THE PROCESS SLICE BY SLICE
		if( _intSpan == 23 )
		{
			stringstream outSlice;
			outSlice<<_intSpan;
			string sSlice = outSlice.str();
			string filenameSlicebySlice = "output/BySlice/out_"+sSlice+"__"+_prefix+".tif";
			if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingSumVotes, filenameSlicebySlice.c_str() )){
				cout<<endl<<"\tProblema escribiendo";
			}

	//int nxxX = _votingBySlice->GetLargestPossibleRegion().GetSize()[0];
	//int nyxX = _votingBySlice->GetLargestPossibleRegion().GetSize()[1];
	//int nzxX = _votingBySlice->GetLargestPossibleRegion().GetSize()[2];
	//int npixX = nxxX*nyxX*nzxX;
	//cout<<endl<<"PASDFASD: "<<npixX;

			// Copy one slice
			nftkVot::InputImageType_3D::PixelType * votingBySliceArray = _votingBySlice->GetBufferPointer();
			for( unsigned int yyyo = _slice*nx*ny; yyyo<(_slice+1)*nx*ny; ++yyyo )
			{
				votingBySliceArray[yyyo] = votingSumArray[yyyo-_slice*nx*ny];
			}
			if( (_slice == 20) || (_slice == 50) || (_slice == 100) || (_slice == 150) || (_slice == 200) || (_slice == 274) )
			{
			stringstream outSlice_3D;
			outSlice_3D<<_intSpan;
			string sSlice_3D = outSlice_3D.str();
			string filenameSlicebySlice_3D = "output/BySlice_3D/out_"+sSlice_3D+"__"+_prefix+".tif";
			if( nftkVot::writeImage_3D< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_16 >(_votingBySlice, filenameSlicebySlice_3D.c_str() )){
				cout<<endl<<"\tProblema escribiendo";
			}
			}
				
			// To end the iteration
			_intSpan=0;
		}



		if( _intSpan==0)
			break;

	}

	// FOR YOUSEF SEGMENTATION
	// Here we have good direction for the votes, now what we can do is actually "diffuse in that direction" so that we get a good result for yousef segmentation
	// Restore the initial span of the voting direction
	_voteDirec = _voteDirec_copy;
	// VOTE ONE TIME
	memset(votingMaskArray, 0,  npix*sizeof(VotingDirType::PixelType));
	for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
	{
		votar(votingMaskArray+it->pos,*it,it->angIndex);
	}
	string filenameSlicebySliceLastbigSapan = "output/BySliceLastbigSapan/out___"+_prefix+".tif";
	if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingMaskVotes, filenameSlicebySliceLastbigSapan.c_str() )){
		cout<<endl<<"\tProblema escribiendo";
	}

			// Copy one slice
			nftkVot::InputImageType_3D::PixelType * votingBySliceLastbigSapanArray = _votingBySliceLastbigSapan->GetBufferPointer();
			for( unsigned int yyyo = _slice*nx*ny; yyyo<(_slice+1)*nx*ny; ++yyyo )
			{
				votingBySliceLastbigSapanArray[yyyo] = votingMaskArray[yyyo-_slice*nx*ny];
			}
			if( (_slice == 50) || (_slice == 100) || (_slice == 150) || (_slice == 200) || (_slice == 274) )
			{
			string filenameSlicebySliceLastbigSapan_3D = "output/BySliceLastbigSapan_3D/out___"+_prefix+".tif";
			if( nftkVot::writeImage_3D< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_16 >(_votingBySliceLastbigSapan, filenameSlicebySliceLastbigSapan_3D.c_str() )){
				cout<<endl<<"\tProblema escribiendo";
			}
			}

	// Using probability
	_voteDirec_prob = _voteDirec_prob_copy;
	memset(votingMaskArray, 0,  npix*sizeof(VotingDirType::PixelType));
	for( int gg=0; gg<=_hmax-_hmin; gg++ )
	{
		for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++)
		{
			votar_prob(votingMaskArray+it->pos,*it,it->angIndex,gg);
		}
	}
	string filenameSlicebySliceLastbigSapanProb = "output/BySliceLastbigSapanProb/out___"+_prefix+".tif";
	if( nftkVot::writeImage< VotingDirType, OutputImageType >(_votingMaskVotes, filenameSlicebySliceLastbigSapanProb.c_str() )){
		cout<<endl<<"\tProblema escribiendo";
	}

			// Copy one slice
			nftkVot::InputImageType_3D::PixelType * votingBySliceLastbigSapanProbArray = _votingBySliceLastbigSapanProb->GetBufferPointer();
			for( unsigned int yyyo = _slice*nx*ny; yyyo<(_slice+1)*nx*ny; ++yyyo )
			{
				votingBySliceLastbigSapanProbArray[yyyo] = votingMaskArray[yyyo-_slice*nx*ny];
			}
			if( (_slice == 50) || (_slice == 100) || (_slice == 150) || (_slice == 200) || (_slice == 274) )
			{
			string filenameSlicebySliceLastbigSapanProb_3D = "output/BySliceLastbigSapanProb_3D/out___"+_prefix+".tif";
			if( nftkVot::writeImage_3D< nftkVot::InputImageType_3D, nftkVot::InputImageType_3D_16 >(_votingBySliceLastbigSapanProb, filenameSlicebySliceLastbigSapanProb_3D.c_str() )){
				cout<<endl<<"\tProblema escribiendo";
			}
			}







	/*memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType::PixelType));*/

	//while (r>= 0) 
	//{
	//	cout<<endl<<"	r: "<<r;

	//	computeCones(_hmin, _hmax, r); // Eror pues los acabamos de calcular

	//	for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) {	
	//			getMaxDirection(*it);    
	//	}    

	//	memcpy(votingMaskArray, votingSumArray, npix*sizeof(VotingDirType::PixelType));
	//	memset(votingSumArray, 0,  npix*sizeof(VotingDirType::PixelType));


	//	 //NO ENTIENDO PARA QUE ES ESTO
	//	for(vector<VPoint2D>::iterator it=voting_points_begin; it!=voting_points_end; it++) 
	//	{
	//			_cones[it->angIndex].vote(votingSumArray+it->pos, votingMaskArray+it->pos, *it);
	//	}

	//	if (r<=0)
	//		break;

	//	r = nextConeRadius(_hmax, r);

	//}
	////computeCones(_hmin, _hmax, 0);

	//cout<<endl<<"Vote: End";
}


//	_sum = doubleImage(nx, ny);
//	_mask = doubleImage(nx, ny);
//
//	double* pmask = _mask.ptr();
//	double* psum = _sum.ptr();
//
//	vector<VPoint2D>::iterator voting_points_begin = _voting_points.begin();
//	vector<VPoint2D>::iterator voting_points_end = _voting_points.end();
//
//	for(vector<VPoint2D>::iterator it=voting_points_begin; 
//		it!=voting_points_end; it++) {	
//			it->pos = _sum.pos(it->x, it->y); // Poscion donde esta x,y en nuestra recien creada _sum (Image)
//	}
//
//	computeCones(_hmin, _hmax, _radius);
//
//	for(vector<VPoint2D>::iterator it=voting_points_begin; 
//		it!=voting_points_end; it++) {
//			_cones[it->angIndex].vote(psum+it->pos, *it);
//	}
//
//	Image* roi;
//	int loop = 0;
//	int r = _radius;
//
//	while (r>= 0) {
//
//		//    cout << "Radius = " << r<< endl;
//		computeCones(_hmin, _hmax, r);
//		// _reg_type = BrightReg;
//
//		if (_save_history) {
//			roi = Roi::getSubImage(_sum, bw, bw, nx-bw-1, ny-bw-1);
//			string vote_image_number = toString(r);
//			if (r<10) vote_image_number = "0"+ vote_image_number;
//			/* comment for saving 16bit ics file instead
//			writePGM(_prefix+"vote_"+vote_image_number+".pgm",
//			convertToByteImage(*roi, true));
//			*/
//			image_io::writeICS(_prefix+"vote_"+vote_image_number,*roi);
//			delete roi;
//		}
//
//		for(vector<VPoint2D>::iterator it=voting_points_begin; 
//			it!=voting_points_end; it++) {	
//				getMaxDirection(*it);    
//		}    
//
//		memcpy(pmask, psum, npix*sizeof(double));
//		memset(psum, 0,  npix*sizeof(double));
//
//		for(vector<VPoint2D>::iterator it=voting_points_begin; 
//			it!=voting_points_end; it++) {
//				_cones[it->angIndex].vote(psum+it->pos, pmask+it->pos, *it);
//		}
//
//		if (r<=0)
//			break;
//
//		r = nextConeRadius(_hmax, r);
//
//	}
//	computeCones(_hmin, _hmax, 0);
//}
//
//// ############################################################################################################################################################################
void inline ftkVoting::votar(VotingDirType::PixelType * p, const VPoint2D& vp, int angl_indx)
{
	//if( angl_indx==1)
	//{
	//	cout<<endl<<_voteDirec[angl_indx].first<<" "<<_voteDirec[angl_indx].second<<", aca";
	//}
	//cout<<endl<<"\t angl_indx "<<angl_indx<<","<<_voteDirec[angl_indx].first<<","<<_voteDirec[angl_indx].second;
	for( int tt=_voteDirec[angl_indx].first; tt!=_voteDirec[angl_indx].second; tt=(tt+1)%ntheta )
	{
		_conesPru[tt].vote(p,vp);
	}
}

//// ############################################################################################################################################################################
void inline ftkVoting::votar(VotingDirType::PixelType * p, const VPoint2D& vp, int angl_indx, int dist)
{
	//if( angl_indx==1)
	//{
	//	cout<<endl<<_voteDirec[angl_indx].first<<" "<<_voteDirec[angl_indx].second<<", aca";
	//}
	//cout<<endl<<"\t angl_indx "<<angl_indx<<","<<_voteDirec[angl_indx].first<<","<<_voteDirec[angl_indx].second;
	for( int tt=_voteDirec[angl_indx].first; tt!=_voteDirec[angl_indx].second; tt=(tt+1)%ntheta )
	{
		_conesPru[tt].vote(p,vp,dist);
	}
}

////// ############################################################################################################################################################################
//void inline ftkVoting::votar_prob(VotingDirType::PixelType * p, const VPoint2D& vp, int angl_indx)
//{
//	//if( angl_indx==1)
//	//{
//	//	cout<<endl<<_voteDirec[angl_indx].first<<" "<<_voteDirec[angl_indx].second<<", aca";
//	//}
//	//cout<<endl<<"\t angl_indx "<<angl_indx<<","<<_voteDirec[angl_indx].first<<","<<_voteDirec[angl_indx].second;
//	for( int tt=_voteDirec_prob[angl_indx].first.first; tt!=_voteDirec_prob[angl_indx].first.second; tt=(tt+1)%ntheta )
//	{
//		_conesPru_prob[tt].vote(p,vp);
//	}
//}

//// ############################################################################################################################################################################
void inline ftkVoting::votar_prob(VotingDirType::PixelType * p, const VPoint2D& vp, int angl_indx, int dist)
{
	//if( angl_indx==1)
	//{
	//	cout<<endl<<_voteDirec[angl_indx].first<<" "<<_voteDirec[angl_indx].second<<", aca";
	//}
	//cout<<endl<<"\t angl_indx "<<angl_indx<<","<<_voteDirec[angl_indx].first<<","<<_voteDirec[angl_indx].second;
	for( int tt=_voteDirec_prob[angl_indx].first.first; tt!=_voteDirec_prob[angl_indx].first.second; tt=(tt+1)%ntheta )
	{
		_conesPru_prob[tt].vote(p,vp,dist,_voteDirec_prob[angl_indx].second[dist]);
	}
}

//// ############################################################################################################################################################################
void inline ftkVoting::votar_dir(vector<vector<int> >& p_dir,VotingDirType::PixelType * p, const VPoint2D& vp, int angl_indx, int dist, int& offset_1)
{
	//if( angl_indx==1)
	//{
	//	cout<<endl<<_voteDirec[angl_indx].first<<" "<<_voteDirec[angl_indx].second<<", aca";
	//}
	//cout<<endl<<"\t angl_indx "<<angl_indx<<","<<_voteDirec[angl_indx].first<<","<<_voteDirec[angl_indx].second;
	for( int tt=_voteDirec[angl_indx].first; tt!=_voteDirec[angl_indx].second; tt=(tt+1)%ntheta )
	{
		_conesPru[tt].vote_dir(p_dir,p,vp,dist,offset_1);
	}
}

// ############################################################################################################################################################################
void ftkVoting::updateCones()
{
	if(_intSpan==0)
	{
		int rr;
		cout<<endl<<"Error en updatecones";
		cin>>rr;
	}

	int reduction = 2;
	while(_intSpan-reduction<0) //Si el span se reduce a menos de 0, entonces reduzca hasta que quede en cero
	{
		reduction--;
	}
	_intSpan = _intSpan-reduction;
	int temp;
	//_voteDirec = vector< pair< int,int > > (ntheta);
	for( int tt=0; tt<ntheta; tt++ )
	{
		_voteDirec[tt].first = (_voteDirec[tt].first+reduction)%ntheta;
		temp = _voteDirec[tt].second-reduction;
		if(temp<0)
		{
			_voteDirec[tt].second = (temp+ntheta)%ntheta;
		}
		else
		{
			_voteDirec[tt].second = (temp)%ntheta;
		}
	}

}